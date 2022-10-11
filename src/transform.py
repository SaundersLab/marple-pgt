import tempfile
from collections import defaultdict
from os import makedirs
from os.path import abspath, basename, join
from re import finditer, sub
from typing import Dict, Iterable, List, Optional, Tuple
import subprocess
import matplotlib.pyplot as plt
import pandas as pd
from Bio import Phylo
from Bio.SeqIO import parse


from utils import (darken_color, file, get_sample_name_and_extenstion, pushd,
                   run, string_to_color, write_fasta)


# Make a pileup and return the path to it
def reads_to_pileup(
    fastq: str,
    reference: str,
    out_dir: str,
    threads=1,
    trim=True,
) -> str:
    makedirs(out_dir, exist_ok=True)
    sample_name, sample_ext = get_sample_name_and_extenstion(fastq, 'fastq')
    threads = str(threads)
    if trim:
        trimmed = join(out_dir, f'{sample_name}_porechopped{sample_ext}')
        print('trimming', end=' ', flush=True)
        run(['porechop', '--threads', threads, '-i', fastq], trimmed)
    else:
        trimmed = fastq
    print('aligning', end=' ', flush=True)
    run(['bwa', 'index', reference])
    aligned_sam = join(out_dir, f'{sample_name}.sam')
    run(['bwa', 'mem', '-t', threads, reference, trimmed], aligned_sam)
    aligned_bam = join(out_dir, f'{sample_name}_unsorted.bam')
    run(['samtools', 'view', '-@', threads, '-S', '-b', aligned_sam], aligned_bam)
    sorted_bam = join(out_dir, f'{sample_name}.bam')
    run(['samtools', 'sort', '-@', threads, aligned_bam], sorted_bam)
    run(['samtools', 'faidx', reference])
    pileup = join(out_dir, f'{sample_name}.pileup')
    run(['samtools', 'mpileup', '-f', reference, sorted_bam], pileup)

    # Let the caller know where to find the pileup file
    return pileup

def snp_ratios_and_snp_freq_to_consensus(
    snp_ratios_path: str,
    snp_freq_path: str,
    reference_path: str,
    out_dir: str,
    min_match_depth: int,
):
    sample_name, _ = get_sample_name_and_extenstion(snp_ratios_path, 'snp_ratios')
    
    snp_freq = pd.read_csv(snp_freq_path, sep='\t', header=None)
    snp_freq.columns = [
        'gene', 'pos', 'pos1', 'ref', 'depth', 'allele', 'genotype', 'ratios'
    ]
    # Drop the position where 3+ bases qualified (e.g. T ACT ? 0.2,0.6,0.2)
    snp_freq = snp_freq[snp_freq.genotype != '?']

    snp_ratios = pd.read_csv(snp_ratios_path, sep='\t', header=None)
    snp_ratios.columns = ['gene', 'pos', 'ref', 'depth', 'alts', 'ratios']

    allele_to_code = {
        'AA': 'A',  'CC': 'C',  'GG': 'G',  'TT': 'T',  'UU': 'U',
        'AT': 'W',  'CG': 'S',  'AC': 'M',  'GT': 'K',  'AG': 'R',  'CT': 'Y',
        'TA': 'W',  'GC': 'S',  'CA': 'M',  'TG': 'K',  'GA': 'R',  'TC': 'Y',
    }

    # Start with an empty consensus
    consensus = {str(r.id): ['N'] * len(r.seq) for r in parse(reference_path, 'fasta')}

    # Fill in the consensus using alternate bases where coverage >= min_depth
    for gene, pos, _, allele in snp_freq[['gene', 'pos', 'genotype', 'allele']].values:
        consensus[gene][pos] = allele_to_code.get(allele, consensus[gene][pos])

    # Fill in the consensus using reference matches where coverage >= min_match_depth
    ref_matches = snp_ratios[snp_ratios.alts == snp_ratios.ref]
    ref_matches_enough_depth = ref_matches.query(f'depth >= {min_match_depth}')
    for gene, pos, ref in ref_matches_enough_depth[['gene', 'pos', 'ref']].values:
        # snp_ratios pos comes directly from pileup which is 1-based index
        consensus[gene][pos - 1] = ref

    # Convert the sequence lists into strings
    consensus_fasta = {gene: ''.join(l) for gene, l in consensus.items()}

    # Write out the consensus
    consensus_path = join(out_dir, f'{sample_name}.fasta')
    write_fasta(consensus_fasta, consensus_path)

    # Let the caller know where to find the consensus file
    return consensus_path


# Make a consensus and return the path to it


def pileup_to_consensus(
    pileup_path: str,
    reference_path: str,
    out_dir: str,
    min_snp_depth: int,
    min_match_depth: int,
    hetero_min=.25,
    hetero_max=.75,
) -> str:

    sample_name, _ = get_sample_name_and_extenstion(pileup_path, 'pileup')

    snp_ratios_path = join(out_dir, f'{sample_name}_snp_ratios.tsv')
    pileup_to_snp_ratios(pileup_path, snp_ratios_path)

    snp_freq_path = join(out_dir, f'{sample_name}_snp_freq.tsv')
    snp_ratios_to_snp_freq(
        snp_ratios_path, snp_freq_path, min_snp_depth, hetero_min, hetero_max
    )

    # Let the caller know where to find the consensus file
    return snp_ratios_and_snp_freq_to_consensus(
        snp_ratios_path=snp_ratios_path,
        snp_freq_path=snp_freq_path,
        reference_path=reference_path,
        out_dir=out_dir,
        min_match_depth=min_match_depth,
    )

def _subtract_indels_from_counts(upper_reads: str, base_counts: Dict[str, int]) -> None:
    for indel in finditer(r'([+|-])(\d+)(\w+)', upper_reads):
        indel_size = int(indel.group(2))
        indel_sequence = indel.group(3)[:indel_size]
        for base in base_counts:
            base_counts[base] -= indel_sequence.count(base)


def _base_ratios_from_reads(reads: str, depth: int, ref: str) -> Dict[str, float]:
    upper_reads = reads.upper()

    # ^G. marks the start of a read segment where the ASCII minus 33 of G is the mapping quality
    # and the character AFTER those 2 is the actual read - in this case a match to the reference.
    # Drop the first 2 characters as they are not actual reads
    if '^' in upper_reads:
        upper_reads = sub('\^.', '', upper_reads)

    base_counts = {base: upper_reads.count(base) for base in 'ACTG'}
    if ref in 'ACTG':
        base_counts[ref] += upper_reads.count('.') + upper_reads.count(',')
    if '+' in upper_reads or '-' in upper_reads:
        _subtract_indels_from_counts(upper_reads, base_counts)
    # Only calculate ratios for present bases
    base_ratios = {
        base: round(base_counts[base] / depth, 3)
        for base in base_counts if base_counts[base]
    }
    return base_ratios


def _pileup_row_to_snp_ratio_row(row: str) -> str:
    contig, pos_str, ref, depth_str, reads, *_ = row.split('\t')
    depth = int(depth_str)
    base_to_ratio = _base_ratios_from_reads(reads, depth, ref)
    # sort by base for backward compatability
    base_ratios = sorted(base_to_ratio.items(),
                         key=lambda base_and_ratio: base_and_ratio[0])
    bases_str = ','.join(base for base, _ in base_ratios)
    ratios_str = ','.join(str(ratio) for _, ratio in base_ratios)
    return '\t'.join((contig, pos_str, ref, depth_str, bases_str, ratios_str))


def pileup_to_snp_ratios(pileup_path: str, snp_ratios_path: str) -> None:
    with file(pileup_path, 'rt') as f_in, file(snp_ratios_path, 'wt') as f_out:
        for row in f_in:
            f_out.write(_pileup_row_to_snp_ratio_row(row) + '\n')

# Determine the genotype at each position and filter by read coverage


def _get_genotype_and_valid_bases_and_valid_ratios(
    ref: str, bases: List[str], ratio_strs: list, hetero_min: float, hetero_max: float
) -> Optional[Tuple[str, str, str]]:

    """Use the heterozygosity thresholds to filter the base ratios and determine the genotype:
    Genotype	Condition	                                        Example
    0/0	        pref ≥ hetero_max or only refbase qualified	        G GG 0/0 1
    0/1	        ref base and 1 alt base qualified	                T CT 0/1 0.217,0.783
    1/1	        palt 1 ≥ hetero_max or only 1 alt base qualified	T CC 1/1 0.822
    1/2	        2 alt bases qualified and ref base did not	        C AT 1/2 0.476,0.524
    ?	        3+ bases qualified	                                T ACT ? 0.2,0.6,0.2
    """

    ratios = map(float, ratio_strs)

    genotype = ''
    valid_bases = ''
    valid_ratios = ''

    for (base, ratio, ratio_str) in zip(bases, ratios, ratio_strs):
        if ratio < hetero_min:
            continue
        valid_ratios += ratio_str + ','
        if ratio >= hetero_max:
            genotype = '0/0' if base == ref else '1/1'
            valid_bases = base + base
            break
        valid_bases += base

    if not valid_bases:
        return None

    if not genotype:
        if len(valid_bases) == 1:
            if valid_bases == ref:  # If only ref has min <= freq <= max
                return None
            genotype = '1/1'
            valid_bases = valid_bases + valid_bases
        elif len(valid_bases) == 2:
            genotype = '0/1' if ref in valid_bases else '1/2'
        else:
            genotype = '?'

    valid_ratios = valid_ratios.rstrip(',')
    return genotype, valid_bases, valid_ratios


def _snp_ratio_row_to_snp_freq_row(
    row: str, min_depth: int, hetero_min: float, hetero_max: float
) -> Optional[str]:

    contig, pos, ref, depth, bases_str, ratios = row.split('\t')

    if not bases_str or int(depth) < min_depth:
        return None

    bases = bases_str.split(',')
    ratio_strs = ratios.rstrip().split(',')

    genotype_and_valid_bases_and_valid_ratios = _get_genotype_and_valid_bases_and_valid_ratios(
        ref, bases, ratio_strs, hetero_min, hetero_max
    )

    if not genotype_and_valid_bases_and_valid_ratios:
        return None

    genotype, valid_bases, valid_ratios = genotype_and_valid_bases_and_valid_ratios

    return '\t'.join((contig, str(int(pos) - 1), pos, ref, depth, valid_bases, genotype, valid_ratios))


def snp_ratios_to_snp_freq(
    snp_ratios_path: str,
    snp_freq_path: str,
    min_depth: int,
    hetero_min: float,
    hetero_max: float,
) -> None:
    with file(snp_ratios_path, 'rt') as f_in, file(snp_freq_path, 'wt') as f_out:
        for row in f_in:
            snp_freq_row = _snp_ratio_row_to_snp_freq_row(
                row, min_depth, hetero_min, hetero_max
            )
            if snp_freq_row:
                f_out.write(snp_freq_row + '\n')


def reads_to_consensus(
    fastq: str,
    reference: str,
    out_dir: str,
    min_snp_depth: int = 20,
    min_match_depth: int = 2,
    hetero_min: float = .25,
    hetero_max: float = .75,
    threads=1,
    trim=True,
) -> str:
    pileup = reads_to_pileup(fastq, reference, out_dir, threads=threads, trim=trim)
    consensus = pileup_to_consensus(
        pileup, reference, out_dir, min_snp_depth,
        min_match_depth, hetero_min, hetero_max
    )
    return consensus


def reads_to_fastqc(fastq: str, out_dir: str) -> Tuple[str, str]:
    makedirs(out_dir, exist_ok=True)
    sample_name, _ = get_sample_name_and_extenstion(fastq, 'fastq')
    fastqc_page = join(out_dir, f'{sample_name}_fastqc.html')
    fastqc_data = join(out_dir, f'{sample_name}_fastqc.zip')
    run(['fastqc', '--quiet', '-o', out_dir, fastq])
    return fastqc_page, fastqc_data


def alignment_to_flagstat(alignment: str, out_dir: str) -> str:
    makedirs(out_dir, exist_ok=True)
    sample_name, _ = get_sample_name_and_extenstion(alignment, 'alignment')
    flagstat = join(out_dir, f'{sample_name}.txt')
    run(['samtools', 'flagstat', alignment], flagstat)
    return flagstat

def seq_to_pct_coverage(seq: str) -> float:
    seq = str(seq).upper()
    n_unknown = seq.count('?') + seq.count('N')
    pct_unknown = 100 * n_unknown / len(seq)
    return 100 - pct_unknown

def consensus_to_coverage(consensus, out_dir, step=1):
    makedirs(out_dir, exist_ok=True)
    sample_name, _ = get_sample_name_and_extenstion(consensus, 'fasta')
    out_path = join(out_dir, f'{sample_name}.csv')
    cov_pcts = [seq_to_pct_coverage(r.seq) for r in parse(consensus, 'fasta')]

    pct_coverage_thresholds = list(range(0, 101, step))
    number_of_genes_with_pct_coverage_ge_thresholds = [
        len([pct for pct in cov_pcts if pct >= min_pct])
        for min_pct in pct_coverage_thresholds
    ]
    pd.DataFrame({
        'pct_exon_positions_covered': pct_coverage_thresholds,
        'number_of_genes': number_of_genes_with_pct_coverage_ge_thresholds,
    }).to_csv(out_path, index=None, header=None)
    return out_path

# This is a work around to ensure that multiqc can create the plot specified in
# config/multiqc_config.yaml.
# Unless there is at least one _mqc.yaml in the directory, multiqc will not 
# pick up any of the custom_data specified in the --config file. To ensure it
# does pick up those files, we create a file for each sample with only id.
def create_empty_config_required_for_gene_coverage_mqc(sample: str, out_dir: str):
    config_name = f'{sample}_empty_config_required_for_gene_coverage_mqc'
    with open(join(out_dir, f'{config_name}.yaml'), 'w') as f:
        f.write(f'id: "{config_name}"' + '\n')

def sample_report(sample_dir: str):
    sample_name = basename(sample_dir)
    out_dir = join(sample_dir, 'report')
    # TODO: make this more flexible
    reads_to_fastqc(join(sample_dir, f'{sample_name}.fastq'), out_dir)
    alignment_to_flagstat(join(sample_dir, f'{sample_name}.bam'), out_dir)
    consensus_to_coverage(join(sample_dir, f'{sample_name}.fasta'), out_dir)
    create_empty_config_required_for_gene_coverage_mqc(sample_name, out_dir)


def reads_list_to_consensus_with_report(
    fastq_paths: List[str],
    reference: str,
    out_dirs: List[str],
    multiqc_config: str,
    threads=1,
    trim=True,
):
    for fastq_index, (fastq, out_dir) in enumerate(zip(fastq_paths, out_dirs)):
        sample_name = get_sample_name_and_extenstion(fastq, 'fastq')[0]
        print(f'{fastq_index + 1}/{len(fastq_paths)} {sample_name}:', end=' ', flush=True)
        reads_to_consensus(
            fastq=fastq,
            reference=reference,
            out_dir=out_dir,
            threads=threads,
            trim=trim,
        )
        print('assessing', flush=True)
        sample_report(out_dir)
    print('Report: compiling')
    report_dirs = [join(out_dir, 'report') for out_dir in out_dirs]
    run(['multiqc', '--config', multiqc_config] + report_dirs, out='/dev/null')


