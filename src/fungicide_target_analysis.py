from Bio.Data.IUPACData import ambiguous_dna_values
from itertools import product
from typing import Dict, List, Tuple
from Bio.Seq import translate, Seq
import pandas as pd
from Bio.SeqIO import parse
from utils import get_sample_name_and_extenstion

def disambiguate(code: str):
    """
    Examples
    --------
    >>> disambiguate('W')
    'AT'
    """
    return ambiguous_dna_values[code]

def disambiguated_translate(ambiguous_codon: str) -> str:
    """
    Examples
    --------
    >>> disambiguated_translate('WGW')
    '*CRS'
    """
    assert len(ambiguous_codon) == 3
    return ''.join(set(
        translate(''.join(unambiguous_codon))
        for unambiguous_codon in product(*map(disambiguate, ambiguous_codon))
    ))

def fill_N_with_ref(ref: str, alt: str) -> str:
    """
    Examples
    --------
    >>> fill_N_with_ref('ATG', 'ATN')
    ATG'
    """
    return ''.join(r if a == 'N' else a for r, a in zip(ref, alt))

def codon_diff(ref_codon: str, alt_codon: str) -> Dict[str, str]:
    """
    Examples
    --------
    >>> codon_diff('ATG', 'WTG')
    {'alt_codon_filled': 'WTG', 'ref_acid': 'M',
    'alt_possible_acids': 'ML', 'synonymous': 'maybe'}
    """
    alt_codon_filled = fill_N_with_ref(ref_codon, alt_codon)
    ref_acid = str(translate(ref_codon))
    alt_possible_acids = disambiguated_translate(alt_codon_filled)

    if ref_acid == alt_possible_acids:
        synonymous = 'yes'
    elif ref_acid in alt_possible_acids:
        synonymous = 'maybe'
    else:
        synonymous = 'no'

    return {
        'alt_codon_filled': alt_codon_filled,
        'ref_acid': ref_acid,
        'alt_possible_acids': alt_possible_acids,
        'synonymous': synonymous,
    }

def ranges_to_indices(ranges: List[Tuple[int]]) -> List[int]:
    """
    Examples
    --------
    >>> ranges_to_indices([(1, 4), (6, 9)])
    [1, 2, 3, 6, 7, 8]
    """
    indices = set()
    for s, e in ranges:
        indices.update(set(range(s, e)))
    return sorted(indices)

def gene_seqs_to_cds_diffs(
    ref_gene_seq: str,
    alt_gene_seq: str,
    cds_ranges: List[Tuple[int]],
) -> List[dict]:
    """
    Examples
    --------
    >>> gene_seqs_to_cds_diffs(
            ref_gene_seq = 'AATGGGCATTAAG',
            alt_gene_seq = 'GWTGGGCAATNAT',
            cds_ranges = [(1, 4), (6, 12)],
        )
    [{'alt_codon_filled': 'WTG', 'ref_acid': 'M', 'alt_possible_acids': 'ML',
    'synonymous': 'maybe', 'ref_codon': 'ATG', 'alt_codon': 'WTG', 'codon_index': 0,
    'indices': '1,2,3'},
    {'alt_codon_filled': 'CAA', 'ref_acid': 'H', 'alt_possible_acids': 'Q',
    'synonymous': 'no', 'ref_codon': 'CAT', 'alt_codon': 'CAA', 'codon_index': 1,
    'indices': '6,7,8'}]
    """
    diffs = []
    cds_indices = ranges_to_indices(cds_ranges)
    assert len(cds_indices) % 3 == 0
    for codon_index in range(len(cds_indices) // 3):
        codon_indices = cds_indices[codon_index * 3: (codon_index + 1) * 3]
        ref_codon = ''.join(ref_gene_seq[i] for i in codon_indices)
        alt_codon = ''.join(alt_gene_seq[i] for i in codon_indices)
        if alt_codon == ref_codon:
            continue
        diff = codon_diff(ref_codon, alt_codon)
        if diff['alt_codon_filled'] == ref_codon:
            continue
        diff['ref_codon'] = ref_codon
        diff['alt_codon'] = alt_codon
        diff['codon_index'] = codon_index
        diff['indices'] = ','.join(map(str, codon_indices))
        diffs.append(diff)
    return diffs

def gene_seqs_to_non_cds_diffs(
    ref_gene_seq: str,
    alt_gene_seq: str,
    cds_ranges: List[Tuple[int]],
) -> List[dict]:
    """
    Examples
    --------
    >>> gene_seqs_to_non_cds_diffs(
        ref_gene_seq = 'AATGGGCATTAAG',
        alt_gene_seq = 'GWTGGGCAATNAT',
        cds_ranges = [(1, 4), (6, 12)],
    )
    [{'pos': 0, 'ref': 'A', 'alt': 'G'},
    {'pos': 12, 'ref': 'G', 'alt': 'T'}]
    """
    diffs = []
    cds_indices = ranges_to_indices(cds_ranges)
    for i, (ref, alt) in enumerate(zip(ref_gene_seq, alt_gene_seq)):
        if (i in cds_indices) or (alt == 'N'):
            continue
        if ref != alt:
            diffs.append({'pos': i, 'ref': ref, 'alt': alt})
    return diffs

def fungicide_target_analysis(
    reference: str,
    gff_path: str,
    consensus: str,
    cds_diffs_out: str,
    non_cds_diffs_out: str,
    snp_summary_out: str,
):
    
    sample_name, _ = get_sample_name_and_extenstion(consensus, 'fasta')
    
    gff = pd.read_csv(gff_path, sep='\t', header=None)
    gff.columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    
    assert set(gff.strand) == {'+'}
    assert gff[gff['type'] == 'mRNA'].groupby('seqid').size().max() == 1

    gene_to_cds_ranges = {
        gene: list(zip(cds_features.start - 1, cds_features.end))
        for gene, cds_features in gff.query('type == "CDS"').groupby('seqid')
    }

    gene_to_mrna = {
        'SdhA': 'mRNA_9243',
        'SdhB': 'mRNA_44',
        'SdhC': 'mRNA_1581',
        'SdhD': 'mRNA_9290',
    }

    ref = {r.id: r.seq for r in parse(reference, 'fasta')}
    alt = {r.id: r.seq for r in parse(consensus, 'fasta')}

    cds_diffs = []
    non_cds_diffs = []
    for gene, cds_ranges in gene_to_cds_ranges.items():
        if gene not in gene_to_mrna:
            continue
        mrna = gene_to_mrna[gene]
        gene_cds_diffs = gene_seqs_to_cds_diffs(ref[gene], alt[gene], cds_ranges)
        gene_non_cds_diffs = gene_seqs_to_non_cds_diffs(ref[gene], alt[gene], cds_ranges)
        for diff in gene_cds_diffs + gene_non_cds_diffs:
            diff['gene'] = gene
        cds_diffs += gene_cds_diffs
        non_cds_diffs += gene_non_cds_diffs

    cds_diffs = pd.DataFrame(
        cds_diffs,
        columns=[
            'alt_codon_filled', 'ref_acid', 'alt_possible_acids', 'synonymous',
            'ref_codon', 'alt_codon', 'codon_index', 'indices', 'gene'
        ],
    )
    non_cds_diffs = pd.DataFrame(
        non_cds_diffs,
        columns=['pos', 'ref', 'alt', 'gene'],
    )
    non_cds_diffs[['gene', 'pos', 'ref', 'alt']].to_csv(non_cds_diffs_out, sep='\t', index=None)
    cds_diffs[[
        'gene', 'codon_index', 'indices', 'synonymous', 'ref_acid',
        'alt_possible_acids', 'ref_codon', 'alt_codon', 'alt_codon_filled',
    ]].to_csv(cds_diffs_out, sep='\t', index=None)

    pd.DataFrame({
        'sample': sample_name,
        'Amino Acid Changes': (cds_diffs.synonymous == 'no').sum(),
        'Codons With SNPs': cds_diffs.shape[0],
        'SNPs Outside CDS': non_cds_diffs.shape[0],
    }, index=[0]).to_csv(snp_summary_out, sep='\t', index=None)
