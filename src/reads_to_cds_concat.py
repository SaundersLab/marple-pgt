from transform import reads_list_to_cds_concat_with_report
from os.path import join, dirname, realpath
import argparse

if __name__ == '__main__':
    src_dir = dirname(realpath(__file__))
    marple_dir = dirname(src_dir)
    data_dir = join(marple_dir, 'data')
    reference_dir = join(data_dir, 'reference')

    parser = argparse.ArgumentParser(
        description='Align reads to reference genes then extract and concatenate CDS',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        'relative_fastq_paths',
        type=str,
        nargs='+',
        help='Read files for aligning to the reference genes. Output for each read will be created in the same directory as the fastq file',
    )
    parser.add_argument(
        '--ref',
        help='Reference FASTA file to align reads to',
        default=join(reference_dir, '280_genes.fa'),
    )
    parser.add_argument(
        '--gff',
        help='Annotation file giving the position of exons within the reference genes. Genes should be specificied as if everything was on the + strand.',
        default=join(reference_dir, '280_genes.gff3'),
    )
    parser.add_argument(
        '--multiqc_config',
        help='MultiQC config file for creating the report',
        default=join(marple_dir, 'config', 'multiqc_config.yaml'),
    )
    parser.add_argument(
        '--primers',
        help='Primers spreadhseet with gene name and l_boundary_padding and r_boundary_padding',
        default=join(marple_dir, 'data', 'primers', '282_pairs_280_genes_4_pools.xlsx'),
    )
    parser.add_argument('--threads', type=int, help='Number of threads to use', default=2)
    parser.add_argument('--trim', help='Should FASTQ files be trimmed (yes/no)', default='yes')
    parser.add_argument('--max_read_length', type=int, help='Maxmium length of read in FASTQ to keep', default=4_000)
    parser.add_argument('--min_snp_depth', type=int, help='Minimum depth required to call a SNP', default=20)
    args = parser.parse_args()
    fastq_paths = [realpath(path) for path in args.relative_fastq_paths]
    out_dirs = [dirname(path) for path in fastq_paths]
    reads_list_to_cds_concat_with_report(
        fastq_paths=fastq_paths,
        reference=args.ref,
        gff=args.gff,
        out_dirs=out_dirs,
        multiqc_config=args.multiqc_config,
        primers_path=args.primers,
        threads=args.threads,
        trim=args.trim.lower().startswith('y'),
        max_read_length=args.max_read_length,
        min_snp_depth=args.min_snp_depth,
    )
