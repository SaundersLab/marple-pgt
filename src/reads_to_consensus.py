from transform import reads_list_to_consensus_with_report
from os.path import join, dirname, realpath
import argparse

if __name__ == '__main__':
    src_dir = dirname(realpath(__file__))
    marple_dir = dirname(src_dir)
    data_dir = join(marple_dir, 'data')
    reference_dir = join(data_dir, 'reference')

    parser = argparse.ArgumentParser(description='Align reads to reference genes then extract and concatenate exons')
    parser.add_argument(
        'relative_fastq_paths',
        type=str,
        nargs='+',
        help='Read files for aligning to the reference genes. Output for each read will be created in the same directory as the fastq file',
    )
    parser.add_argument(
        '--ref',
        help='Reference FASTA file to align reads to',
        default=join(reference_dir, 'crl75_205_genes.fasta'),
    )
    parser.add_argument(
        '--multiqc_config',
        help='MultiQC config file for creating the report',
        default=join(marple_dir, 'config', 'multiqc_config.yaml'),
    )
    parser.add_argument('--threads', type=int, help='Number of threads to use', default=2)
    parser.add_argument('--trim', help='Should FASTQ files be trimmed (yes/no)', default='yes')
    args = parser.parse_args()
    fastq_paths = [realpath(path) for path in args.relative_fastq_paths]
    out_dirs = [dirname(path) for path in fastq_paths]
    reads_list_to_consensus_with_report(
        fastq_paths=fastq_paths,
        reference=args.ref,
        out_dirs=out_dirs,
        multiqc_config=args.multiqc_config,
        threads=args.threads,
        trim=args.trim.lower().startswith('y'),
    )
