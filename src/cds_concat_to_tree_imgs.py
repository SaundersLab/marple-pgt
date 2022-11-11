from os.path import dirname, join, realpath
from transform import cds_concat_paths_to_tree_imgs
import argparse

if __name__ == '__main__':
    src_dir = dirname(realpath(__file__))
    marple_dir = dirname(src_dir)
    data_dir = join(marple_dir, 'data')

    parser = argparse.ArgumentParser(description='Create and visualise tree with new samples')
    parser.add_argument(
        'relative_cds_concat_paths',
        type=str,
        nargs='*',
        help='Concatenated coding sequences file for new samples to add to the tree',
        default=[],
    )
    parser.add_argument(
        '--meta',
        help='Path to spreadsheet containing isolate metadata and styles',
        default=join(data_dir, 'metadata_92.xlsx')
    )
    parser.add_argument(
        '--start',
        help='Starting tree input, i.e. concatenated CDS file for previously sequenced samples',
        default=join(data_dir, '72_samples_276_genes_cds.fa.gz')
    )
    parser.add_argument('--name', help='Name to use as a prefix for all created tree files',)
    parser.add_argument('--out_dir', help='Directory to create tree files in', default=realpath('.'))
    parser.add_argument('--threads', type=int, help='Number of threads to use for RAxML', default=1)
    parser.add_argument('--img_fmt', help='Format to output tree images as', default='pdf')

    args = parser.parse_args()
    cds_concat_paths_to_tree_imgs(
        cds_concat_paths=[realpath(path) for path in args.relative_cds_concat_paths],
        starting_tree_input=args.start,
        tree_name=args.name,
        out_dir=args.out_dir,
        metadata_path=args.meta,
        n_threads=args.threads,
        img_fmt=args.img_fmt,
    )