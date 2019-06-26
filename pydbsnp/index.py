#===============================================================================
# index.py
#===============================================================================

"""Index dbSNP VCFs by rsid"""



# Imports ======================================================================

import argparse
import gzip
import os
import os.path
import pysam
import subprocess
import tempfile

from functools import partial
from multiprocessing import Pool

from pydbsnp.env import VCF_GRCH37, VCF_GRCH38, BED_GRCH37, BED_GRCH38




# Functions ====================================================================

def reformat_sort(
    input_vcf_path,
    output_bed_path,
    quiet=False,
    temp_dir=None
):
    if not quiet:
        print(
            f'Reformatting database {input_vcf_path} and sorting by RSID. '
            'This will take a while. Reformatted and sorted data will be '
            f'written to {output_bed_path}.'
        )
    with open(output_bed_path, 'wb') as f:
        with subprocess.Popen(
            ('zcat', input_vcf_path), stdout=subprocess.PIPE
        ) as zcat:
            with subprocess.Popen(
                (
                    'awk', '-v', 'OFS=\t',
                    (
                        '!/##/ && !/#CHROM/ '
                        '{sub(/rs/, "", $4); print "rs", $4 - 1, $4, $1, $3}'
                    ),
                ),
                stdin=zcat.stdout,
                stdout=subprocess.PIPE
            ) as awk:
                with subprocess.Popen(
                    (
                        'sort',
                        '-k3,3',
                        '-n',
                        '-T', temp_dir or os.path.dirname(output_bed_path)
                    ),
                    stdin=awk.stdout,
                    stdout=subprocess.PIPE
                ) as sort:
                    subprocess.run(
                        ('bgzip', '--stdout'),
                        stdin=sort.stdout,
                        stdout=f
                    )
                

def index(bed_file_path, quiet=False):
    if not quiet:
        print(f'Indexing {bed_file_path}.')
    subprocess.run(
        (
            'tabix',
            '--csi',
            '--sequence', '1',
            '--begin', '2',
            '--end', '3', 
            bed_file_path
        )
    )


def reformat_sort_index(
    input_vcf_path,
    output_bed_path,
    quiet=False,
    temp_dir=None
):
    reformat_sort(
        input_vcf_path,
        output_bed_path,
        quiet=quiet,
        temp_dir=temp_dir
    )
    index(output_bed_path, quiet=False)
    

def parse_arguments():
    parser = argparse.ArgumentParser(
        description='index dbSNP VCF data by rsid'
    )
    parser.add_argument(
        '--quiet',
        action='store_true',
        help='suppress printed status updates'
    )
    parser.add_argument(
        '--processes',
        choices=(1, 2),
        default=1,
        help='set to 2 to index GRCh37 and GRCh38 data in parallel'
    )
    parser.add_argument(
        '--tmp-dir',
        metavar='<path/to/tmp/dir/>',
        help='directory for temporary files'
    )
    return parser.parse_args()


def main():
    args = parse_arguments()
    with Pool(processes=args.processes) as pool:
        pool.starmap(
            partial(
                reformat_sort_index,
                quiet=args.quiet,
                temp_dir=args.tmp_dir
            ),
            (
                (vcf, bed)
                for vcf, bed
                in ((VCF_GRCH37, VCF_GRCH38), (BED_GRCH37, BED_GRCH38))
                if os.path.isfile(vcf)
            )
        )
