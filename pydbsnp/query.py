#===============================================================================
# query.py
#===============================================================================

"""query the dbSNP VCF data"""




# Imports ======================================================================

import re

from argparse import ArgumentParser
from pysam import TabixFile, VariantFile

from pydbsnp.env import BUILD_TO_VCF, BUILD_TO_RSID




# Constants ====================================================================

RSID_REGEX = re.compile('rs[1-9][0-9]+$')
COORD_REGEX = re.compile('.+:[0-9]+$')




# Functions ====================================================================

def parse_arguments():
    parser = ArgumentParser(description='query dbSNP VCF data')
    parser.add_argument(
        'variants',
        nargs='+',
        metavar='<rsid or chr:pos>',
        help='variant for which to query database'
    )
    parser.add_argument(
        '-r',
        '--reference-build',
        choices=('GRCh37', 'hg19', 'GRCh38', 'hg38'),
        default='GRCh38',
        help='reference build for coordinates'
    )
    return parser.parse_args()


def main():
    args = parse_arguments()
    variant_file = VariantFile(BUILD_TO_VCF[args.reference_build])
    print(variant_file.header)
    for variant in args.variants:
        if COORD_REGEX.match(variant):
            chrom, pos = id.split(':')
            pos = int(pos)
            tbx = TabixFile(BUILD_TO_VCF[args.reference_build])
            for row in tbx.fetch(chrom, pos, pos):
                print(str(row))
        elif RSID_REGEX.match(variant):
            for coord in rsid_to_coordinates():
        else:
            raise RuntimeError('Improperly formatted query')
    pass
