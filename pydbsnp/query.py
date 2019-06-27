#===============================================================================
# query.py
#===============================================================================

"""query the dbSNP VCF data"""




# Imports ======================================================================

import pysam
import re

from argparse import ArgumentParser

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
    for variant in args.variants:
        if RSID_REGEX.match(variant):
            query_by_rsid(variant)
        elif COORD_REGEX.match(variant):
            query_by_coord(variant)
        else:
            raise RuntimeError('Improperly formatted query')
    pass
