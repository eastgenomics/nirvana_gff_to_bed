#!/usr/bin/python
"""
Convert a Nirvana Refseq gff file to bed file

The output bed only contains regions that satisify the following critera:
       Type is 'CDS'
       transcript_type is 'coding'
       transcript_id starts with 'NM_'

Adapted from ccbg_toolbox import_gff created by Kim Brugger
This is the code which previously populated the genetics ark database
with regions in the gff. Database related code has been removed and
instead we print a bed to stdout

Pipe output into
| sort -k1V -k2n -k3n -k4
for a sorted bed

Matt Garner 200504
"""

import sys
import os
import re
import argparse
import pprint
pp = pprint.PrettyPrinter(indent=4)
import gzip


def gff_to_bed_file(gff_file, create_exon_file, flank=0):
    """The output bed only contains regions that satisify the following critera:
       Type is 'CDS'
       transcript_type is 'coding'
       transcript_id starts with 'NM_'

    Args:
      gff_file (str): name of gff file
      gff_file (int): number of flanking bases to append to each region (default: 0)

    Returns:
      None

    """

    # Just so we can handle both compressed and uncompressed gff files
    # The reason not to use pysam for this is the gff needs to be
    # index for it to work.

    assert int(flank) >= 0, "Flank must be an integer >= 0"

    if ( re.search(".gz", gff_file)):
         gff_fh = gzip.open(gff_file, 'rb')
    else:
         gff_fh = open(gff_file, 'r')

    for line in gff_fh:
        if type(line) is bytes:
            line = line.decode("utf-8")

        line = line.rstrip("\n").rstrip("\r")
        line = line.rstrip("; ")

        fields = line.split("\t")
        fields[ 0 ] = re.sub(r'^chr', '', fields[ 0 ])
        fields[ 3 ] = int ( fields[ 3 ])  # Start
        fields[ 4 ] = int ( fields[ 4 ])  # End

        # Only interested in coding regions
        if ( fields[ 2 ] != 'CDS'):
            continue

        # Make a dict of annotations
        annots = fields[ 8 ].split("; ")
        annot_dict = {}

        for annot in annots:
            key, value = annot.split(" ", 1)
            value = re.sub(r'\"', '', value)
            annot_dict[ key ] = value

        # Only interested in coding and NM_ accessions (not XM_ etc)
        if ( annot_dict[ 'transcript_type'] != 'protein_coding' or 
             not re.match(r'^NM_', annot_dict[ 'transcript_id'])):
            continue

        chrom = fields[0]
        start = int(fields[3]) - 1 - int(flank)  # -1 since bed is 0 based
        end = int(fields[4]) + int(flank)

        if create_exon_file:
            print("\t".join(map(str, [
                chrom,
                start,
                end,
                annot_dict["gene_name"],
                annot_dict['transcript_id'],
                annot_dict["exon_number"]
            ])))
        else:
            print("\t".join(map(str,
                            [chrom, start, end, annot_dict['transcript_id']])))


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description='Convert a Nirvana Refseq gff file to bed file')
    parser.add_argument(
        'gff', metavar='gff', nargs=1,  help="gff file")
    parser.add_argument(
        'flank', metavar='flank', nargs='?', type=int, default=0, help="flank")
    parser.add_argument(
        "-e", "--create_exon",
        default=False,
        action="store_true",
        required=False,
        help="Pass this option if you want to create the exon file"
    )

    args = parser.parse_args()

    gff_to_bed_file(args.gff[0], args.create_exon, args.flank)
