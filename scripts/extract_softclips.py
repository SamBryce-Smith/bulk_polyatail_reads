#!/usr/bin/env python3
import pysam
import sys
import os
import argparse


def main(input_bam: str,
         output_bam: str,
         min_clip_len: int = 2):
    
    # number of primary alignments considered
    primary_count = 0
    # read has a soft-clip at either end (e.g. possible tail read)
    sclip_count = 0

    with pysam.AlignmentFile(input_bam, "rb") as bamfile:
        
        out_bam = pysam.AlignmentFile(output_bam, "wb", template=bamfile)

        i = 0

        for alignment in bamfile.fetch():

            i += 1
            if i % 1000000 == 0: print(f'Number of records parsed: {i}')

            # skip if not primary alignment
            if alignment.is_secondary: continue
            primary_count += 1
            
            # skip if no soft clips
            if not "S" in alignment.cigarstring: continue

            # extract tuple of (op, len) for each element in cigar string
            cigtup = alignment.cigartuples

            # [(op, len)] - check if 1st or last op is a soft-clip & passes minimum length criteria
            # 4 = soft-clip operation key
            sclip_start = cigtup[0][0] == 4 and cigtup[0][1] >= min_clip_len
            sclip_end = cigtup[-1][0] == 4 and cigtup[-1][1] >= min_clip_len

            # skip read if neither end has a soft clip of min length
            if not sclip_start and not sclip_end: continue
            sclip_count += 1

            # write alignment to output file
            out_bam.write(alignment)

    

    # all reads iterated - print stats
    print("Number of primary alignments considered:", primary_count)
    print("Read has a soft-clip at either end of at least min length (e.g. possible tail read):", sclip_count)
    print("Fraction of primary alignments with soft-clips of at least minimum length", sclip_count / primary_count)
    print("Indexing output bam file:", output_bam)

    # index output BAM
    out_bam.close()
    pysam.index(output_bam)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extract reads containing terminal soft clips from short-read RNA-seq alignments',
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter, # Add defaults to end of help strings
                                      )

    parser.add_argument('-b', '--bam', type=str, help='Path to BAM file')
    parser.add_argument('-m','--min_clip_len', type=int, default=2, help='Minimum soft-clip length for read to be written')
    parser.add_argument('-o','--output_bam', type=str, help='Name of output file containing soft-clipped reads (primary alignments, at leas min_clip_len)')

    if len(sys.argv) == 1:
        parser.print_help()
        parser.exit()

    args = parser.parse_args()
    main(args.bam, args.output_bam, args.min_clip_len)