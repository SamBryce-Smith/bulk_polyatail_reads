#!/usr/bin/env python3

import pandas as pd
import numpy as np
import pysam
import sys
import os
import argparse
from collections import Counter, defaultdict
from shutil import rmtree

# TODO: consider rename - more accurate to say it gives you transcribed strand?
def get_transcribed_strand(read: pysam.AlignedSegment,
                     strandedness: str = ["--rf", "--fr", "--unstranded"]):
    
    assert strandedness in ["--rf", "--fr", "--unstranded"]

    if strandedness == "--unstranded":
        # can't determine for this library prep
        return 0

    # determine alignment strand
    elif strandedness == "--rf":
        # reverse stranded - if read 1 then maps to strand of origin, otherwise opposite strand
        if read.is_reverse:
            # maps to minus strand - if read1 then opposite strand of origin, otherwise same strand
            if read.is_read1:
                return 1
            else:
                return -1
    
        # read maps to plus strand - if read1 then opposite strand of origin, otherwise same strand
        elif read.is_read1:
            return -1
    
        elif read.is_read2:
            return 1
        
    else:
        # forward stranded - if read1 then maps to strand of origin, otherwise opposite
        if read.is_read1:
            return -1 if read.is_reverse else 1
        
        else:
            return 1 if read.is_reverse else -1


def get_soft_clip_sequence(read: pysam.AlignedSegment,
                           softclip_len: int,
                           end: str = ["l", "r"]) -> str:
    '''Extract soft-clipped read sequence at specified terminus of read

    Parameters
    ----------
    read : pysam.AlignedSegment
        _description_
    softclip_len : int
        length of soft clipped region to extract
    end : str, optional
        extract sequence from left-hand ('l') or right-hand ('r') side of read, by default ["l", "r"]

    Returns
    -------
    str
        _description_
    '''

    assert end in ["l", "r"]

    if end == "r":
        # cut from end of sequence
        # query_sequence returns read sequence as stored in BAM file including soft clip
        # If aligned to reverse strand then read sequence is reverse complemented
        # Desired to match with orientation of cigar string and simplify counting of tail nucleotides (regardless of alignment orientation have same nucleotide)
        ohang_seq = read.query_sequence[-softclip_len:]

    else:
        # cut from start of sequence
        ohang_seq = read.query_sequence[:softclip_len]

    return ohang_seq


def count_tail_soft_clip(softclip_seq: str,
                         tail_nucleotide: str = ["A", "T"],
                         ):
    '''Calculate the number of polyA tail nucleotides in a soft-clipped region of read sequence

    Parameters
    ----------
    softclip_seq : str
        soft-clipped region of read sequence, extracted by e.g. get_soft_clip_sequence()
    tail_nucleotide : str, optional
        which nucleotide to count as the expected tail nucleotide, by default ["A", "T"]
    '''

    assert tail_nucleotide in ["A", "T"]

    return softclip_seq.count(tail_nucleotide)


def main(bam_path: str,
         output_bam: bool,
         output_prefix: str,
         output_partition_cols: list,
         strandedness: str = ["--rf", "--fr", "--unstranded"],
         min_clip_len: int = 2,
         min_frac_a: float = 0.0,
         force_overwrite: bool = True,
         ):
    
    # Check existence of output parquet file in case of restrictions with overwriting
    if os.path.exists(output_prefix + ".parquet"):
        if force_overwrite:
            print(f"Removing existing .parquet file for provided output prefix: {output_prefix + '.parquet'}")
            try:
                os.remove(output_prefix + ".parquet")
            except IsADirectoryError:
                rmtree(output_prefix + ".parquet")
        else:
            raise FileExistsError(f"Output .parquet file - {output_prefix + '.parquet'} - already exists. Provide a different prefix or permit overwriting")

    # number of primary alignments considered
    primary_count = 0

    # For some aligners have found that read.cigarstring returns a TypeError - track these counts here
    cigar_incompat = 0

    # read has a soft-clip at either end (e.g. possible tail read)
    sclip_count = 0
    # read has a soft-clip at both ends
    sclip_both_count = 0
    # terminal soft clip incompatible with alignment & inferred transcribed strand
    sclip_incompat = 0

    # number of reads failing minimum A content filter
    sclip_min_a = 0

    # number of reads with valid soft clip
    sclip_valid = 0

    # valid soft clips - track counts of num soft clipped reads
    # length of soft-clip/overhang - key
    # num As, frac As - Counter
    # .update([val]) to update with incoming value
    sclip_stats = defaultdict(lambda: (Counter(), Counter()))

    # Output column name order for parquet file (BED-like format)
    outcol_order =  ["Chromosome", "Start", "End", "read_name", "overhang_len", "Strand", "frac_As", "overhang_seq"]

   
    with pysam.AlignmentFile(bam_path, "rb") as bamfile:

        if output_bam:
            out_bam = pysam.AlignmentFile(output_prefix + ".bam", "wb", template=bamfile)

        # track global number of records parsed
        i = 0

        ref_names = bamfile.references
        # num_refs = len(ref_names)
        # print(f"Number of reference sequence names in BAM file - {num_refs}")
        
        for ref in ref_names:
            
            # Initialise dict storing soft clip read info - store final aligned position, overhang len, tail nuc content & overhang seq
            sclip_reads = defaultdict(list)

            for alignment in bamfile.fetch(contig=ref):

                i += 1
                if i % 1000000 == 0: print(f'Number of records parsed: {i}')

                # skip if not primary alignment
                if alignment.is_secondary: continue
                primary_count += 1
                
                # skip if no soft clips
                try:
                    if not "S" in alignment.cigarstring: continue
                
                # some aligners return TypeError: argument of type 'NoneType' is not iterable
                except TypeError:
                    sclip_incompat +=1
                    # either no cigar or can't infer - skip
                    continue

                # extract tuple of (op, len) for each element in cigar string
                cigtup = alignment.cigartuples

                # [(op, len)] - check if 1st or last op is a soft-clip & passes minimum length criteria
                # 4 = soft-clip operation key
                sclip_start = cigtup[0][0] == 4 and cigtup[0][1] >= min_clip_len
                sclip_end = cigtup[-1][0] == 4 and cigtup[-1][1] >= min_clip_len

                # skip read if neither end has a soft clip
                if not sclip_start and not sclip_end: continue
                sclip_count += 1


                # now check inferred transcribed strand of the read (1 = plus, -1 = minus, 0 = unstranded)
                tx_strand = get_transcribed_strand(alignment, strandedness)
    
                # check if clips at both end of read
                # If stranded protocol, can use inferred txibed strand to get putative 3'end
                # If unstranded, can't infer txibed strand to get the 3'end.
                # Also don't want to drastically change logic by checking both ends for softclips (then a case of how to decide in favour of reporting 1 or the other?)
                if sclip_start and sclip_end:
                    sclip_both_count +=1
                    if tx_strand == 0:
                        continue

                # use alignment strand and transcribed strand to get what should be the correct soft-clip region
                # Transcribed strand/strand of origin tells us which side of the alignment/fragment should be 3'most
                # Plus strand origin = RHS = 3'most
                # Minus strand origin = LHS = 3'most

                # Check if transcribed strand/strand of origin is consistent with soft-clip region
                # (if unstranded, this cannot be determined and is effectively skipped)
                if not (tx_strand == 0) and not (tx_strand == 1 and sclip_end) and not (tx_strand == -1 and sclip_start):
                    # soft clip incompatible with inferred strand of origin + being at 3'end)
                    sclip_incompat +=1
                    continue
                
                # Now use strand of origin to define side of alignment + length of softclip to extract 
                sclip_end = "r" if (tx_strand in [0, 1] and sclip_end) else "l"
                sclip_len = cigtup[-1][1] if sclip_end == "r" else cigtup[0][1]

                # Strand of origin/transcribed strand then tells us what nucleotide we should expect soft clipped tail to contain
                # (ASSUMING WE'RE CONSIDERING THE 'QUERY SEQUENCE' as stored in the BAM file (reverse complemented if aligned to minus strand))
                # plus strand aligned = As
                # minus strand aligned = Ts
                # unstranded - A = right-most clip, T if left-most clip
                tail_nuc = "A" if tx_strand == 1 or (tx_strand == 0 and sclip_end == "r") else "T"

                # Count number of tail nucleotides in soft clip region
                ohang_seq = get_soft_clip_sequence(alignment, sclip_len, sclip_end)
                n_bases = count_tail_soft_clip(ohang_seq, tail_nuc)
                frac_bases = n_bases / sclip_len

                # update stats (think good to track for all reads irrespective of filter)
                sclip_stats[str(sclip_len)][0].update([n_bases])
                sclip_stats[str(sclip_len)][1].update([frac_bases])

                # Check passing minimum A content filter
                if not frac_bases >= min_frac_a:
                    sclip_min_a +=1
                    continue

                # Passed all filters - now store required read-level info
                sclip_valid +=1

                # update read dict
                read_num = str(1) if alignment.is_read1 else str(2)

                # extract 3'most aligned position of read to mark putative cleavage site
                # plus strand = right-most position (reference end points 1 past last aligned residue, -1 to get to aligned position)
                # minus strand = left-most position
                start_coord = alignment.reference_end - 1 if tx_strand == 1 else alignment.reference_start
                end_coord =  alignment.reference_end if tx_strand == 1 else alignment.reference_start + 1

                sclip_reads["Start"].append(start_coord)
                sclip_reads["End"].append(end_coord)
                sclip_reads["Strand"].append("+" if tx_strand == 1 else "-")
                sclip_reads["read_name"].append(alignment.query_name + ";" + read_num)
                sclip_reads["overhang_len"].append(sclip_len)
                sclip_reads["frac_As"].append(frac_bases)
                sclip_reads["overhang_seq"].append(ohang_seq)

                if output_bam:
                    out_bam.write(alignment)

                i += 1
            
            # finished for ref name
   
            if len(sclip_reads) == 0:
                # no extracting reads for ref, skip to the next one
                continue
            # convert sclip reads dict to df & add reference name
            reads_df = pd.DataFrame.from_dict(sclip_reads, orient="columns")
            reads_df["Chromosome"] = ref

            # to support partitioning by overhang length (but not massively amplify number of partitions)
            # assign all overhangs above certain length (here = 10) to a common 'large' value
            if "overhang_len_bin" in output_partition_cols:

                reads_df["overhang_len_bin"] = np.where(reads_df["overhang_len"].gt(10), 10, reads_df["overhang_len"])
                outcol_order.append("overhang_len_bin")

            if "frac_As_bin" in output_partition_cols:
                n_bins = 10
                reads_df["frac_As_bin"] = pd.cut(reads_df["frac_As"], bins=n_bins, labels=[str((i*5 / 100)) for i in range(n_bins)], right=False, include_lowest=True)
                outcol_order.append("frac_As_bin")

            # output in format compatible with bgzip and tabix
            # i.e. TSV file, no header (/ header with comment line) - keep cols in BED-like format
            # position sorted
            reads_df = reads_df.sort_values(by=["Start"], ascending=True)[outcol_order]

            # Output to parquet, using fastparquet to append if necessary
            try:
                reads_df.to_parquet(output_prefix + ".parquet", engine="fastparquet",  index=False, partition_cols=output_partition_cols, append=True)
            except FileNotFoundError:
                # first write so output as is
                reads_df.to_parquet(output_prefix + ".parquet", engine="fastparquet",  index=False, partition_cols=output_partition_cols)
        
        # all reads iterated - print stats
        print("Number of primary alignments considered:", primary_count)
        print("Read has a soft-clip at either end (e.g. possible tail read):", sclip_count)
        print("Read has a soft-clip at both ends:", sclip_both_count)
        print("Terminal soft clip incompatible with alignment and inferred transcribed strand (valid for stranded flags only):", sclip_incompat)
        print(f"Number of reads failing minimum fraction A content threshold of {min_frac_a}: {sclip_min_a}")
        print("Number of reads with valid soft clip:", sclip_valid)
        print("Fraction of primary alignments with valid soft-clips of minimum length:", sclip_valid / primary_count)
        print("Number of primary alignments with an incompatible cigar string:", cigar_incompat)
        print("Fraction of primary alignments with an incompatible cigar string:", cigar_incompat / primary_count)
        
        # write stats to TSV
        with open(output_prefix + ".alignment_stats.tsv", "w") as outfile:
            outfile.write("stat\tvalue\n")
            outfile.write("primary_alignments\t" + str(primary_count) + "\n")
            outfile.write("soft_clip_any\t" + str(sclip_count) + "\n")
            outfile.write("soft_clip_both\t" + str(sclip_both_count) + "\n")
            outfile.write("soft_clip_incompatible\t" + str(sclip_incompat) + "\n")
            outfile.write("soft_clip_min_A\t" + str(sclip_min_a) + "\n")
            outfile.write("soft_clip_valid\t" + str(sclip_valid) + "\n")
            outfile.write("fraction_valid\t" + str(sclip_valid / primary_count) + "\n")
            outfile.write("cigar_incompatible\t" + str(cigar_incompat) + "\n")
            outfile.write("fraction_cigar_incompatible\t" + str(cigar_incompat / primary_count) + "\n")

        # Convert counter of number of overhang lens + their different tail contents
        sclip_stats_counts_df = (pd.concat({k: pd.DataFrame(list(v[0].items()),
                                        columns=["number_tail_nucleotides", "count"])
                                        for k,v in sclip_stats.items()},
                                        keys=sclip_stats.keys(),
                                        names=["overhang_length", None]
                                        )
                                .reset_index("overhang_length")
                                .reset_index(drop=True)
                                )

        # index output BAM
        if output_bam:
            out_bam.close()
            pysam.index(output_prefix + ".bam")
            print(pysam.index.get_messages())
        
        # output overhang length + tail content global statistics
        print("Outputting global soft clip length and tail content statistics to TSV...")
        sclip_stats_counts_df = sclip_stats_counts_df.astype({"overhang_length": int})
        sclip_stats_counts_df.sort_values(by=["overhang_length", "number_tail_nucleotides"]).to_csv(output_prefix + ".global_softclip_tail_counts.tsv", sep="\t", header=True, index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extract polyA tail containing reads from short-read RNA-seq alignments',
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter, # Add defaults to end of help strings
                                      )

    parser.add_argument('-b', '--bam', type=str, default=argparse.SUPPRESS, help='Path to BAM file')
    parser.add_argument('-s','--strandedness', type=str, choices=["rf", "fr", "unstranded"], default="rf", help='Strandedness option. If unstranded, reads are searched for As on rightmost or Ts on leftmost softclips, and reads with soft clips on both ends are discarded')
    parser.add_argument('-m','--min_clip_len', type=int, default=2, help='Minimum soft-clip length for read to be considered')
    parser.add_argument('-f', '--min-fracA', dest="min_frac_a", type=float, default=0, help='Minimum A content to retain as a valid soft clip (in fractional terms)')
    parser.add_argument('-p', '--partition-cols', nargs='+', choices=["Chromosome", "overhang_len", "Strand", "frac_As_bin", "overhang_len_bin"], default=['Chromosome'], help='Columns on which to partition in output parquet file. If multiple, pass consecutively and space-separated')
    parser.add_argument("--no-output-bam", dest="output_bam", action="store_false", default=argparse.SUPPRESS, help="Do not output a BAM file containing soft clipped reads. If not provided, a BAM file is produced")
    parser.add_argument('-o','--output_prefix', type=str, default="polyA_tail_reads", help='Prefix for output file names')

    if len(sys.argv) == 1:
        parser.print_help()
        parser.exit()

    args = parser.parse_args()
    main(args.bam, args.output_bam, args.output_prefix, args.partition_cols, "--" + args.strandedness, args.min_clip_len, args.min_frac_a, True)


