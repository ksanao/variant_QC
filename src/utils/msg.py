import numpy as np
import seaborn as sns
import pandas as pd
import warnings


def get_sam_stat(df, filename, attr):
    val = df.loc[df.Filename.isin([filename]) & df.Attribute.isin([attr]),'Value'].tolist()[0]
    return val


def print_samstats_message(stats):
    for filename in stats.Filename.unique():
        n_paired = get_sam_stat(stats, filename, 'reads_paired:')
        insert_size = get_sam_stat(stats, filename, 'insert_size_average:')
        ins_dev = get_sam_stat(stats, filename, 'insert_size_standard_deviation:')
        read_len = get_sam_stat(stats, filename, 'average_length:')

        print(f"\nSequencing library {filename}")
        if n_paired>0:
            print(f"Paired {read_len}bp reads with {insert_size}+/-{ins_dev}bp average insert size")
        else:
            print(f"Single reads with {read_len}bp average read length")


def print_samstats_mismatch(stats):
    for filename in stats.Filename.unique():
        mapped_bp = get_sam_stat(stats, filename, 'bases_mapped:')
        mismatch_bp = get_sam_stat(stats, filename, 'mismatches:')
        err = get_sam_stat(stats, filename, 'error_rate:')
        avg_qual = get_sam_stat(stats, filename, 'average_quality:')
        dups = get_sam_stat(stats, filename, 'reads_duplicated:')

        print(f"\nSequencing library {filename}")
        print(f"Mapped: {mapped_bp}bp\nMismatched: {mismatch_bp}bp")
        print(f"Error-rate: {err}")
        print(f"Average read quality: {avg_qual}")
        print(f"Duplicate reads: {dups}")


def print_sams_pairstats(stats):
    for filename in stats.Filename.unique():
        total_pairs = get_sam_stat(stats, filename, 'reads_paired:')
        proper_pairs = get_sam_stat(stats, filename, 'reads_properly_paired:')
        pairs_diff_chr = get_sam_stat(stats, filename, 'pairs_on_different_chromosomes:')
        inw = get_sam_stat(stats, filename, 'inward_oriented_pairs:')
        out = get_sam_stat(stats, filename, 'outward_oriented_pairs:')
        other = get_sam_stat(stats, filename, 'pairs_with_other_orientation:')
        print(f"\nSequencing library {filename}")
        print(f"Total paired reads: {total_pairs}")
        print(f"Proper pairs: {proper_pairs} ({np.round((proper_pairs/total_pairs)*100)}%)")
        print(f"Pairs on different chromsomes {pairs_diff_chr}")
        print(f"Inward orientation: {inw}")
        print(f"Outward orientation: {out}")
        print(f"Other orientation: {other}")


