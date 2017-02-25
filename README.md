# vcf-utils

## tools

**filter_events**: takes an input file, a file that contains a list of events (like “chr1:10:A:AT”) and removes all events that are NOT in the list, and writes the result to a third file.
**find_duplicates**: takes an input file, writes the list of duplicates to std::out in “chr1:10:A:AT” format
**find_mlma**: takes an input file, and writes a list of 'duplicate loci' (so events that have same chromosome and position, but may have different alt (or ref!) alleles in “chr1:10” fomat). Useful for pindel output VCFs.
**fuse**: fuses two VCF files that describe the same sample(s) (actually merges them).
**indel_split**: splits an indel file into an insertion file and a deletion file (only pure insertions and deletions, no replacements!).
**remove_double_alts**: takes an input file (and output file), removes events containing GATK-like double alts (“A TA,TAA”) from the input file
**remove_events**: takes an input file, a file that contains a list of events (like “chr1:10:A:AT” or “chr1:10”), removes all events that are in the list, and writes the result to a third file.
**remove_homref**: takes an input file, removes hom refs (whether encoded like “.” or like “0/0”)
**size_ass**: takes an input file, outputs a file containing rows in the form of “1 10232 Pindel deletion”
**sort**: sorts a VCF file into the sequence chr1, chr2...chr22, chrX, chrY, chrM

