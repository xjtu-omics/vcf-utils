# vcf-utils

## tools

**compare**: finds all events in the first file that have a comparable event (same type and size, similar location) in the second file, and returns it as the output of a third file

**count_events**: tool to count the number of events (insertions, deletions, replacements, SNPs, whatever) in a VCF file

**del_corr**: to be used for VCF files that have an uncommon deletion notation, like "chr1 10 AT C", which is weird since you would expect the alt of a deletion to be identical to the first base of the reference. If the basic problem is unconventional VCF-creating software, then del_corr may help to transform the deletions into the more generic format of "chr1 9 CA C"

**eventizer**: takes an input VCF file, and produces an output text file that can be used for selecting or removing specific events or loci, for example as filter_events does.

**filter_events**: takes an input file, a file that contains a list of events (like “chr1:10:A:AT”) and removes all events that are NOT in the list, and writes the result to a third file.

**filter_eventtypes**: takes an input VCF file and returns all events that are of a certain class (INS/DEL/SNP/ALL) and have certain minimum and maximum sizes.

**find_duplicates**: takes an input file, writes the list of duplicates to std::out in “chr1:10:A:AT” format

**find_mlma**: takes an input file, and writes a list of 'duplicate loci' (so events that have same chromosome and position, but may have different alt (or ref!) alleles in “chr1:10” fomat). Useful for pindel output VCFs

**find_uncrowded**: finds events that have no other events in a certain window around them, such as no other event within the 100 bases before and after the event. This can help eliminate events in repetitive regions, which are relatively likely to be inaccurate.

**fuse**: fuses two VCF files that describe the same sample(s) (actually merges them)

**indel_split**: splits an indel file into an insertion file and a deletion file (only pure insertions and deletions, no replacements!)

**left_align**: aligns the events in a VCF file to the left (not all pipelines produce left-aligned events)

**min_bedmaker**: Turns a VCF file into a BED file, defining the start of each interval as the start position in the VCF, and the end as the start position + the event length (1 for SNPs, 2 for 1-base indels, etc)

**read_reference**: not strictly something that manipulates a VCF, this is a quick tool to check the sequence of the reference genome at a certain position, useful for finding the context of an event in the VCF

**remove_double_alts**: takes an input file (and output file), removes events containing GATK-like double alts (“A TA,TAA”) from the input file

**remove_events**: takes an input file, a file that contains a list of events (like “chr1:10:A:AT” or “chr1:10”), removes all events that are in the list, and writes the result to a third file

**remove_homref**: takes an input file, removes hom refs (whether encoded like “.” or like “0/0”)

**size_ass**: takes an input file, outputs a file containing rows in the form of “1 10232 Pindel deletion”

**sort**: sorts a VCF file into the sequence chr1, chr2...chr22, chrX, chrY, chrM

**standardize**: basically helps transform a 'normal' PacBio file (with \<INS\> and \<DEL\> alt labels) into something with explicit REF and ALT fields. Note that if the file also has a weird format for deletions, it is better to use del_corr instead

**uniquify**: If any event (chrom-pos-REF-ALT) occurs multiple times in a file (for example after merging VCF files), only keeps one copy of the event, ensuring that all events reported by the VCF-file are unique.

**uniquify_loci**: If any locus (chrom-pos) occurs in multiple events in the VCF, only keep the first event at the locus. Note that this will eliminate more events than the similar utility 'uniquify'.

**unravel_alts**: splits a mixed alt in a VCF file (like A AT,AGC) into separate lines. Can be useful when processing GATK VCF files.

## usage

To compile/create the set of utilities, use
./runme

Running any tool without arguments (for example "./remove_homref") shows its usage instructions.

