#!/usr/bin/env python3
'''
The "subsample_test_global_folder" should contain each value of n "max mapped reads" as its own nested folder,
and "subsample_test_global_folder/nested_full_folder" should be the folder that contains the non-subsampled runs.

./subsample_test.py Indel_test_for_Karthik_20210322_210309_A00953_0250_BH2M2TDRXY_results 2021-03-23_22-08-33_pe2021-03-23_22-08-33_pe
'''
# global constants
IGNORE_FIRST = 100 # ignore the first IGNORE_START positions of reference genome
IGNORE_LAST  =  50 # ignore the last IGNORE_END positions of reference genome
REF_FIRST_ORF_START_0BASED = 265
REF_LAST_ORF_END_0BASED = 29673
REF_LENGTH = 29903 # total length of reference
ORF_START_TO_SUBSTRING_START = REF_FIRST_ORF_START_0BASED - IGNORE_FIRST # inclusive
ORF_END_TO_SUBSTRING_END = REF_LENGTH - IGNORE_LAST - REF_LAST_ORF_END_0BASED # exclusive

# imports
from csv import reader as csv_reader
from glob import glob
from json import load as json_load
import matplotlib.pyplot as plt
from os.path import isdir
from sys import argv

# check user arguments
if len(argv) != 3:
    print("USAGE: %s <subsample_test_global_folder> <nested_full_folder>" % argv[0]); exit(1)

# see which samples passed the full run
try:
    full_summary_csv_path = list(glob('%s/%s/*_quality_control/*-summary.csv' % (argv[1], argv[2])))[0]
except:
    full_summary_csv_path = None
if full_summary_csv_path is None:
    passed_IDs = None
else:
    passed_IDs = set()
    flag_inds_true = set()
    flag_inds_false = set()
    fields = dict() # fields[name of field] = index
    for row in csv_reader(open(full_summary_csv_path)):
        # first row, so find flag columns
        if len(fields) == 0:
            for col,name in enumerate(row):
                field_name = name.strip().lower()
                fields[field_name] = col
                if field_name in {'is_accepted'}: # should be TRUE
                    flag_inds_true.add(col)
                elif field_name in {}:#'indels_flagged'}: # should be FALSE
                    flag_inds_false.add(col)

        # not first row, so check if this sample passed all flags
        else:
            passed_all = True
            for col in flag_inds_true:
                if row[col].strip().lower() != 'true':
                    passed_all = False; break
            for col in flag_inds_false:
                if row[col].strip().lower() != 'false':
                    passed_all = False; break
            if passed_all:
                passed_IDs.add(row[fields['sample']])

# measure the deltas
x = list() # total number of mapped reads
y = list() # proportion of passed_IDs that yielded identical consensus sequence in non-ignored range
for run_folder in glob('%s/*' % argv[1]):
    # ignore non-folders and "full" folder
    if not isdir(run_folder) or run_folder.split('/')[-1] == argv[2] or run_folder.split('/')[-1].lower().startswith('full_'):
        continue

    # compute the "total number of successfully mapped reads" limit of this run
    try:
        n = max(int(float(row[fields['mapped_reads']])) for row in csv_reader(open(list(glob('%s/*_quality_control/*-summary.csv' % run_folder))[0])) if row[fields['mapped_reads']].lower().strip() not in {'','mapped_reads'})
        summary_present = True
    except:
        n = 0.; summary_present = False

    # check each sample
    identical = 0.
    for sample_folder in glob('%s/*' % run_folder): # TODO this changed across runs
        sample_ID = sample_folder.split('/')[-1].strip()
        if passed_IDs is None or sample_ID in passed_IDs:
            # load subsampled and full alignments
            sub_alignment = json_load(open(list(glob('%s/*.align.json' % sample_folder))[0]))
            full_alignment = json_load(open(list(glob('%s/%s/%s*/*.align.json' % (argv[1], argv[2], sample_ID)))[0])) # TODO this changed across runs

            # extract subsampled and full consensus sequences
            sub_consensus = ''.join(l.strip() for l in open(list(glob('%s/*.consensus.fa' % sample_folder))[0]).read().splitlines() if not l.startswith('>'))
            full_consensus = ''.join(l.strip() for l in open(list(glob('%s/%s/%s*/*.consensus.fa' % (argv[1], argv[2], sample_ID)))[0]).read().splitlines() if not l.startswith('>')) # TODO this changed across runs

            # get start indices (inclusive) of non-ignored windows
            sub_substring_start = sub_alignment['cons_first_orf_start_0based'] - ORF_START_TO_SUBSTRING_START
            full_substring_start = full_alignment['cons_first_orf_start_0based'] - ORF_START_TO_SUBSTRING_START

            # get end indices (exclusive) of non-ignored windows
            sub_substring_end = sub_alignment['cons_last_orf_end_0based'] + ORF_END_TO_SUBSTRING_END
            full_substring_end = full_alignment['cons_last_orf_end_0based'] + ORF_END_TO_SUBSTRING_END

            # get trimmed substrings
            sub_substring = sub_consensus[sub_substring_start : sub_substring_end]
            full_substring = full_consensus[full_substring_start : full_substring_end]

            # if identical, success
            if sub_substring == full_substring:
                identical += 1
            else: # debug
                pass#print(sub_substring); print(full_substring) # debug

            # if no summary.csv is present, keep track of the number of samples
            if not summary_present:
                n += 1

    # add this point to the output
    if summary_present:
        x.append(n); y.append(identical/len(passed_IDs))
    else:
        x.append(run_folder); y.append(identical/n)

# plot scatterplot
plt.scatter(x, y)
plt.title("Subsampling Experiment")
plt.xlabel("Maximum Allowed Number of Mapped Reads")
plt.ylabel("Prop. Samples w/ Identical Consensus\n(ignore first %d and last %d positions of reference)" % (IGNORE_FIRST, IGNORE_LAST))
plt.show()
