"""
Port the VOC searching to use readiff output, 
Make this portable to pypy so use cffi for readiff and pandas only for the voc searching

NOTE: Using the cffi with pypy may require some tweaking of pythons garbage collection
TODO: replace pandas calls to utilize pypy
TODO: Verify that for indels only exact matches produce a positive result
TODO: Add in MNP support
TODO: As testing if item is in list, convert lists to sets
TODO: Maybe more efficient to implement regex seaerching instead of set matching

2022-02-08: Matthew Wells
"""
import enum
from re import T
import numpy as np
#import pandas as pd
import time
import functools
import sys
import copy
import cProfile
import sys

##---GLOBAL---
_VOC_DICT_ = dict()


def timer(func):
    """Print the runtime of the decorated function
        from: https://realpython.com/primer-on-python-decorators/#timing-functions
    """
    @functools.wraps(func)
    def wrapper_timer(*args, **kwargs):
        start_time = time.perf_counter()    # 1
        value = func(*args, **kwargs)
        end_time = time.perf_counter()      # 2
        run_time = end_time - start_time    # 3
        print(f"Finished {func.__name__!r} in {run_time:.4f} secs")
        return value
    return wrapper_timer



#def read_voc_table(fp_tsv, sep="\t"):
#    """
#    Parse the tsv voc file, this will allow for pypy to be utilized. resulting in speed ups
#    of the loops used below for the variant tabulation
#    """
#    columns  = []
#    mutation_types = ["Sub", "Ins", "Del"]
#    vocs = {}
#    with open(fp_tsv, 'r', encoding='utf-8') as voc_file:
#        header = next(voc_file)
#        columns = header.strip("\n").split(sep)
#        group_var = columns[0] # VOC is the gathering val
#        for row in voc_file:
#            row_ = row.strip("\n").split(sep)
#            group = row_[0] # the grouping value
#            if vocs.get(group) is None:
#                vocs[group]

        


#df = pd.read_csv(VOCTABLE, sep="\t")
cols = ['VOC', 'PangoLineage', 'NextStrainClade', 'NucName', 'AAName', 'Key',
       'SignatureSNV', 'Position', 'Type', 'Length', 'Ref', 'Alt']

grouping_variable = cols[0]
mutation_types = ["Sub", "Ins", "Del"]
# Example deletion
# B.1.243	B.1.243	B.1.243	[22296-22304]del	S:246-248del	TRUE	TRUE	22295	Del	9	CATAGAAGTT	C
# Example insertion
# Panama	A.2.5	TBD	22205CGGCAGGCT	S:ins215AAG	TRUE	TRUE	22205	Ins	9	G	GCGGCAGGCT


# commenting out while testing other versions
#groups = df[grouping_variable].unique()
#print("Preparing mutation strings for the following groups", flush=True, file=sys.stderr)
#groups_lvl = len(groups)
#for i in range(groups_lvl):
#    print(f"\t{i + 1}) {groups[i]}", flush=True, file=sys.stderr)
#
#voc_mut_dict = dict()
#print("Preparing search strings", flush=True, file=sys.stderr)
#for group in groups:
#    temp_df = df.loc[df[grouping_variable] == group]
#    voc_mut_dict[group] = []
#    sub_only_df = temp_df.loc[temp_df[cols[8]] == mutation_types[0]]
#    ins_only_df = temp_df.loc[temp_df[cols[8]] == mutation_types[1]]
#    del_only_df = temp_df.loc[temp_df[cols[8]] == mutation_types[2]]
#
#    for i in sub_only_df.index: # leaving loop in preperation of pandas removal
#        # converting combos to keys would yield a speed boost as well
#        voc_mut_dict[group].append(f"{sub_only_df.at[i, cols[10]]}{int(sub_only_df.at[i, cols[7]])}{sub_only_df.at[i, cols[-1]]}")  
#
#    for i in ins_only_df.index:
#        pos = int(ins_only_df.at[i, cols[7]]) + 1
#        alt = ins_only_df.at[i, cols[-1]][1:]
#        ins_val = "".join([f"-{pos + k}{v}" for k, v in enumerate(alt)])
#        voc_mut_dict[group].append(ins_val)
#
#    for i in del_only_df.index:
#        # indexing of values reflects the changes made in the vcfpaser sheet
#        pos = int(del_only_df.at[i, cols[7]]) + 1
#        ref = del_only_df.at[i, cols[10]][1:] # skip reference base
#        #alt = del_only_df.at[i, cols[-1]] alt is empty
#        del_value = ''.join([f"{v}{pos + k}-" for k, v in enumerate(ref)])
#        voc_mut_dict[group].append(del_value)
#    voc_mut_dict[group] = set(voc_mut_dict[group])
#

cols = ['VOC', 'PangoLineage', 'NextStrainClade', 'NucName', 'AAName', 'Key',
       'SignatureSNV', 'Position', 'Type', 'Length', 'Ref', 'Alt']

def inssubdel_preperation(mut_type, pos, ref, alt):
    """
    The previous implementation used pandas to parse the dataframe, but pandas
    support in pypy is fickle. So this func is to create the indels as was done in the 
    pandas version

    """
    mut_type = mut_type.lower()
    mutation_types = ["sub", "ins", "del"]
    mutations = []
    if mut_type == mutation_types[0]: # snp
        mutations.append(f"{ref}{pos}{alt}")
    elif mut_type == mutation_types[1]: # Insertion
        pos = int(pos) + 1 # position given is 0 based
        alt = alt[1:] # as reference base is included in mut string
        mutation_types.append(''.join([f"-{pos + k}{v}" for k, v in enumerate(alt)]))
    elif mut_type == mutation_types[2]: # deletion
        pos = int(pos) + 1 # once again dels are given 0 based as well
        ref = ref[1:] # skip reference base in string
        mutation_types.append(''.join([f"{v}{pos + k}-" for k, v in enumerate(ref)]))
    else:
        print(f"Mutation type of {mut_type} is currently not supported.", flush=True,
        file=sys.stderr)
    return mutations



def read_voc_table(fp_tsv, sep="\t"):
    """
    Parse the tsv voc file, this will allow for pypy to be utilized. resulting in speed ups
    of the loops used below for the variant tabulation
    """
    columns  = []
    vocs = {} # may set this globally to help with logging
    with open(fp_tsv, 'r', encoding='utf-8') as voc_file:
        header = next(voc_file)
        columns = header.strip("\n").split(sep)
        for row in voc_file:
            row_ = row.strip("\n").split(sep)
            group = row_[0] # the grouping value
            pos = row_[7] # position
            mut_type = row_[-4] # mutation type
            ref = row_[-2]
            alt = row_[-1]
            if vocs.get(group) is None:
                vocs[group] = inssubdel_preperation(mut_type, pos, ref, alt)
            else:
                vocs[group].extend(inssubdel_preperation(mut_type, pos, ref, alt))
        #print(vocs[group])
        vocs[group] = set(vocs[group]) # faster then frozen set
    return vocs



@timer
def tabulate_combos(r_diff, voc_dict):
    for combo in r_diff.index:
        for key in voc_dict.keys(): # preeration of voc mut combos based on proximity will yield a speed boost
            snps = r_diff.iat[combo, 1].split("|") # second columm contains the combination
            mut_count = 0
            muts = voc_dict[key]
            for snp in snps:
                if snp in muts:
                    mut_count += 1
            out_df.at[key, mut_count] += 1



def process_indels(rd_list: list):
    """
    Take the input from readdiff and reformat indels
    TODO: evaluate performance of assinging rd_list[rd_idx] to stat valus
    TODO: Convert contiguous changes into an mnp, maybe better to do this globally
    :param rd_str: The read diff string split on '|'
    """

    # Mulitple lists is a major bottle neck, this is just for debug purposes
    ins_strings = []
    del_strings = []
    snp_strings = []
    mnp_strings = []
    len_rd = len(rd_list)
    rd_idx = 0
    while rd_idx < len_rd:
        if rd_list[rd_idx].startswith("-"):
            #pos = int("".join([x for x in rd_list[rd_idx] if x.isdigit()])) # can speed this up with a slice index is known
            pos = int(rd_list[rd_idx][1:-1])
            ins_str = [f"-{pos}{rd_list[rd_idx][-1]}"]
            tabulator = 1
            rd_idx += 1
            while True and rd_idx < len_rd:
                #pos_ = int("".join([x for x in rd_list[rd_idx] if x.isdigit()]))
                pos_ = int(rd_list[rd_idx][1:-1])
                comp = pos_ - tabulator
                if rd_list[rd_idx].startswith("-") and pos == comp:
                    ins_str.append(f"-{pos_}{rd_list[rd_idx][-1]}")
                    tabulator += 1
                    rd_idx += 1
                else:
                    ins_strings.append("".join(ins_str))
                    #rd_idx -= 1
                    break
        elif rd_list[rd_idx].endswith("-"):
            #pos = int("".join([x for x in rd_list[rd_idx] if x.isdigit()]))
            pos = int(rd_list[rd_idx][1:-1])
            del_str = [f"{rd_list[rd_idx][0]}{pos}-"]
            tabulator = 1
            rd_idx += 1
            while True and rd_idx < len_rd:
                #pos_ = int("".join([x for x in rd_list[rd_idx] if x.isdigit()]))
                pos_ = int(rd_list[rd_idx][1:-1])
                comp = pos_ - tabulator
                if rd_list[rd_idx].endswith("-") and pos == comp:
                    del_str.append(f"{rd_list[rd_idx][0]}{pos_}-")
                    tabulator += 1
                    rd_idx += 1
                else:
                    del_strings.append("".join(del_str))
                    #rd_idx -= 1
                    break
        else:
            snp_strings.append(rd_list[rd_idx])
            rd_idx += 1
    
    # probably faster to do this in the loop, but the snp strings can be screened for contigous
    # nucleotides out of the loop so doing it here now
    mnp_start = ""
    mnp_end = ""
    mnp_pos = 0
    remove_mnp_l = []
    k = 0 # k is the tracker fro the list
    len_snps = len(snp_strings)
    if len_snps >= 2: # 2 or more values needed to call an mnp
        while k < (len_snps - 1): # prevent breach of index, as there is the +1 step
            pos = int(snp_strings[k][1:-1])
            npos = int(snp_strings[k + 1][1:-1])
            if pos == (npos - 1):
                mnp_pos = pos
                remove_mnp_l.append(k)
                remove_mnp_l.append(k + 1)
                mnp_start += snp_strings[k][0] + snp_strings[k + 1][0] 
                mnp_end += snp_strings[k][-1] + snp_strings[k + 1][-1] 
                tracker = 1
                k += 2 # skip npos and go to the other
                for val in snp_strings[k:]:
                    if (int(val[1:-1]) - tracker) == npos:
                        mnp_start += val[0]
                        mnp_end += val[-1]
                        remove_mnp_l.append(k)
                        tracker += 1
                        k += 1
                    else:
                        mnp_strings.append(f"{mnp_start}{mnp_pos}{mnp_end}")
                        mnp_end = ""
                        mnp_start = ""
                        mnp_pos = 0
                        k += 1
                        break
                        #k -= 1 # TODO think on if this is nesseccary later
            else:
                k += 1
       
    snp_strings.extend(del_strings)
    snp_strings.extend(ins_strings)
    for i in remove_mnp_l:
        snp_strings[i] = None
    #if mnp_strings != []:
    #    print(mnp_strings, snp_strings)
    snp_strings.extend(mnp_strings)
    return set(snp_strings)

@timer
def tabulate_combos_nopd(fp, voc_dict, out_df_):
    """Tabulate mutation combos on file stream may not be issue if using cffi"""
    with open(fp, 'r') as combos:
        next(combos) # skip the first line (header) in the file
        for combo in combos:
            vals = combo.split("\t")
            for key in voc_dict.keys():
                snps = vals[1].strip().split("|")[1:] # strings start with pipe, remove it
                snps = process_indels(snps)
                count_ = int(vals[0])
                mut_count = 0
                muts = voc_dict[key]
                for snp in snps:
                    if snp in muts:
                        mut_count += 1
                out_df_[key][mut_count] += (1 * count_)
    return out_df_




#VOCTABLE = sys.argv[1]
#VOCTABLE = "../data/VCFParser_20220204.txt"
VOCTABLE = "/mnt/w/Projects/Project_Chrystal/2020_SARS-CoV-2_Sewage_Project/Bioinformatics/ReportAutomation/Dojo/ReportablesDataGeneration/VCFParser_20220204_new_yaml_testing.txt"

_VOC_DICT_ = read_voc_table(VOCTABLE)

print("Prepared mutation strings for the following groups", flush=True, file=sys.stderr)
groups = [i for i in _VOC_DICT_.keys()]
groups_lvl = len(groups)
for i in range(groups_lvl):
    print(f"\t{i + 1}) {groups[i]}", flush=True, file=sys.stderr)


# arbitrary value for column size as just starting
#out_df = pd.DataFrame(0, index=_VOC_DICT_.keys(), columns=range(12))
out_df = {}

# default dict can probably do this faster
out_df_cols = {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0}
for key in _VOC_DICT_.keys():
    #out_df[key] = copy.deepcopy(out_df_cols) # profile time of deep copy look for alternative
    out_df[key] = {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0}

print("Tabulating co-occuring mutations", flush=True, file=sys.stderr)
#file_to_analyse = sys.argv[1]
#filled_df = tabulate_combos_nopd(file_to_analyse, _VOC_DICT_, out_df)
filled_df = tabulate_combos_nopd("data/22_WPG9_NE_0301_1.rdf.tsv", _VOC_DICT_, out_df)
#filled_df = tabulate_combos_nopd('../data/22_WAT2_Leslie_0110_test_readdiff.tsv', _VOC_DICT_, out_df)
#cProfile.run("'../data/22_WAT2_Leslie_0110_test_readdiff.tsv', _VOC_DICT_, out_df")

# print this nicely later

def dict_to_table(out_dict, df_cols):
    """
    Print a 2d dict as a table
    """
    
    longest_key = 0 
    for key in out_dict.keys():
        key_len = len(key)
        if longest_key < key_len:
            longest_key = key_len
    longest_key = longest_key + 2 # add in some white space to format things
    columns = [""] # adding dummy value
    columns.extend([str(i) for i in df_cols.keys()])
    print("{: <40} {: >6} {: >6} {: >6} {: >6} {: >6} {: >6} {: >6} {: >6} {: >6} {: >6} {: >6} {: >6}".format(*columns))
    
    for key in out_dict.keys():
        whitespace_add = longest_key - len(key)
        row_white = ' ' * whitespace_add
        row_vals = []
        for val in out_dict[key].keys():
            row_vals.append(str(out_dict[key][val]))
        row_print = [key]
        row_print.extend(row_vals)
        print("{: <40} {: >6} {: >6} {: >6} {: >6} {: >6} {: >6} {: >6} {: >6} {: >6} {: >6} {: >6} {: >6}".format(*row_print))
        #print(key + row_white + "    ".join(row_vals))

dict_to_table(filled_df, out_df_cols)


def helper_deep():
    out_df_cols = {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0}
    for key in _VOC_DICT_.keys():
        out_df[key] = copy.deepcopy(out_df_cols) # profile time of deep copy look for alternative
        #out_df[key] = {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0}

def helper_assign():
    out_df_cols = {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0}
    for key in _VOC_DICT_.keys():
        #out_df[key] = copy.deepcopy(out_df_cols) # profile time of deep copy look for alternative
        out_df[key] = {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0}

helper_assign()
#cProfile.run("helper_deep()")
#cProfile.run("helper_assign()")
