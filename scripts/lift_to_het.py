
import pandas as pd
import argparse

header_list = ["chr", "bp", "hg38", "monkey", "hg38_strand", "monkey_strand", "chrom", "b_p",
               "monkey_ref", "monkey_alt", "AC", "DP", "MQ", "FS", "QD", "SOR", "GQ", "ADFrac", "ADP",
               "f_monkey_ref", "f_monkey_alt"]

parser = argparse.ArgumentParser()

parser.add_argument('-f', help='file', type=str)
parser.add_argument('-i', help='path to infile', type=str)
parser.add_argument('-o', help='outpath', type=str)
args = parser.parse_args()

name = args.f
data_path = args.i
out_path = args.o

name_elements = name.split(".")
ind = name_elements[0]
species = name_elements[1][:-10]
match_type = name_elements[4]
data = pd.read_csv(data_path+name,
                compression="gzip", sep="\t", names=header_list)
subset_list = []
data_dict = {}
subset_list.append((data.loc[(data.chr != "chrX") & (data.chr != "chrY")], "autosome"))
subset_list.append((data.loc[(data.chr == "chrX") & (data.bp <= 2700000)], "PAR"))
subset_list.append((data.loc[(data.chr == "chrX") & (data.bp > 2700000)], "nonPAR"))
subset_list.append((data.loc[(data.chr == "chrY")], "chrY"))
data_dict["PGDP_ID"], data_dict["species"], data_dict["match_type"] = ind, species, match_type
for subset, n in subset_list:
    data_dict["lines_"+n] = len(subset)
    data_dict["GQ_"+n] = len(subset.loc[subset.GQ == 99])
    data_dict["ADFrac_"+n] = len(subset.loc[(subset.ADFrac < 0.7) & (subset.ADFrac > 0.3)])
    data_dict["both_"+n] = len(subset.loc[(subset.ADFrac < 0.7) & (subset.ADFrac > 0.3) & (subset.GQ == 99)])
df = pd.DataFrame(data_dict, index=[0])
df.to_csv(out_path+"{}.{}.txt".format(ind, match_type), index=False)
