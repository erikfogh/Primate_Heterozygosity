
import pandas as pd
import argparse

header_list = ["chr", "bp", "hg38", "monkey", "hg38_strand", "monkey_strand", "chrom", "b_p",
               "monkey_ref", "monkey_alt", "AC", "DP", "MQ", "FS", "QD", "SOR", "GQ", "ADFrac", "ADP",
               "f_monkey_ref", "f_monkey_alt"]

parser = argparse.ArgumentParser()

parser.add_argument('-f', help='file', type=str)
parser.add_argument('-i', help='path to infile', type=str)
parser.add_argument('-o', help='outpath', type=str)
parser.add_argument('-w', help='window size', type=int)
args = parser.parse_args()

name = args.f
data_path = args.i
out_path = args.o
mask = "/home/eriks/primatediversity/people/erik/data/conservation_beds/{}.bed"
window_size = args.w

name_elements = name.split(".")
ind = name_elements[0]
species = name_elements[1][:-10]
match_type = name_elements[4]
data_full = pd.read_csv(data_path+name,
                compression="gzip", sep="\t", names=header_list)
df_list = []


for chrom in data_full.chr.unique():
    data_dict = {}
    data_dict["lines"] = []
    data_dict["GQ_pass"] = []
    data_dict["ADFrac_pass"] = []
    data_dict["both_pass"] = []
    data_dict["mask_percentage"] = []
    window_l = []
    chrom_data = data_full.loc[data_full.chr == chrom]
    chrom_mask = pd.read_csv(mask.format(chrom), sep='\t', comment='t', header=None)
    mask_header = ['chrom', 'start', 'end', 'name']
    chrom_mask.columns = mask_header
    print(chrom)
    for i in range(0, chrom_data.bp.iloc[-1], window_size):
        subset = chrom_data.loc[(chrom_data.bp >= i) & (chrom_data.bp < i+window_size)]
        mask_subset = chrom_mask.loc[(chrom_mask.end >= i) & (chrom_mask.start < i+window_size)]
        start_sum, end_sum = sum(mask_subset.start), sum(mask_subset.end)
        if mask_subset.start[0] < i:
            start_sum-mask_subset.start[0]+i
        if mask_subset.end[-1] > i+window_size:
            end_sum-mask_subset.end[-1]+i+window_size
        data_dict["mask_percentage"].append((end_sum-start_sum)/window_size)
        
        data_dict["lines"].append(len(subset))
        data_dict["GQ_pass"].append(len(subset.loc[subset.GQ == 99]))
        data_dict["ADFrac_pass"].append(len(subset.loc[(subset.ADFrac < 0.7) & (subset.ADFrac > 0.3)]))
        data_dict["both_pass"].append(len(subset.loc[(subset.ADFrac < 0.7) & (subset.ADFrac > 0.3) & (subset.GQ == 99)]))
        window_l.append(i)
    df = pd.DataFrame(data_dict)
    df["PGDP_ID"], df["species"], df["match_type"] = ind, species, match_type
    df["chr"], df["window"] = chrom, window_l
    df_list.append(df)

df_out = pd.concat(df_list)
df_out.to_csv(out_path+"{}.{}.window{}.txt".format(ind, match_type, window_size), index=False)
