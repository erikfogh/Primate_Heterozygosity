'''
------------------------------------------------------------------------------------------------------------------------
This workflow is for the baboon data.
------------------------------------------------------------------------------------------------------------------------

------------------------------------------------------------------------------------------------------------------------
Original Author: Juraj Bergman
Date: 02/02/2021
------------------------------------------------------------------------------------------------------------------------
Modifications by Erik Fogh Sørensen
Date: 19/03/2021
------------------------------------------------------------------------------------------------------------------------
'''

from gwf import Workflow
import os
import pandas as pd

# Absolute path to the input data (gvcfs)
gvcfs_path = "/faststorage/project/primatediversity/data/gVCFs_baboons_16_03_2021/"
outpath = "data/callmasks/"  # Making a relate path right now, should probably be altered.
meta_data_samples = pd.read_table("data/metadata_with_x_missing.txt", sep=" ")

gwf = Workflow()


########################################################################################################################
############################################# ---- CALLABILITY MASK ---- ###############################################
########################################################################################################################


def callMask(infile, outfile, min_het, gq, path, stat_path):
    """Individual callability."""
    inputs = [path + infile]
    outputs = [outfile + ".bed"]
    options = {'cores': 2, 'memory': "16g", 'walltime': "06:00:00", "account": 'primatediversity'}

    spec = """
    bcftools stats -d  2,500,1 {gvcf} | grep 'DP' | grep -iv -e '#' -e '<' -e '>' | sort -k 6 -V -r | head -1 | awk {awk} > {stat_path}_modcov.txt

    modcov=$(<{stat_path}_modcov.txt)
    min_cov=$((modcov/2))
    max_cov=$((modcov*2))
    echo $max_cov
    echo $min_cov

    bcftools filter -e "(GT='./.') | (GT='het' & FMT/AD[*:*] < {MIN_HET_AD} ) | FMT/DP <= $min_cov | FMT/DP >= $max_cov | FMT/GQ <= {GQ} " {gvcf} |\
    grep -v '#' | \
    awk 'BEGIN{{OFS="\\t"}}{{ print $1, $2-1, $2 }}' - | \
    bedtools merge | \
    sort -k1,1 -k2,2n | \
    bedtools merge > {outfile}_temp

    mv {outfile}_temp {outfile}.bed
    """.format(awk="'{print $3}'", MIN_HET_AD=min_het, GQ=gq,
               gvcf=path + infile, outfile=outfile, stat_path=stat_path)
    return (inputs, outputs, options, spec)


def callMask_spec_filter_third(infile, outfile, path, stat_path, f):
    """Individual callability."""
    inputs = [path + infile]
    outputs = [outfile + ".bed"]
    options = {'cores': 2, 'memory': "16g", 'walltime': "06:00:00", "account": 'primatediversity'}
    spec = """
    bcftools stats -d  2,500,1 {gvcf} | grep 'DP' | grep -iv -e '#' -e '<' -e '>' | sort -k 6 -V -r | head -1 | awk {awk} > {stat_path}_modcov.txt

    modcov=$(<{stat_path}_modcov.txt)
    min_cov=$((modcov/3))
    max_cov=$((modcov*2))
    echo $max_cov
    echo $min_cov

    {f} {gvcf}  | \
    grep -v '#' | \
    awk 'BEGIN{{OFS="\\t"}}{{ print $1, $2-1, $2 }}' - | \
    bedtools merge | \
    sort -k1,1 -k2,2n | \
    bedtools merge > {outfile}_temp

    mv {outfile}_temp {outfile}.bed
    """.format(awk="'{print $3}'", f=f,
               gvcf=path + infile, outfile=outfile, stat_path=stat_path)
    return (inputs, outputs, options, spec)


########################################################################################################################
############################################# ---- MERGING BED FILES ---- ##############################################
########################################################################################################################

def multiinter_and_trim(input_list, outfile, p):
    """Merging bed files to create one"""
    # print(path + infile)
    inputs = input_list
    outputs = [outfile+"_merged.bed"]
    if p == 0:
        cutoff = len(input_list)
    else:
        c = 100/p
        cutoff = len(input_list)-len(input_list)//c
    options = {'cores': 2, 'memory': "16g", 'walltime': "02:00:00", "account": 'primatediversity'}

    spec = """
    multiIntersectBed -i {i} > {o}

    awk 'BEGIN{{OFS="\\t"}} $4 >= {cutoff} {{print $1, $2, $3}}' {o}  > {o}.bed

    bedtools merge -i {o}.bed > {o}_merged.bed

    rm {o}

    """.format(i=" ".join(input_list), o=outfile, cutoff=cutoff)
    return (inputs, outputs, options, spec)


########################################################################################################################
################################################ ---- RUN PIPELINE ---- ################################################
########################################################################################################################


### get individual callability with Lukas protocol


os.makedirs("steps/depth_stats", exist_ok=True)
os.makedirs(outpath, exist_ok=True)

# Finding the names

names = []
for r, subdir, files in os.walk(gvcfs_path):
    for sdir in subdir:
        names.append(sdir)  # Gives it without order, but order is not needed.
chromosomes = ['chr{}'.format(x) for x in range(1, 21)] + ['chrX']

# Generating a df containing the species relationship for each ID
d = {}
s_list = []
f_list = []
all_individuals = []
for s in meta_data_samples.Species.unique():
    s_meta = meta_data_samples.loc[meta_data_samples.Species == s]
    if s == "gelada":
        continue
    all_list = []
    for i, row in s_meta.iterrows():
        if row.PGDP_ID not in names:
            continue  # skip hypothetical individuals which are removed.
        if row.PGDP_ID.startswith("PD"):   
            all_list.append(row.PGDP_ID)
            all_individuals.append(row.PGDP_ID)
            if row.Sex == "F":
                f_list.append(row.PGDP_ID)
    if s == 'ursinus (grayfoot)':
        s = 'ursinus'
    s_list.append(s)
    d[s] = all_list

for i in range(len(names)):
    os.makedirs(outpath + names[i], exist_ok=True)
    for chrom in chromosomes:
        infile = names[i] + "/output.raw.snps.indels.{}.genotyped.g.vcf.gz".format(chrom)
        if not os.path.exists(gvcfs_path+infile):
            # print("Missing ", infile)
            infile = names[i] + "/output.raw.snps.indels.{}.exclude.genotyped.g.vcf.gz".format(chrom)
        gwf.target_from_template('callMask_' + names[i]+'_'+chrom,
                                 callMask(infile=infile,
                                          outfile=outpath + names[i]+"/{}_mask".format(chrom),
                                          min_het="3",
                                          gq="30",
                                          path=gvcfs_path,
                                          stat_path="steps/depth_stats/{}_{}".format(names[i], chrom)))

f = """bcftools filter -e "(GT='./.') | (GT='het' & FMT/AD[*:*] < 3 ) | FMT/DP <= $min_cov | FMT/DP >= $max_cov | FMT/GQ <= 30 " """
for i in range(len(names)):
    os.makedirs(outpath + names[i], exist_ok=True)
    for chrom in chromosomes:
        infile = names[i] + "/output.raw.snps.indels.{}.genotyped.g.vcf.gz".format(chrom)
        if not os.path.exists(gvcfs_path+infile):
            # print("Missing ", infile)
            infile = names[i] + "/output.raw.snps.indels.{}.exclude.genotyped.g.vcf.gz".format(chrom)
        gwf.target_from_template('callMask_third_' + names[i]+'_'+chrom,
                                 callMask_spec_filter_third(infile=infile,
                                          outfile=outpath + names[i]+"/{}_mask_third".format(chrom),
                                          f=f,
                                          path=gvcfs_path,
                                          stat_path="steps/depth_stats/{}_{}".format(names[i], chrom)))

percentage_list = [50, 20, 15, 10, 5, 2.5, 1, 0.5, 0]

for p in percentage_list:
    for chrom in chromosomes:
        input_list = [outpath+'{}/{}_mask.bed'.format(x, chrom) for x in all_individuals]
        gwf.target_from_template("callMask_all_"+chrom+"_cutoff_"+str(p),
                                 multiinter_and_trim(input_list=input_list,
                                                     p=p,
                                                     outfile=outpath+'cutoff_{}_all_baboons_{}_min_half'.format(p, chrom)
                                 ))

for p in percentage_list:
    for chrom in chromosomes:
        input_list = [outpath+'{}/{}_mask_third.bed'.format(x, chrom) for x in all_individuals]
        gwf.target_from_template("callMask_all_third_"+chrom+"_cutoff_"+str(p),
                                 multiinter_and_trim(input_list=input_list,
                                                     p=p,
                                                     outfile=outpath+'cutoff_{}_all_baboons_{}_min_third'.format(p, chrom)
                                 ))
                        
for p in percentage_list:
    input_list = [outpath+'{}/{}_mask_third.bed'.format(x, "chrX") for x in f_list]
    gwf.target_from_template("callMask_all_third_f_"+chrom+"_cutoff_"+str(p),
                                 multiinter_and_trim(input_list=input_list,
                                                     p=p,
                                                     outfile=outpath+'cutoff_{}_all_baboons_{}_f_only_min_third'.format(p, "chrX")
                                 ))

for p in percentage_list:
    input_list = [outpath+'{}/{}_mask.bed'.format(x, "chrX") for x in f_list]
    gwf.target_from_template("callMask_f_"+chrom+"_cutoff_"+str(p),
                                 multiinter_and_trim(input_list=input_list,
                                                     p=p,
                                                     outfile=outpath+'cutoff_{}_all_baboons_{}_f_only_min_half'.format(p, "chrX")
                                 ))

# Testing the degree of filtering. 1 male and 1 female from each species.


def callMask_spec_filter(infile, outfile, path, stat_path, f):
    """Individual callability."""
    inputs = [path + infile]
    outputs = [outfile + ".bed"]
    options = {'cores': 2, 'memory': "16g", 'walltime': "04:00:00", "account": 'primatediversity'}
    spec = """
    bcftools stats -d  2,500,1 {gvcf} | grep 'DP' | grep -iv -e '#' -e '<' -e '>' | sort -k 6 -V -r | head -1 | awk {awk} > {stat_path}_modcov.txt

    modcov=$(<{stat_path}_modcov.txt)
    min_cov=$((modcov/2))
    max_cov=$((modcov*2))
    echo $max_cov
    echo $min_cov

    {f} {gvcf}  | \
    grep -v '#' | \
    awk 'BEGIN{{OFS="\\t"}}{{ print $1, $2-1, $2 }}' - | \
    bedtools merge | \
    sort -k1,1 -k2,2n | \
    bedtools merge > {outfile}_temp

    mv {outfile}_temp {outfile}.bed
    """.format(awk="'{print $3}'", f=f,
               gvcf=path + infile, outfile=outfile, stat_path=stat_path)
    return (inputs, outputs, options, spec)


names = ['PD_0271', 'PD_0768', 'PD_0222', 'PD_0508', 'PD_0789', 'PD_0793',
         'PD_0708', 'PD_0718', 'PD_0390', 'PD_0401', 'PD_0695', 'PD_0692']

test_name = "only_GT_miss/"
f = """bcftools filter -e "(GT='./.')" """
for i in range(len(names)):
    os.makedirs(outpath + test_name + names[i], exist_ok=True)
    for chrom in chromosomes:
        infile = names[i] + "/output.raw.snps.indels.{}.genotyped.g.vcf.gz".format(chrom)
        if not os.path.exists(gvcfs_path+infile):
            infile = names[i] + "/output.raw.snps.indels.{}.exclude.genotyped.g.vcf.gz".format(chrom)
        gwf.target_from_template('callMask_' + names[i]+'_'+ chrom+ test_name[:-1],
                                 callMask_spec_filter(infile=infile,
                                                      outfile=outpath + test_name + names[i]+"/{}_mask".format(chrom),
                                                      f=f,
                                                      path=gvcfs_path,
                                                      stat_path="steps/depth_stats/{}_{}".format(names[i], chrom)))

test_name = "min_het/"
f = """bcftools filter -e "(GT='het' & FMT/AD[*:*] < 3 )"  """
for i in range(len(names)):
    os.makedirs(outpath + test_name + names[i], exist_ok=True)
    for chrom in chromosomes:
        infile = names[i] + "/output.raw.snps.indels.{}.genotyped.g.vcf.gz".format(chrom)
        if not os.path.exists(gvcfs_path+infile):
            infile = names[i] + "/output.raw.snps.indels.{}.exclude.genotyped.g.vcf.gz".format(chrom)
        gwf.target_from_template('callMask_' + names[i]+'_'+ chrom+ test_name[:-1],
                                 callMask_spec_filter(infile=infile,
                                                      outfile=outpath + test_name + names[i]+"/{}_mask".format(chrom),
                                                      f=f,
                                                      path=gvcfs_path,
                                                      stat_path="steps/depth_stats/{}_{}".format(names[i], chrom)))

test_name = "min_cov/"
f = """bcftools filter -e "FMT/DP <= $min_cov" """
for i in range(len(names)):
    os.makedirs(outpath + test_name + names[i], exist_ok=True)
    for chrom in chromosomes:
        infile = names[i] + "/output.raw.snps.indels.{}.genotyped.g.vcf.gz".format(chrom)
        if not os.path.exists(gvcfs_path+infile):
            infile = names[i] + "/output.raw.snps.indels.{}.exclude.genotyped.g.vcf.gz".format(chrom)
        gwf.target_from_template('callMask_' + names[i]+'_'+ chrom + test_name[:-1],
                                 callMask_spec_filter(infile=infile,
                                                      outfile=outpath + test_name + names[i]+"/{}_mask".format(chrom),
                                                      f=f,
                                                      path=gvcfs_path,
                                                      stat_path="steps/depth_stats/{}_{}".format(names[i], chrom)))

test_name = "max_cov/"
f = """bcftools filter -e "FMT/DP >= $max_cov" """
for i in range(len(names)):
    os.makedirs(outpath + test_name + names[i], exist_ok=True)
    for chrom in chromosomes:
        infile = names[i] + "/output.raw.snps.indels.{}.genotyped.g.vcf.gz".format(chrom)
        if not os.path.exists(gvcfs_path+infile):
            infile = names[i] + "/output.raw.snps.indels.{}.exclude.genotyped.g.vcf.gz".format(chrom)
        gwf.target_from_template('callMask_' + names[i]+'_'+ chrom + test_name[:-1],
                                 callMask_spec_filter(infile=infile,
                                                      outfile=outpath + test_name + names[i]+"/{}_mask".format(chrom),
                                                      f=f,
                                                      path=gvcfs_path,
                                                      stat_path="steps/depth_stats/{}_{}".format(names[i], chrom)))

test_name = "gq/"
f = """bcftools filter -e " FMT/GQ <= 30" """
for i in range(len(names)):
    os.makedirs(outpath + test_name + names[i], exist_ok=True)
    for chrom in chromosomes:
        infile = names[i] + "/output.raw.snps.indels.{}.genotyped.g.vcf.gz".format(chrom)
        if not os.path.exists(gvcfs_path+infile):
            infile = names[i] + "/output.raw.snps.indels.{}.exclude.genotyped.g.vcf.gz".format(chrom)
        gwf.target_from_template('callMask_' + names[i]+'_'+ chrom + test_name[:-1],
                                 callMask_spec_filter(infile=infile,
                                                      outfile=outpath + test_name + names[i]+"/{}_mask".format(chrom),
                                                      f=f,
                                                      path=gvcfs_path,
                                                      stat_path="steps/depth_stats/{}_{}".format(names[i], chrom)))

test_name = "all_filter/"
f = """bcftools filter -e "(GT='./.') | (GT='het' & FMT/AD[*:*] < 3 ) | FMT/DP <= $min_cov | FMT/DP >= $max_cov | FMT/GQ <= 30 " """
for i in range(len(names)):
    os.makedirs(outpath + test_name + names[i], exist_ok=True)
    for chrom in chromosomes:
        infile = names[i] + "/output.raw.snps.indels.{}.genotyped.g.vcf.gz".format(chrom)
        if not os.path.exists(gvcfs_path+infile):
            infile = names[i] + "/output.raw.snps.indels.{}.exclude.genotyped.g.vcf.gz".format(chrom)
        gwf.target_from_template('callMask_' + names[i]+'_'+ chrom + test_name[:-1],
                                 callMask_spec_filter(infile=infile,
                                                      outfile=outpath + test_name + names[i]+"/{}_mask".format(chrom),
                                                      f=f,
                                                      path=gvcfs_path,
                                                      stat_path="steps/depth_stats/{}_{}".format(names[i], chrom)))


test_name = "min_cov_third/"
f = """bcftools filter -e "FMT/DP <= $min_cov" """
for i in range(len(names)):
    os.makedirs(outpath + test_name + names[i], exist_ok=True)
    for chrom in chromosomes:
        infile = names[i] + "/output.raw.snps.indels.{}.genotyped.g.vcf.gz".format(chrom)
        if not os.path.exists(gvcfs_path+infile):
            infile = names[i] + "/output.raw.snps.indels.{}.exclude.genotyped.g.vcf.gz".format(chrom)
        gwf.target_from_template('callMask_' + names[i]+'_'+ chrom + test_name[:-1],
                                 callMask_spec_filter_third(infile=infile,
                                                      outfile=outpath + test_name + names[i]+"/{}_mask".format(chrom),
                                                      f=f,
                                                      path=gvcfs_path,
                                                      stat_path="steps/depth_stats/{}_{}".format(names[i], chrom)))

test_name = "all_filter_third/"
f = """bcftools filter -e "(GT='./.') | (GT='het' & FMT/AD[*:*] < 3 ) | FMT/DP <= $min_cov | FMT/DP >= $max_cov | FMT/GQ <= 30 " """
for i in range(len(names)):
    os.makedirs(outpath + test_name + names[i], exist_ok=True)
    for chrom in chromosomes:
        infile = names[i] + "/output.raw.snps.indels.{}.genotyped.g.vcf.gz".format(chrom)
        if not os.path.exists(gvcfs_path+infile):
            infile = names[i] + "/output.raw.snps.indels.{}.exclude.genotyped.g.vcf.gz".format(chrom)
        gwf.target_from_template('callMask_' + names[i]+'_'+ chrom + test_name[:-1],
                                 callMask_spec_filter_third(infile=infile,
                                                      outfile=outpath + test_name + names[i]+"/{}_mask".format(chrom),
                                                      f=f,
                                                      path=gvcfs_path,
                                                      stat_path="steps/depth_stats/{}_{}".format(names[i], chrom)))

for s in d:
    for chrom in chromosomes:
        input_list = [outpath+'{}/{}_mask.bed'.format(x, chrom) for x in d[s]]
        gwf.target_from_template("callMask_"+s+"_"+chrom,
                                 multiinter_and_trim(input_list=input_list,
                                                     p=10,
                                                     outfile=outpath+'{}_{}'.format(s, chrom)
                                 ))

# for chrom in chromosomes:
#     input_list = [outpath+'{}_{}_merged.bed'.format(x, chrom) for x in s_list]
#     gwf.target_from_template("callMask_all_"+chrom,
#                              multiinter_and_trim(input_list=input_list,
#                                                  p=10,
#                                                  outfile=outpath+'all_baboons_{}'.format(chrom)
#                              ))
