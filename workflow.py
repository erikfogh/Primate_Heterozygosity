from gwf import Workflow, AnonymousTarget
import os
import glob
import pandas as pd
import genominterv
from groups import Group

gwf = Workflow(defaults={"account": "primatediversity"})

#Inputs
vcf_dir = "/faststorage/project/primatediversity/data/variants/"
vcf_suffix = ".variable.filtered.HF.snps.vcf.gz"
metainfo = "data/New_Papio.xlsx"
chromosomes = list(range(1, 21))+["X"]
meta_df = pd.read_excel(metainfo)

def vcf_tools(chrom, vcf, out):
    inputs = vcf
    outputs = out+"/chrom{}".format(chrom)
    options = {
        "cores": 2,
        "memory": "8g",
        "walltime": "1:00:00"
    }
    spec = """
    vcftools --gzvcf {} --out {} --chr chr{} --counts2 --min-alleles 2 --max-alleles 2
    """.format(vcf, outputs, chrom)
    print(spec)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)




for i, row in meta_df.iterrows():
    ID = row.PGDP_ID
    sex = row.Sex
    vcf = vcf_dir+ID+vcf_suffix
    out = "steps/"+ID
    os.makedirs(out, exist_ok=True)
    with Group(gwf, suffix=ID) as g:
        if ID == "PD_0793":
            print(ID, sex, vcf)
            chr_gen = g.map(vcf_tools, chromosomes, name="vcf_tools", extra={
                "vcf": vcf, "out": out
            })