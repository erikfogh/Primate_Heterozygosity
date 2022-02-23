'''
------------------------------------------------------------------------------------------------------------------------
This workflow is for the primatediversity lifted data.
------------------------------------------------------------------------------------------------------------------------

------------------------------------------------------------------------------------------------------------------------
Modifications by Erik Fogh Sørensen
Date: 05/09/2021
------------------------------------------------------------------------------------------------------------------------
'''

from gwf import Workflow
import os
import pandas as pd

# Absolute path to the input data (gvcfs)
data_path = "/home/eriks/primatediversity/data/lifted_tables_human_coordinates_21_04_21/"
outpath = "steps/het_stats/"
files = os.listdir(data_path)
window = 1000000

gwf = Workflow()


########################################################################################################################
############################################# ---- CALLABILITY MASK ---- ###############################################
########################################################################################################################


def het_count(infile, inpath, outpath):
    """Individual het stats."""
    f_names = infile.split(".")
    inputs = [inpath + infile]
    outputs = [outpath+"{}.{}.txt".format(f_names[0], f_names[4])]
    options = {'cores': 2, 'memory': "24g", 'walltime': "01:00:00", "account": 'primatediversity'}

    spec = """
    python scripts/lift_to_het.py -f {} -i {} -o {}
    """.format(infile, inpath, outpath)
    return (inputs, outputs, options, spec)


def het_count_windows(infile, inpath, outpath, window_size):
    """Individual het stats."""
    f_names = infile.split(".")
    inputs = [inpath + infile]
    outputs = [outpath+"{}.{}.window{}.txt".format(f_names[0], f_names[4], window_size)]
    options = {'cores': 2, 'memory': "24g", 'walltime': "01:00:00", "account": 'primatediversity'}

    spec = """
    python scripts/lift_to_het_windows.py -f {} -i {} -o {} -w {}
    """.format(infile, inpath, outpath, window_size)
    return (inputs, outputs, options, spec)


########################################################################################################################
################################################ ---- RUN PIPELINE ---- ################################################
########################################################################################################################


os.makedirs(outpath+"temp/", exist_ok=True)

l = []
for f in files:
     if f.endswith(".gz"):
         l.append(f)


#gwf.map(het_count, l, extra= {"inpath": data_path, "outpath": outpath+"temp/"})

gwf.map(het_count_windows, l, extra= {"inpath": data_path, "outpath": outpath+"temp/", "window_size": window})