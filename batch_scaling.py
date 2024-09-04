import warnings
warnings.filterwarnings("ignore")
from itertools import combinations
# import semi-core packages
import matplotlib.pyplot as plt
from matplotlib import colors
plt.style.use('seaborn-poster')
import numpy as np
import pandas as pd
# import open2c libraries
import bioframe
import cooler
import cooltools
import os

mm10_arms = pd.read_csv('mm10_chr_arms.bed', sep='\t')
mm10_arms=mm10_arms[0:20]

folder_path = '/data/zuowu/CTCFL/Hi-C/hepatocytes/merge/'

for filename in os.listdir(folder_path):
    if filename.endswith('.cool'):
        file_path = os.path.join(folder_path, filename)
        file_prefix = os.path.splitext(filename)[0]
        print(f'Processing {file_prefix}')
        coolfile = cooler.Cooler(file_path)
        coolfile_exp = cooltools.expected_cis(
            clr=coolfile,
            view_df=mm10_arms,
            smooth=True,
            aggregate_smoothed=True,
            nproc=12 
        )
        coolfile_exp.to_csv(f'{file_prefix}_expected_smoothed.tsv.csv' ,sep='\t', index=False, na_rep='nan')
