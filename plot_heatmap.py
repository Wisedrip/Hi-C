import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import multiprocess as mp
import bioframe
import cooler
import itertools
import click
import cooltools
#import cooltools.eigdecomp
#import cooltools.expected
#import cooltools.saddle
from dask.distributed import Client, LocalCluster
from scipy.linalg import toeplitz

folder_path = '/data/zuowu/CTCFL/Hi-C/hepatocytes/merge/1Mb'

for coolfile in os.listdir(folder_path):
    if coolfile.endswith('.cool'):
        prefix = os.path.splitext(coolfile)[0]
        print(f'Processing {prefix}')
        c = cooler.Cooler(coolfile)
        obs_mat = c.matrix()[:]
        scale='linear'
        out=''.join([prefix,'_all_by_all_linear_yor_obs.pdf'])
        dpi= 100
        colormap='YlOrRd'
        row_matrix= prefix
        col_matrix= prefix
        zmin=0.00000
        zmax=0.0004
        plt.figure(figsize=(10,10))
        plt.gcf().canvas.manager.set_window_title("Contact matrix".format())
        plt.title("")
        plt.imshow(obs_mat, interpolation="none",vmin=zmin,vmax=zmax, cmap=colormap)
        plt.ylabel("{} coordinate".format(row_matrix))
        plt.xlabel("{} coordinate".format(col_matrix))
        cb = plt.colorbar()
        cb.set_label({"linear": "relative contact frequency", "log2": "log 2 ( relative contact frequency )",
          "log10": "log 10 ( relative contact frequency )",
          }[scale])
        plt.savefig(out, dpi=dpi, format='pdf')

# cis region heatmap"
chroms=['chr1','chr2']
for coolfile in os.listdir(folder_path):
    if coolfile.endswith('.cool'):
        prefix = os.path.splitext(coolfile)[0]
        print(f'Processing {prefix}')
        c = cooler.Cooler(coolfile)
        for chr in chroms:
            mat = c.matrix().fetch(chr)
            row_chrom= chr
            col_chrom= chr
            scale='linear'
            out=''.join([chr,'_', prefix,'_obs_1Mb.pdf'])
            dpi= 300
            colormap='YlOrRd'
            zmin= 0
            zmax= 0.004
            plt.figure(figsize=(10,10))
            plt.gcf().canvas.manager.set_window_title("Contact matrix".format())
            plt.title("")
            plt.imshow(mat, interpolation="none",vmin=zmin,vmax=zmax, cmap=colormap)
            plt.ylabel("{} coordinate".format(row_chrom))
            plt.xlabel("{} coordinate".format(col_chrom))
            cb = plt.colorbar()
            cb.set_label({"linear": "relative contact frequency",
            "log2": "log 2 ( relative contact frequency )",
            "log10": "log 10 ( relative contact frequency )",
             }[scale])
            plt.savefig(out, dpi=dpi, format='pdf') 

chroms=['chr1','chr2']
for coolfile in os.listdir(folder_path):
    if coolfile.endswith('.cool'):
        prefix = os.path.splitext(coolfile)[0]
        print(f'Processing {prefix}')
        c = cooler.Cooler(coolfile)
        for chr in chroms:
            mat = c.matrix().fetch(chr)
            row_chrom= chr
            col_chrom= chr
            scale='log2'
            out=''.join([chr,'_', prefix,'_log2_obs_1Mb.pdf'])
            dpi= 300
            colormap='YlOrRd'
            zmin= -15
            zmax= -6
            plt.figure(figsize=(10,10))
            plt.gcf().canvas.manager.set_window_title("Contact matrix".format())
            plt.title("")
            plt.imshow(np.log2(mat), interpolation="none",vmin=zmin,vmax=zmax, cmap=colormap)
            plt.ylabel("{} coordinate".format(row_chrom))
            plt.xlabel("{} coordinate".format(col_chrom))
            cb = plt.colorbar()
            cb.set_label({"linear": "relative contact frequency",
            "log2": "log 2 ( relative contact frequency )",
            "log10": "log 10 ( relative contact frequency )",
             }[scale])
            plt.savefig(out, dpi=dpi, format='pdf') 


# trans region heatmap"
for coolfile in os.listdir(folder_path):
    if coolfile.endswith('.cool'):
        prefix = os.path.splitext(coolfile)[0]
        print(f'Processing {prefix}')
        c = cooler.Cooler(coolfile)
        mat = c.matrix().fetch('chr1','chr2')
        row_chrom='chr1'
        col_chrom='chr2'
        zmin=0
        zmax=0.00035
        out=''.join([prefix,'_chr1_chr2_obs_1Mb_YlOrRd.pdf'])
        dpi= 300
        scale='linear'
        colormap='YlOrRd'
        plt.figure(figsize=(10,10))
        plt.gcf().canvas.manager.set_window_title("Contact matrix".format())
        plt.title("")
        plt.imshow(mat, interpolation="none", vmin=zmin,vmax=zmax,cmap=colormap)
        plt.ylabel("{} coordinate".format(row_chrom))
        plt.xlabel("{} coordinate".format(col_chrom))
        cb = plt.colorbar()
        cb.set_label({"linear": "relative contact frequency", "log2": "log 2 ( relative contact frequency )",
          "log10": "log 10 ( relative contact frequency )",
          }[scale])
        plt.savefig(out, dpi=dpi, format='pdf')
