import json
import numpy as np
from glob import glob
import os
import yaml
import pandas as pd
import matplotlib.pyplot as plt

plt.rc('font', family='serif', serif='Times New Roman')
plt.rc('xtick', labelsize='small')
plt.rc('ytick', labelsize='small')
plt.rc('text', usetex=True)
plt.rc('savefig', dpi=500)
plt.rc('xtick', labelsize='small')
plt.rc('ytick', labelsize='small')
plt.rcParams.update({'font.size': 11})
params = {'text.latex.preamble': [r'\usepackage{amsmath}']}
plt.rcParams.update(params)



def get_Ns(dr):
    return int(os.path.basename(os.path.dirname(dr)).split('_')[1])

def get_sample_index(dr):
    return int(os.path.basename(os.path.dirname(dr)).split('_')[2])


def get_dlm_runtime(filename, steps):

    with open(filename) as fh:
        for line in fh.readlines():
            line = line.strip().split()

            if len(line) == 13 and (line[0] == 'Iteration' and int(line[1]) == steps):
                assert line[-1] == 's'
                s = float(line[-2])
                assert line[-4] == 'm'
                m = float(line[-5])
                assert line[-7] == 'h'
                h = float(line[-8])
                return s + 60 * (m + 60 * h)


def get_dlm_accept_count(filename):

    with open(filename) as fh:
        for line in fh.readlines():
            if line.startswith(' successful no of atom moves & acceptance'):
                line = line.strip().split()
                return int(line[-2])



if __name__ == '__main__':
    
    L = 12
    dlm_out = glob('DLM*/OUTPUT*')
    
    with open('input.json') as fh:
        input_data = json.loads(fh.read())

    steps = input_data['steps']

    
    data = {
        'method': [],
        'N' : [],
        'sample_index': [],
        'T' : [],
        'accept_ratio': [],
    }

    for filex in dlm_out:

        N = get_Ns(filex) ** 3
        T = get_dlm_runtime(filex, steps)
        sample_index = get_sample_index(filex)
        naccept = get_dlm_accept_count(filex)
        
        data['method'].append('dlm')
        data['N'].append(N)
        data['T'].append(1000 * T / steps)
        data['sample_index'].append(sample_index)
        data['accept_ratio'].append(naccept/steps)
 

    lmmm = glob('DLM*/timing_data*.json')
    for filex in lmmm:

        sample_index = get_sample_index(filex)
        with open(filex) as fh:
            file_data = json.loads(fh.read())
        
        if file_data['L'] != L: continue

        data['method'].append(file_data['method'])
        data['N'].append(file_data['N'])
        data['T'].append(1000 * file_data['time_taken'] / steps)
        data['sample_index'].append(sample_index)
        data['accept_ratio'].append(file_data['accept_count']/steps)


    df = pd.DataFrame(data)

    print(df)

    exp_table = df.loc[(df.method == 'lm') | (df.method == 'mm') ]
    print("Mean acceptance ratio:", np.mean(exp_table.accept_ratio))

    
    mean_data = {
        'method': [],
        'N': [],
        'T': [],
    }

    for mx in df.method.unique():
        mx_table = df.loc[(df.method == mx)]
        for nx in mx_table.N.unique():
            nx_mx_table = mx_table.loc[mx_table.N == nx]
            mean_time = np.mean(nx_mx_table['T'])

            mean_data['method'].append(mx)
            mean_data['N'].append(nx)
            mean_data['T'].append(mean_time)
    

    df_mean = pd.DataFrame(mean_data)
    df_mean = df_mean.sort_values(by=['method', 'N'])

    print(df_mean)


    labels = {
        'lm': 'Local',
        'mm': 'Multipole',
        'dlm': 'DL\_MONTE',
    }

    linestyles = {
        'lm': '-.',
        'mm': '--',
        'dlm': '-',
    }

    linewidths = {
        'lm': 2,
        'mm': 1,
        'dlm': 1,
    }

    linecolours = {
        'lm': 'k',
        'mm': 'r',
        'dlm': 'b',
    }


    fig = plt.figure(figsize=(4, 3))
    ax = fig.add_subplot(111)

    for dx in ('lm', 'mm', 'dlm'):
        tx = df_mean.loc[df_mean.method == dx]
        print(tx)
        ax.plot(
            np.array(tx['N']), 
            np.array(tx['T']), 
            label=labels[dx],
            linestyle=linestyles[dx],
            linewidth=linewidths[dx],
            color=linecolours[dx]
        )
    

    ax.set_xscale('log')
    ax.set_yscale('log')


    ax.legend(loc='upper left', bbox_to_anchor=(0.01, 1.01), fancybox=False)
    ax.set_xlabel(r'Number of charges $N$')
    ax.set_ylabel(r'Time taken per step (ms)')

    fig.savefig('dlm_timing_plot.pdf', bbox_inches='tight')   







