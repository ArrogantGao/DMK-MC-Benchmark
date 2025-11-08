'''
  =======================================================================
  Fit parameters for theoretical curves that model the time spent in the
  propose/accept stage for each of the two methods (multipole/local).
  =======================================================================

  Times are assumed to be of the form

   t(N) = tau[X][Y] + gamma[X][Y}*log_10(N)

  where X in {propose,accept} and Y in {local,multipole}.

  All data is taken from the mc_workspace repository, make sure you edit the
  paths accordingly.

  Results of the fit are used to calculate theoretical speedup curves for a
  given acceptance ratio nu. Generates a plot with those theoretical curves
  and the measured data.

'''

# Path to the git repository containing the data
mc_workspace_dir = '../../../mc_workspace/'

import json
import pandas
import os
from glob import *
import numpy as np
from matplotlib import pyplot as plt
import matplotlib

plt.rc('font', family='serif', serif='Times New Roman')
plt.rc('xtick', labelsize='small')
plt.rc('ytick', labelsize='small')
plt.rc('text', usetex=True)
plt.rc('savefig', dpi=500)
plt.rc('xtick', labelsize='small')
plt.rc('ytick', labelsize='small')
plt.rc('legend', edgecolor='lightgray')
plt.rcParams.update({'font.size': 18})
params = {'text.latex.preamble': [r'\usepackage{amsmath}']}
plt.rcParams.update(params)


def fit_data(filenames):
    '''
    Fit the time per propose/accept for the local/multipole method.

    Returns the fit coefficients tau and gamma as dictionaries
    
    :arg filenames: dictonary {'local':filename_local,
                               'multipole':filename_multipole}
         containing the names of files with the data
    '''
    # Fit parameters
    tau = {'accept':{'local':None,
                     'multipole':None},
           'propose':{'local':None,
                      'multipole':None}
    }

    gamma = {'accept':{'local':None,
                       'multipole':None},
             'propose':{'local':None,
                        'multipole':None}
    }

    for expansion in ('local','multipole'):
        data = np.load(filenames[expansion])
        ncharge     = np.array(data[:, 0], dtype=np.int)
        log_ncharge = np.log10(np.array(data[:, 0], dtype=np.int))
        t_prop      = 1.E6*data[:, 1].view()
        t_accept    = 1.E6*data[:, 2].view()
        plt.clf()
        ax = plt.gca()
        ax.set_xlim(1.E3,1.E6)
        ax.set_ylim(0,1800)
        ax.set_xscale('log')
        ax.set_xlabel(r'Number of charges $N$')
        ax.set_ylabel(r'time [$\mu\operatorname{s}]$')
        plt.plot(ncharge,t_prop,
                 linewidth=2,
                 marker='s',
                 markersize=6,
                 color='blue',
                 markerfacecolor='white',
                 markeredgecolor='blue',                 
                 markeredgewidth=2,
                 label='propose')
        plt.plot(ncharge,t_accept,
                 linewidth=2,
                 marker='o',
                 markersize=6,
                 color='red',
                 markerfacecolor='white',
                 markeredgecolor='red',
                 markeredgewidth=2,
                 label='accept')
        p_fit_prop=np.polyfit(log_ncharge, t_prop, 1)
        p_fit_accept=np.polyfit(log_ncharge, t_accept, 1)
        tau['propose'][expansion] = p_fit_prop[1]
        gamma['propose'][expansion] = p_fit_prop[0]
        tau['accept'][expansion] = p_fit_accept[1]
        gamma['accept'][expansion] = p_fit_accept[0]
        plt.plot(ncharge,p_fit_prop[1]+p_fit_prop[0]*log_ncharge,
                 linewidth=2,
                 linestyle='--',
                 color='blue')
        plt.plot(ncharge,p_fit_accept[1]+p_fit_accept[0]*log_ncharge,
                 linewidth=2,
                 linestyle='--',
                 color='red')
        plt.plot([1,1],[1,1],
                 linewidth=2,
                 linestyle='--',
                 color='black',
                 label='fit')
        plt.legend(loc='upper left')
        th = 100
        lsz=14
        if (expansion=='local'):
            ax.axvline(x=12525, color='slategray', linestyle='--')
            ax.axvline(x=100207, color='slategray', linestyle='--')
            ax.axvline(x=801663, color='slategray', linestyle='--')
            ax.text(3500,   th, r"$L=3$", ha="center", va="center", size=lsz)
            ax.text(35000,   th, r"$4$", ha="center", va="center", size=lsz)
            ax.text(310000,  th, r"$5$", ha="center", va="center", size=lsz)
            ax.text(900000, th, r"$6$", ha="center", va="center", size=lsz)
        if (expansion=='multipole'):
            ax.axvline(x=29681 , color='slategray', linestyle='--')
            ax.axvline(x=237449, color='slategray', linestyle='--')
            ax.text(7000,   th, r"$L=3$", ha="center", va="center", size=lsz)
            ax.text(120000,   th, r"$4$", ha="center", va="center", size=lsz)
            ax.text(590000,  th, r"$5$", ha="center", va="center", size=lsz)
        plt.savefig('times_'+expansion+'.pdf',bbox_inches='tight')
    tpl_dict = {}
    for stage in ('propose','accept'):
        for expansion in ('local','multipole'):
            tpl_dict['TAU_'+stage+'_'+expansion] = ('%5.1f' % tau[stage][expansion])
            tpl_dict['GAMMA_'+stage+'_'+expansion] = ('%5.1f' % gamma[stage][expansion])
    # Print LaTeX which can be copied-and-pasted
    with open('fit_coefficients.tpl') as f:
        for line in f.readlines():
            print (line.strip() % tpl_dict)
    return tau, gamma


def speedup_theoretical(tau,gamma,nu,N):
    '''
    Calculate theoretical speedup.

    This uses the fit parameters tau and gamma to calculate the speedup
    of the local method relative to the multipole method.

    :arg tau: fit coefficient tau
    :arg gamma: fit coefficient gamma
    :arg nu: acceptance rate
    :arg N: problem size
    '''
    def time(stage,expansion,N):
        '''
        Theoretical runtime for a particular stage and expansion method

        :arg stage: one of {'propose','accept'}
        :arg expansion: one of {'local','multipole'}
        :arg N: problem size (=number of charges)
        '''
        return tau[stage][expansion] + gamma[stage][expansion]*np.log10(N)

    t_multipole = 1./nu*time('propose','multipole',N) + time('accept','multipole',N)
    t_local = 1./nu*time('propose','local',N) + time('accept','local',N)
    return t_multipole/t_local

def speedup_theoretical_max(tau,gamma,nu):
    '''
    Calculate maximal achievable speedup.

    Computes the theoretical speedup in the limit of infinite system
    size.

    :arg tau: fit coefficient tau
    :arg gamma: fit coefficient gamma        
    :arg nu: acceptance ratio
    '''
    
    t_multipole = 1./nu*gamma['propose']['multipole'] + gamma['accept']['multipole']
    t_local = 1./nu*gamma['propose']['local'] + gamma['accept']['local']
    return t_multipole/t_local

def plot_speedup(tau,gamma,df):
    '''
    Plot both measured and theoretical speedup 

    Uses the measured data in the pandas data frame df and fit parameters
    to plot the measured and theoretical speedup as a function
    of the problem size

    :arg tau: fit coefficient tau
    :arg gamma: fit coefficient gamma
    :arg df: pandas dataframe with measured times and metadata
    '''
    
    nu = np.mean(df.accept_ratio)
    print("Mean acceptance ratio:", nu)

    print (gamma)
    N_min = 1.E3
    N_max = 1.E6
    plt.clf()
    ax = plt.gca()
    N = np.exp(np.arange(np.log(N_min),np.log(N_max),0.1))
    ax.set_xlim(N_min,N_max)
    ax.set_ylim(0,2.25)
    ax.set_xlabel('Number of charges $N$')
    ax.set_xscale('log')
    ax.set_ylabel('Speedup $S(N)$')
    plt.plot(N,speedup_theoretical(tau,gamma,nu,N),
             linewidth=2,
             linestyle='-',
             color='blue',
             label=r'theory [$\nu = '+('%4.3f' % nu)+r'$]')
    S_max = speedup_theoretical_max(tau,gamma,nu)
    plt.plot([N_min,N_max],[S_max,S_max],
             linewidth=2,
             linestyle='--',
             color='blue',
             label='asymptotic speedup')
    plt.annotate(r'$S_\infty = '+('%2.1f' % speedup_theoretical_max(tau,gamma,nu))+r'$         ',
                 (N_max*0.9,S_max-0.075),
                 horizontalalignment='right',
                 verticalalignment='top',
                 size=16,
                 color='blue')
    plt.plot([N_min,N_max],[1,1],
             linewidth=2,
             color='black',
             linestyle=':')
    # Measured speedups
    print (df)
    ncharge = np.array(df.N.loc[df.method == 'multipole'])
    t_meas = {}
    for expansion in df.method.unique():
        t_meas[expansion] = np.array(df.time.loc[df.method == expansion])
    plt.plot(ncharge,t_meas['multipole']/t_meas['local'],
             linewidth=2,
             color='red',
             marker='o',
             markersize=6,
             markeredgewidth=2,
             markerfacecolor='white',
             markeredgecolor='red',
             label='measured')
    plt.legend(loc='lower right')
    plt.savefig('speedup.pdf',bbox_inches='tight')

def read_measured_times(data_dir):
    '''
    Read the measured runtimes (and metadata) from a particular
    directory and return pandas dataframe

    :arg data_dir: directory to process
    '''
    with open(data_dir+'/input.json') as fh:
        input_data = json.loads(fh.read())
        steps = input_data['steps']
    print ('steps = ',steps)

    data = {'method': [],
            'N' : [],
            'time' : [],
            'sample_index': [],
            'accept_ratio': []}

    method_map = {'mm':'multipole','lm':'local'}
    
    lmmm = glob(data_dir+'/DLM*/timing_data*.json')
    for filex in lmmm:
        with open(filex) as fh:
            sample_index = int(os.path.basename(os.path.dirname(filex)).split('_')[2])
            file_data = json.loads(fh.read())
            data['method'].append(method_map[file_data['method']])
            data['N'].append(file_data['N'])
            data['sample_index'].append(sample_index)
            data['time'].append(1000 * file_data['time_taken'] / steps)
            data['accept_ratio'].append(file_data['accept_count']/steps)
    df = pandas.DataFrame(data)

    mean_data = {
        'method': [],
        'N': [],
        'time': [],
        'accept_ratio':[]
    }
    
    
    for mx in df.method.unique():
        mx_table = df.loc[(df.method == mx)]
        for nx in mx_table.N.unique():
            nx_mx_table = mx_table.loc[mx_table.N == nx]
            mean_time = np.mean(nx_mx_table['time'])
            mean_accept = np.mean(nx_mx_table['accept_ratio'])            
            mean_data['method'].append(mx)
            mean_data['N'].append(nx)
            mean_data['accept_ratio'].append(mean_accept)
            mean_data['time'].append(mean_time)

    df_mean = pandas.DataFrame(mean_data)
    df_mean = df_mean.sort_values(by=['method', 'N'])
    
    return df_mean

###################################################################
#              M A I N 
###################################################################
if (__name__=='__main__'):
    # Directory and files containing the data for propose/accpt only
    # timings
    fit_data_dir = '../order_N/ivb_vn_l12_2019_local/'
    filenames = {'local':fit_data_dir+'/ivb_balena_2019_near_tset_local.npy',
                 'multipole':fit_data_dir+'/ivb_balena_2019_near_tset_multipole.npy'}
    
    # Step 1: Fit coefficients using the propose/accept-only runs
    tau, gamma = fit_data(filenames)
    
    # Step 2: Extract data from full MC runs
    data_dir = '../dlmonte_compare/balena_ivb_4sample_2019_2/'
    df_mean = read_measured_times(data_dir)
    
    # Step 3: Plot theoretical and measured speedup
    plot_speedup(tau,gamma,df_mean)
