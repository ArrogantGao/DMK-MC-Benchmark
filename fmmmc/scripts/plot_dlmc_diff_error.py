import json
import numpy as np
from glob import glob
import os
import yaml
import pandas as pd
import matplotlib.pyplot as plt
from scipy import constants as spc
from math import exp, log

plt.rc("font", family="serif", serif="Times New Roman")
plt.rc("xtick", labelsize="small")
plt.rc("ytick", labelsize="small")
plt.rc("text", usetex=True)
plt.rc("savefig", dpi=500)
plt.rc("xtick", labelsize="small")
plt.rc("ytick", labelsize="small")
plt.rcParams.update({"font.size": 11})
params = {"text.latex.preamble": [r"\usepackage{amsmath}"]}
plt.rcParams.update(params)


def get_L_set():
    with open("input.json") as fh:
        data = json.loads(fh.read())
    Lset = data["L_set"]
    return Lset


if __name__ == "__main__":

    Lset = get_L_set()
    output_files = glob("*/err_data.json")
    energy_files = glob("*/energy_data.json")

    df = pd.concat([pd.read_json(fx) for fx in output_files])
    df_energy = pd.concat([pd.read_json(fx) for fx in energy_files])
    df["rel_err"] = df.apply(lambda row: abs(row["diff"] - row["correct_diff"]) / abs(row["correct_diff"]), axis=1)
    df["abs_err"] = df.apply(lambda row: abs(row["diff"] - row["correct_diff"]), axis=1)
    df["abs_correct_diff"] = df.apply(lambda row: abs(row["correct_diff"]), axis=1)
    

    temperature = 273.0
    ikbT = (-1.0) / ((spc.Boltzmann / spc.elementary_charge) * temperature)
    def prob_rel_err(row):
        return abs(exp(ikbT * row['diff']) - exp(ikbT * row['correct_diff'])) / exp(ikbT * row['correct_diff'])
    df["rel_err_prob"] = df.apply(prob_rel_err, axis=1)


    print(df_energy)
    print(df)

    df_means = df.groupby(["method", "p"]).mean().reset_index()
    dlm_mean_error = float(df_means.loc[df_means.method == "dlm"].rel_err_prob)



    diff_err_lm = df.loc[(df.method == 'lm') & (df.p == 12)]
    diff_err_mm = df.loc[(df.method == 'lm') & (df.p == 12)]
    diff_err_dlm = df.loc[(df.method == 'dlm')]

    

    bins = np.logspace(-8, -2, 50)


    fig, ax = plt.subplots(3,1,figsize=(4,2.0), sharex=True, sharey=False)
    
    tx = 6E-9

    ax[0].hist(diff_err_lm.rel_err_prob,  bins, color='0.50', edgecolor='k')
    ax[0].text(tx, 1800, r"Local")
    ax[1].hist(diff_err_mm.rel_err_prob,  bins, color='0.50', edgecolor='k')
    ax[1].text(tx, 1800, r"Multipole")
    ax[2].hist(diff_err_dlm.rel_err_prob, bins, color='0.50', edgecolor='k')
    ax[2].text(tx, 1800, r"DL\_MONTE")
    for axx in (0, 1, 2):
        ax[axx].set_xscale("log")
        ax[axx].set_ylim(top=2700)
        ax[axx].yaxis.set_major_locator(plt.MaxNLocator(4))


    ax[1].set_ylabel('Number of samples')
    ax[2].set_xlabel('Relative error $\delta P$ in acceptance probability')
    fig.savefig("diff_error_hists.pdf", bbox_inches="tight")


    emin = -9
    emax = -1
    bins = np.logspace(emin, emax, 50)
    mbins = np.array([0.5*(bins[bx+1] + bins[bx]) for bx in range(bins.shape[0]-1)])
    bin_widths = np.array([(bins[bx+1] - bins[bx]) for bx in range(bins.shape[0]-1)])



    cdiff_lm  = np.cumsum(np.multiply(np.histogram(diff_err_lm.rel_err_prob,  bins=bins, density=True)[0], bin_widths))
    cdiff_mm  = np.cumsum(np.multiply(np.histogram(diff_err_mm.rel_err_prob,  bins=bins, density=True)[0], bin_widths))
    cdiff_dlm = np.cumsum(np.multiply(np.histogram(diff_err_dlm.rel_err_prob, bins=bins, density=True)[0], bin_widths))
    

    

    fig, ax = plt.subplots(3,1,figsize=(4,2.0), sharex=True, sharey=False)
    
    tx = 6E-9

    ax[0].bar(mbins, cdiff_lm, width=bin_widths, color='0.50', edgecolor='k')
    #ax[0].text(tx, 1800, r"Local")
    ax[1].bar(mbins, cdiff_mm, width=bin_widths, color='0.50', edgecolor='k')
    #ax[1].text(tx, 1800, r"Multipole")
    ax[2].bar(mbins, cdiff_dlm, width=bin_widths, color='0.50', edgecolor='k')
    #ax[2].text(tx, 1800, r"DL\_MONTE")
    for axx in (0, 1, 2):
        ax[axx].set_xscale("log")
        #ax[axx].yaxis.set_major_locator(plt.MaxNLocator(4))

    ax[2].set_xticks([10**bx for bx in range(emin, emax+1)])
    ax[2].set_xlabel('Relative error $\delta P$ in acceptance probability')
    ax[1].set_ylabel('Cumulative density')
    fig.savefig("cdiff_error_hists.pdf", bbox_inches="tight")

    import pdb; pdb.set_trace()




    fig, ax = plt.subplots(1,3,figsize=(7,2.2), sharey=True)

    size = 0.2
    alpha = 0.5

    ax[0].scatter(diff_err_lm.rel_err_prob, diff_err_lm.abs_correct_diff, size, alpha=alpha, color='k')
    ax[0].set_title('Local')

    ax[1].scatter(diff_err_mm.rel_err_prob, diff_err_mm.abs_correct_diff, size, alpha=alpha, color='k')
    ax[1].set_title('Multipole')

    ax[2].scatter(diff_err_dlm.rel_err_prob, diff_err_dlm.abs_correct_diff, size, alpha=alpha, color='k')
    ax[2].set_title(r'DL\_MONTE 2')

    #ax.scatter(diff_err_mm.abs_correct_diff, diff_err_mm.rel_err)
    #ax.scatter(diff_err_dlm.abs_correct_diff, diff_err_dlm.rel_err)

    ax[0].set_ylabel(r"Correct Energy Difference $\delta U^*$")
    ax[1].set_xlabel('Relative error $\delta P$ in acceptance probability')

    for axx in (0, 1, 2):
        ax[axx].set_yscale("log")
        ax[axx].set_xlim((1E-9, 1E-1))
        ax[axx].set_ylim((0.000001, 1.0))
        ax[axx].set_xscale("log")


    fig.savefig("diff_error_plot.pdf", bbox_inches="tight")
    fig.savefig("diff_error_plot.png", bbox_inches="tight", dpi=600)




    labels = {"lm": "Local", "mm": "Multipole"}
    linestyles = {
        "lm": "-",
        "mm": "--",
    }
    linewidths = {
        "lm": 2,
        "mm": 1,
    }
    linecolours = {
        "lm": "k",
        "mm": "r",
    }
    fig = plt.figure(figsize=(4, 3))
    ax = fig.add_subplot(111)
    for dx in ("lm", "mm"):
        y = np.array(df_means.loc[df_means.method == dx].rel_err_prob)
        print(Lset)
        print(y)
        ax.plot(Lset, y, label=labels[dx], linestyle=linestyles[dx], linewidth=linewidths[dx], color=linecolours[dx])
    ax.axhline(dlm_mean_error, color="k", linestyle=":", label=r"DL\_MONTE 2")
    ax.set_yscale("log")
    ax.set_xticks(Lset)
    ax.legend(loc="upper right", bbox_to_anchor=(1.01, 1.01), fancybox=False)
    ax.set_xlabel(r"Number of Expansion Terms")
    ax.set_ylabel(r"Mean relative Error $\delta P$\\of Acceptance Probability")
    fig.savefig("error_plot.pdf", bbox_inches="tight")


    

