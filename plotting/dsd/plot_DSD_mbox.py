from plot_DSD_common import plot_DSD, plot_DSD_diff
import matplotlib.pyplot as plt

labeldict = {
 0 : "(a)",
 1 : "(b)",
 2 : "(c)",
 3 : "(d)",
 4 : "(e)",
 5 : "(f)",
 6 : "(g)",
 7 : "(h)",
 8 : "(i)",
 9 : "(j)",
 10 : "(k)",
 11 : "(l)",
}

directory_box = "/home/piotr/praca/coal_fluctu_dim/LCM_DSD_fluctuations/data/LCM/Onishi/"
directory_mbox = "/home/piotr/praca/coal_fluctu_dim/LCM_DSD_fluctuations/data/LCM_multi_box/Onishi/"

available_datasets = [
directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_ConstMulti1_Ens10_T300/",
directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e4Tail_Ens1e3_T300/",
directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e3Tail_Ens1e4_T300/",
directory_mbox + "03_07_2023_GA17_Np64e3_nx1e1_eps0.1_dt0.1var_sd1e1_Onishi_HallDavis/",
directory_mbox + "03_07_2023_GA17_Np64e6_nx1_eps0.1_dt0.1var_sd1e2_Onishi_HallDavis/",
directory_mbox + "04_07_2023_GA17_Np64e3_nx1e1_eps0.1_dt0.1var_sd1_Onishi_HallDavis/",
directory_mbox + "GA17_Np1_nx400_eps0.1_dt0.01var_ConstMulti1_Onishi_HallDavis/",
directory_mbox + "06_07_2023_GA17_Np64e3_nx1e1_eps0.1_dt0.1var_SdConc1Tail_Spinup30_Onishi_HallDavis/",
directory_mbox + "07_07_2023_GA17_Np64e3_nx1e1_eps0.1_dt0.1var_SdConc1Tail_DomainInit_Spinup30_Onishi_HallDavis/",
directory_mbox + "07_07_2023_GA17_NpTot64e6_nx2_eps0.1_dt0.1var_SdTot1e4Tail_DomainInit_Spinup30_Onishi_HallDavis/",
directory_mbox + "07_07_2023_GA17_NpTot64e6_nx5_eps0.1_dt0.1var_SdTot1e4Tail_DomainInit_Spinup30_Onishi_HallDavis/",
directory_mbox + "07_07_2023_GA17_NpTot64e6_nx10_eps0.1_dt0.1var_SdTot1e4Tail_DomainInit_Spinup30_Onishi_HallDavis/",
directory_mbox + "07_07_2023_GA17_NpTot64e6_nx20_eps0.1_dt0.1var_SdTot1e4Tail_DomainInit_Spinup30_Onishi_HallDavis/",
directory_mbox + "07_07_2023_GA17_NpTot64e6_nx2_eps0.1_dt0.1var_SdTot1e3Tail_DomainInit_Spinup30_Onishi_HallDavis/",
directory_mbox + "07_07_2023_GA17_NpTot64e6_nx5_eps0.1_dt0.1var_SdTot1e3Tail_DomainInit_Spinup30_Onishi_HallDavis/",
]

data_nx = {}
data_nx[directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_ConstMulti1_Ens10_T300/"] = 1
data_nx[directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e4Tail_Ens1e3_T300/"] = 1
data_nx[directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e3Tail_Ens1e4_T300/"] = 1
data_nx[directory_mbox + "03_07_2023_GA17_Np64e3_nx1e1_eps0.1_dt0.1var_sd1e1_Onishi_HallDavis/"] = 10
data_nx[directory_mbox + "03_07_2023_GA17_Np64e6_nx1_eps0.1_dt0.1var_sd1e2_Onishi_HallDavis/"] = 1
data_nx[directory_mbox + "04_07_2023_GA17_Np64e3_nx1e1_eps0.1_dt0.1var_sd1_Onishi_HallDavis/"] = 10
data_nx[directory_mbox + "GA17_Np1_nx400_eps0.1_dt0.01var_ConstMulti1_Onishi_HallDavis/"] = 400
data_nx[directory_mbox + "06_07_2023_GA17_Np64e3_nx1e1_eps0.1_dt0.1var_SdConc1Tail_Spinup30_Onishi_HallDavis/"] = 10
data_nx[directory_mbox + "07_07_2023_GA17_Np64e3_nx1e1_eps0.1_dt0.1var_SdConc1Tail_DomainInit_Spinup30_Onishi_HallDavis/"] = 10
data_nx[directory_mbox + "07_07_2023_GA17_NpTot64e6_nx2_eps0.1_dt0.1var_SdTot1e4Tail_DomainInit_Spinup30_Onishi_HallDavis/"] = 2
data_nx[directory_mbox + "07_07_2023_GA17_NpTot64e6_nx5_eps0.1_dt0.1var_SdTot1e4Tail_DomainInit_Spinup30_Onishi_HallDavis/"] = 5
data_nx[directory_mbox + "07_07_2023_GA17_NpTot64e6_nx10_eps0.1_dt0.1var_SdTot1e4Tail_DomainInit_Spinup30_Onishi_HallDavis/"] = 10
data_nx[directory_mbox + "07_07_2023_GA17_NpTot64e6_nx20_eps0.1_dt0.1var_SdTot1e4Tail_DomainInit_Spinup30_Onishi_HallDavis/"] = 20
data_nx[directory_mbox + "07_07_2023_GA17_NpTot64e6_nx2_eps0.1_dt0.1var_SdTot1e3Tail_DomainInit_Spinup30_Onishi_HallDavis/"] = 2
data_nx[directory_mbox + "07_07_2023_GA17_NpTot64e6_nx5_eps0.1_dt0.1var_SdTot1e3Tail_DomainInit_Spinup30_Onishi_HallDavis/"] = 5

# number of SDs in the entire box
data_nsd = {}
data_nsd[directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_ConstMulti1_Ens10_T300/"] = 64e6
data_nsd[directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e4Tail_Ens1e3_T300/"] = 1e4
data_nsd[directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e3Tail_Ens1e4_T300/"] = 1e3
data_nsd[directory_mbox + "03_07_2023_GA17_Np64e3_nx1e1_eps0.1_dt0.1var_sd1e1_Onishi_HallDavis/"] = 1e4
data_nsd[directory_mbox + "03_07_2023_GA17_Np64e6_nx1_eps0.1_dt0.1var_sd1e2_Onishi_HallDavis/"] = 1e2
data_nsd[directory_mbox + "04_07_2023_GA17_Np64e3_nx1e1_eps0.1_dt0.1var_sd1_Onishi_HallDavis/"] = 1e3
data_nsd[directory_mbox + "GA17_Np1_nx400_eps0.1_dt0.01var_ConstMulti1_Onishi_HallDavis/"] = 64e6
data_nsd[directory_mbox + "06_07_2023_GA17_Np64e3_nx1e1_eps0.1_dt0.1var_SdConc1Tail_Spinup30_Onishi_HallDavis/"] = 1e3
data_nsd[directory_mbox + "07_07_2023_GA17_Np64e3_nx1e1_eps0.1_dt0.1var_SdConc1Tail_DomainInit_Spinup30_Onishi_HallDavis/"] = 1e3
data_nsd[directory_mbox + "07_07_2023_GA17_NpTot64e6_nx2_eps0.1_dt0.1var_SdTot1e4Tail_DomainInit_Spinup30_Onishi_HallDavis/"] = 1e4
data_nsd[directory_mbox + "07_07_2023_GA17_NpTot64e6_nx5_eps0.1_dt0.1var_SdTot1e4Tail_DomainInit_Spinup30_Onishi_HallDavis/"] = 1e4
data_nsd[directory_mbox + "07_07_2023_GA17_NpTot64e6_nx10_eps0.1_dt0.1var_SdTot1e4Tail_DomainInit_Spinup30_Onishi_HallDavis/"] = 1e4
data_nsd[directory_mbox + "07_07_2023_GA17_NpTot64e6_nx20_eps0.1_dt0.1var_SdTot1e4Tail_DomainInit_Spinup30_Onishi_HallDavis/"] = 1e4
data_nsd[directory_mbox + "07_07_2023_GA17_NpTot64e6_nx2_eps0.1_dt0.1var_SdTot1e3Tail_DomainInit_Spinup30_Onishi_HallDavis/"] = 1e3
data_nsd[directory_mbox + "07_07_2023_GA17_NpTot64e6_nx5_eps0.1_dt0.1var_SdTot1e3Tail_DomainInit_Spinup30_Onishi_HallDavis/"] = 1e3


## linestyles denote dt
##linestyles = ['-', '--']
##ls_dt = {
##  1e10 : 'solid',
##  0.01 : 'dashed', 
##  0.1  : 'dotted', 
##  1    : 'dashdot', 
##  10   : (0, (3, 5, 1, 5, 1, 5))
##}
#ls_dt = {
#  1e10 : 'dotted',
#  0.01 : 'solid', 
#  0.1  : 'solid', 
#  1    : 'solid', 
#  10   : 'solid'
#}

# alternatively - line alpha denotes nx

# colors denote NSD
# colors fom the default
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']
lc_nsd = {
  1e1 : colors[0],
  1e2 : colors[1],
  1e3 : colors[2],
  1e4 : colors[3],
  1e5 : colors[4],
  64e6 : 'black'
}

data_color = {}
data_ls = {}
data_la = {}
data_nsd_conc = {} # initial mean number of SDs per cell
for ds in available_datasets:
  data_nsd_conc[ds] = float(data_nsd[ds]) / float(data_nx[ds])**3
  data_color[ds] = lc_nsd[data_nsd[ds]] 
  data_ls[ds] = 'solid'#ls_dt[data_dt[ds]] 
  data_la[ds] = 0.4 #+ (1-0.4) * () 




## ---- porownanie symulacji multi-box dla roznej liczby SD dla roznej wielkosci komorek ----
fig, ax = plt.subplots(3, 2, figsize=(8,12))

# one-to-one
data_labels = {}
data_labels[directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_ConstMulti1_Ens10_T300/"] =  "one-to-one"
data_labels[directory_mbox + "GA17_Np1_nx400_eps0.1_dt0.01var_ConstMulti1_Onishi_HallDavis/"] = "one-to-one nx400"
#plot_DSD(data_labels, data_color, data_ls, data_la, True, fig, ax[0], 0)
plot_DSD(data_labels, data_color, data_ls, data_la, False, fig, ax[0], 300)

# SD=1e4
data_labels = {}
data_labels[directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e4Tail_Ens1e3_T300/"] =  "$N_\mathrm{SD}^\mathrm{(bin)}=10^4$"
data_labels[directory_mbox + "07_07_2023_GA17_NpTot64e6_nx2_eps0.1_dt0.1var_SdTot1e4Tail_DomainInit_Spinup30_Onishi_HallDavis/"] = "SdTot 1e4 nx 2"
data_labels[directory_mbox + "07_07_2023_GA17_NpTot64e6_nx5_eps0.1_dt0.1var_SdTot1e4Tail_DomainInit_Spinup30_Onishi_HallDavis/"] = "SdTot 1e4 nx 5"
data_labels[directory_mbox + "07_07_2023_GA17_NpTot64e6_nx10_eps0.1_dt0.1var_SdTot1e4Tail_DomainInit_Spinup30_Onishi_HallDavis/"] = "SdTot 1e4 nx 10"
data_labels[directory_mbox + "07_07_2023_GA17_NpTot64e6_nx20_eps0.1_dt0.1var_SdTot1e4Tail_DomainInit_Spinup30_Onishi_HallDavis/"] = "SdTot 1e4 nx 20"
plot_DSD(data_labels, data_color, data_ls, data_la, False, fig, ax[1], 300)

# SD=1e3
data_labels = {}
data_labels[directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e3Tail_Ens1e4_T300/"] =  "$N_\mathrm{SD}^\mathrm{(bin)}=10^3$"
data_labels[directory_mbox + "07_07_2023_GA17_NpTot64e6_nx2_eps0.1_dt0.1var_SdTot1e3Tail_DomainInit_Spinup30_Onishi_HallDavis/"] = "SdTot 1e3 nx 2"
data_labels[directory_mbox + "07_07_2023_GA17_NpTot64e6_nx5_eps0.1_dt0.1var_SdTot1e3Tail_DomainInit_Spinup30_Onishi_HallDavis/"] = "SdTot 1e3 nx 5"
data_labels[directory_mbox + "07_07_2023_GA17_Np64e3_nx1e1_eps0.1_dt0.1var_SdConc1Tail_DomainInit_Spinup30_Onishi_HallDavis/"] = "SdTot 1e3 nx 10"
plot_DSD(data_labels, data_color, data_ls, data_la, False, fig, ax[2], 300)

#ax[0,0].set_xlabel('')
#ax[0,1].set_xlabel('')
#ax[0,0].set_title('$t=0\ s$')
#ax[0,1].set_title('$t=0\ s$')
#ax[1,0].set_title('$t=300\ s$')
#ax[1,1].set_title('$t=300\ s$')
#single legend for the whole figure
#handles, labels = ax[0,0].get_legend_handles_labels()
#lgd = fig.legend(handles, labels, handlelength=4, loc='lower center', bbox_to_anchor=(0.475,0))
# a b c d labels
for i,axs in enumerate(ax.flatten()):
  axs.text(0.1, 0.9, labeldict[i], fontsize=10, transform=axs.transAxes)
fig.tight_layout()
fig.subplots_adjust(bottom=0.28, wspace=0.3)#, hspace=0.25)
fig.savefig("/home/piotr/praca/coal_fluctu_dim/LCM_DSD_fluctuations/img/DSD_multibox_t0_t300.pdf")
