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

available_datasets = [
# const multi
directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_ConstMulti1_Ens10_T300/",
directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_ConstMulti1_Ens10_T300_RngSeedInit44/",

directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e5Tail_Ens100_T300/",

#directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e4Tail_Ens100_T300/",
#directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e4Tail_Ens1e2_T300/",
directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e4Tail_Ens1e3_T300/",

#directory_box + "GA17_Np64e6_nx1_dt0.01_SstpCoal1_Onishi_SdConc1e3Tail_Ens100_T300/",
directory_box + "GA17_Np64e6_nx1_dt0.01_SstpCoal1_Onishi_SdConc1e3Tail_Ens1e4_T300/",
#directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e3Tail_Ens100_T300/",
#directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal10_Onishi_SdConc1e3Tail_Ens100_T300/",
#directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e3Tail_Ens1e3_T300/",
#directory_box + "GA17_Np64e6_nx1_dt1_SstpCoal1_Onishi_SdConc1e3Tail_Ens1e3_T300/",
directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e3Tail_Ens1e4_T300/",
directory_box + "GA17_Np64e6_nx1_dt1_SstpCoal1_Onishi_SdConc1e3Tail_Ens1e4_T300/",
#directory_box + "GA17_Np64e6_nx1_dt10_SstpCoal1_Onishi_SdConc1e3Tail_Ens1e3_T300/",
directory_box + "GA17_Np64e6_nx1_dt10_SstpCoal1_Onishi_SdConc1e3Tail_Ens1e4_T300/",

directory_box + "GA17_Np64e6_nx1_dt0.01_SstpCoal1_Onishi_SdConc1e2Tail_Ens1e4_T300/",
directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e2Tail_Ens1e4_T300/",
directory_box + "GA17_Np64e6_nx1_dt1_SstpCoal1_Onishi_SdConc1e2Tail_Ens1e4_T300/",
directory_box + "GA17_Np64e6_nx1_dt10_SstpCoal1_Onishi_SdConc1e2Tail_Ens1e4_T300/",
directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e2Tail_Ens1e4_T300_RngSeedInit44/",

directory_box + "GA17_Np64e6_nx1_dt0.01_SstpCoal1_Onishi_SdConc1e1Tail_Ens1e5_T300/",
directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e1Tail_Ens1e5_T300/",
#directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e1Tail_Ens1e4_T300/",
directory_box + "GA17_Np64e6_nx1_dt1_SstpCoal1_Onishi_SdConc1e1Tail_Ens1e5_T300/",
directory_box + "GA17_Np64e6_nx1_dt10_SstpCoal1_Onishi_SdConc1e1Tail_Ens1e5_T300/",
]

data_nsd = {}

data_nsd[directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_ConstMulti1_Ens10_T300/"] = 64e6
data_nsd[directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_ConstMulti1_Ens10_T300_RngSeedInit44/"] = 64e6
data_nsd[directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e5Tail_Ens100_T300/"] = 1e5
data_nsd[directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e4Tail_Ens1e3_T300/"] = 1e4
data_nsd[directory_box + "GA17_Np64e6_nx1_dt0.01_SstpCoal1_Onishi_SdConc1e3Tail_Ens1e4_T300/"] = 1e3
data_nsd[directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e3Tail_Ens1e4_T300/"] = 1e3
data_nsd[directory_box + "GA17_Np64e6_nx1_dt1_SstpCoal1_Onishi_SdConc1e3Tail_Ens1e4_T300/"] = 1e3
data_nsd[directory_box + "GA17_Np64e6_nx1_dt10_SstpCoal1_Onishi_SdConc1e3Tail_Ens1e4_T300/"] = 1e3
data_nsd[directory_box + "GA17_Np64e6_nx1_dt0.01_SstpCoal1_Onishi_SdConc1e2Tail_Ens1e4_T300/"] = 1e2
data_nsd[directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e2Tail_Ens1e4_T300/"] = 1e2
data_nsd[directory_box + "GA17_Np64e6_nx1_dt1_SstpCoal1_Onishi_SdConc1e2Tail_Ens1e4_T300/"] = 1e2
data_nsd[directory_box + "GA17_Np64e6_nx1_dt10_SstpCoal1_Onishi_SdConc1e2Tail_Ens1e4_T300/"] = 1e2
data_nsd[directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e2Tail_Ens1e4_T300_RngSeedInit44/"] = 1e2
data_nsd[directory_box + "GA17_Np64e6_nx1_dt0.01_SstpCoal1_Onishi_SdConc1e1Tail_Ens1e5_T300/"] = 1e1
data_nsd[directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e1Tail_Ens1e5_T300/"] = 1e1
data_nsd[directory_box + "GA17_Np64e6_nx1_dt1_SstpCoal1_Onishi_SdConc1e1Tail_Ens1e5_T300/"] = 1e1
data_nsd[directory_box + "GA17_Np64e6_nx1_dt10_SstpCoal1_Onishi_SdConc1e1Tail_Ens1e5_T300/"] = 1e1

data_dt = {}
data_dt[directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_ConstMulti1_Ens10_T300/"] = 1e10 # 1e10 used to mark variable dt, lame...
data_dt[directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_ConstMulti1_Ens10_T300_RngSeedInit44/"] = 1e10
data_dt[directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e5Tail_Ens100_T300/"] = 0.1
data_dt[directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e4Tail_Ens1e3_T300/"] = 0.1
data_dt[directory_box + "GA17_Np64e6_nx1_dt0.01_SstpCoal1_Onishi_SdConc1e3Tail_Ens1e4_T300/"] = 0.01
data_dt[directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e3Tail_Ens1e4_T300/"] = 0.1
data_dt[directory_box + "GA17_Np64e6_nx1_dt1_SstpCoal1_Onishi_SdConc1e3Tail_Ens1e4_T300/"] = 1
data_dt[directory_box + "GA17_Np64e6_nx1_dt10_SstpCoal1_Onishi_SdConc1e3Tail_Ens1e4_T300/"] = 10
data_dt[directory_box + "GA17_Np64e6_nx1_dt0.01_SstpCoal1_Onishi_SdConc1e2Tail_Ens1e4_T300/"] = 0.01
data_dt[directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e2Tail_Ens1e4_T300/"] = 0.1
data_dt[directory_box + "GA17_Np64e6_nx1_dt1_SstpCoal1_Onishi_SdConc1e2Tail_Ens1e4_T300/"] = 1
data_dt[directory_box + "GA17_Np64e6_nx1_dt10_SstpCoal1_Onishi_SdConc1e2Tail_Ens1e4_T300/"] = 10
data_dt[directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e2Tail_Ens1e4_T300_RngSeedInit44/"] = 0.1
data_dt[directory_box + "GA17_Np64e6_nx1_dt0.01_SstpCoal1_Onishi_SdConc1e1Tail_Ens1e5_T300/"] = 0.01
data_dt[directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e1Tail_Ens1e5_T300/"] = 0.1
data_dt[directory_box + "GA17_Np64e6_nx1_dt1_SstpCoal1_Onishi_SdConc1e1Tail_Ens1e5_T300/"] = 1
data_dt[directory_box + "GA17_Np64e6_nx1_dt10_SstpCoal1_Onishi_SdConc1e1Tail_Ens1e5_T300/"] = 10


# linestyles denote dt
#linestyles = ['-', '--']
#ls_dt = {
#  1e10 : 'solid',
#  0.01 : 'dashed', 
#  0.1  : 'dotted', 
#  1    : 'dashdot', 
#  10   : (0, (3, 5, 1, 5, 1, 5))
#}
ls_dt = {
  1e10 : 'dotted',
  0.01 : 'solid', 
  0.1  : 'solid', 
  1    : 'solid', 
  10   : 'solid'
}

# alternatively - line alpha denotes dt
la_dt = {
  1e10 : 1, 
  0.01 : 0.85,
  0.1  : 0.7,
  1    : 0.55,
  10   : 0.4
}

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
for ds in available_datasets:
  data_color[ds] = lc_nsd[data_nsd[ds]] 
  data_ls[ds] = ls_dt[data_dt[ds]] 
  data_la[ds] = la_dt[data_dt[ds]] 

#override line style for const rng seed init runs
data_ls[directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_ConstMulti1_Ens10_T300_RngSeedInit44/"] =    'dashdot' 
data_ls[directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e2Tail_Ens1e4_T300_RngSeedInit44/"] = 'dashdot' 


# wplyw dt
data_labels = {}

#data_labels[directory_box + "GA17_Np64e6_nx1_dt10_SstpCoal1_Onishi_SdConc1e3Tail_Ens1e4_T300/"]  = "const SD, $N_\mathrm{SD}^\mathrm{(bin)}=1e3$, large tail init, $\Delta t_\mathrm{coal} = 10\ \mathrm{s}$"
#data_labels[directory_box + "GA17_Np64e6_nx1_dt1_SstpCoal1_Onishi_SdConc1e3Tail_Ens1e4_T300/"]  = "const SD, $N_\mathrm{SD}^\mathrm{(bin)}=1e3$, large tail init, $\Delta t_\mathrm{coal} = 1\ \mathrm{s}$"
##data_labels[directory_box + "GA17_Np64e6_nx1_dt0.01_SstpCoal1_Onishi_SdConc1e3Tail_Ens100_T300/"] = "DT 0.01 SstpCoal1  Sd1e3Tail"
#data_labels[directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e3Tail_Ens1e4_T300/"]  = "const SD, $N_\mathrm{SD}^\mathrm{(bin)}=1e3$, large tail init, $\Delta t_\mathrm{coal} = 0.1\ \mathrm{s}$"
#data_labels[directory_box + "GA17_Np64e6_nx1_dt0.01_SstpCoal1_Onishi_SdConc1e3Tail_Ens1e4_T300/"] = "const SD, $N_\mathrm{SD}^\mathrm{(bin)}=1e3$, large tail init, $\Delta t_\mathrm{coal} = 0.01\ \mathrm{s}$"

data_labels[directory_box + "GA17_Np64e6_nx1_dt0.01_SstpCoal1_Onishi_SdConc1e2Tail_Ens1e4_T300/"] = "$\Delta t_\mathrm{coal} = 0.01\ \mathrm{s}$"
data_labels[directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e2Tail_Ens1e4_T300/"] =  "$\Delta t_\mathrm{coal} = 0.1\ \mathrm{s}$"
data_labels[directory_box + "GA17_Np64e6_nx1_dt1_SstpCoal1_Onishi_SdConc1e2Tail_Ens1e4_T300/"] =    "$\Delta t_\mathrm{coal} = 1\ \mathrm{s}$"
data_labels[directory_box + "GA17_Np64e6_nx1_dt10_SstpCoal1_Onishi_SdConc1e2Tail_Ens1e4_T300/"] =   "$\Delta t_\mathrm{coal} = 10\ \mathrm{s}$"

fig, ax = plt.subplots(1, 2, figsize=(8,4))
plot_DSD(data_labels, data_color, data_ls, data_la, False, fig, ax, 300)
#single legend for the whole figure
handles, labels = ax[0].get_legend_handles_labels()
lgd = fig.legend(handles, labels, handlelength=4, loc='lower center', bbox_to_anchor=(0.475,0))
# a b c d labels
for i,axs in enumerate(ax.flatten()):
  axs.text(0.1, 0.9, labeldict[i], fontsize=10, transform=axs.transAxes)
fig.tight_layout()
fig.subplots_adjust(bottom=0.4, wspace=0.3)#, hspace=0.25)
#fig.suptitle("$N_\mathrm{SD}^\mathrm{(bin)}=10^2$")
fig.savefig("/home/piotr/praca/coal_fluctu_dim/LCM_DSD_fluctuations/img/DSD_timestep_t300.pdf")



# porownanie symulacji z i bez losowosci w rozkladzie poczatkowym - nie ma znaczenia, losowosc po 300 s i tak jest taka sama
data_labels = {}
data_labels[directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_ConstMulti1_Ens10_T300/"] = "one-to-one"
data_labels[directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_ConstMulti1_Ens10_T300_RngSeedInit44/"] = "one-to-one, no randomness in initial DSD"
#data_labels[directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e1Tail_Ens1e4_T300/"] =               "DT 0.1  SstpCoal1  Sd1e1Tail"
#data_labels[directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e1Tail_Ens100_T300_RngSeedInit44/"] = "DT 0.1  SstpCoal10 Sd1e1Tail RngSeedInit"
data_labels[directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e2Tail_Ens1e4_T300/"] =               "$N_\mathrm{SD}^\mathrm{(bin)}=10^2$, $\Delta t_\mathrm{coal} = 0.1\ \mathrm{s}$"
data_labels[directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e2Tail_Ens1e4_T300_RngSeedInit44/"] = "$N_\mathrm{SD}^\mathrm{(bin)}=10^2$, $\Delta t_\mathrm{coal} = 0.1\ \mathrm{s}$, no randomness in initial DSD"

fig, ax = plt.subplots(2, 2, figsize=(8,8))
plot_DSD(data_labels, data_color, data_ls, data_la, False, fig, ax[0], 0)
plot_DSD(data_labels, data_color, data_ls, data_la, False, fig, ax[1], 300)
ax[0,0].set_xlabel('')
ax[0,1].set_xlabel('')
ax[0,0].set_title('$t=0\ s$')
ax[0,1].set_title('$t=0\ s$')
ax[1,0].set_title('$t=300\ s$')
ax[1,1].set_title('$t=300\ s$')
#single legend for the whole figure
handles, labels = ax[0,0].get_legend_handles_labels()
lgd = fig.legend(handles, labels, handlelength=4, loc='lower center', bbox_to_anchor=(0.475,0))
# a b c d labels
for i,axs in enumerate(ax.flatten()):
  axs.text(0.1, 0.9, labeldict[i], fontsize=10, transform=axs.transAxes)
fig.tight_layout()
fig.subplots_adjust(bottom=0.2, wspace=0.3)#, hspace=0.25)
fig.savefig("/home/piotr/praca/coal_fluctu_dim/LCM_DSD_fluctuations/img/DSD_initial_randomness_t0_t300.pdf")


## porownanie symulacji dla roznej liczby SD
data_labels = {}
data_labels[directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e1Tail_Ens1e5_T300/"] =  "$N_\mathrm{SD}=10^1$"
data_labels[directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e2Tail_Ens1e4_T300/"] =  "$N_\mathrm{SD}=10^2$"
data_labels[directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e3Tail_Ens1e4_T300/"] =  "$N_\mathrm{SD}=10^3$"
data_labels[directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e4Tail_Ens1e3_T300/"] =  "$N_\mathrm{SD}=10^4$"
data_labels[directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e5Tail_Ens100_T300/"] =  "$N_\mathrm{SD}=10^5$"
data_labels[directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_ConstMulti1_Ens10_T300/"] =     "reference"

fig, ax = plt.subplots(1, 2, figsize=(8,5))

plot_DSD(data_labels, data_color, data_ls, data_la, False, fig, ax, 300)
ax[0].set_ylabel('')
ax[1].set_ylabel('')

ax[0].set_title('mean droplet mass distribution')
ax[1].set_title('std. dev. of droplet mass distribution')
#single legend for the whole figure
handles, labels = ax[0].get_legend_handles_labels()
lgd = fig.legend(handles, labels, handlelength=4, loc='lower center', bbox_to_anchor=(0.475,0))
fig.tight_layout()
fig.subplots_adjust(bottom=0.35, wspace=0.3)#, hspace=0.25)
fig.savefig("/home/piotr/praca/coal_fluctu_dim/LCM_DSD_fluctuations/img/DSD_NSD_t300.png", dpi=300)


## porownanie symulacji dla roznej liczby SD
data_labels = {}
data_labels[directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e1Tail_Ens1e5_T300/"] =  "$N_\mathrm{SD}^\mathrm{(bin)}=10^1$"
data_labels[directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e2Tail_Ens1e4_T300/"] =  "$N_\mathrm{SD}^\mathrm{(bin)}=10^2$"
data_labels[directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e3Tail_Ens1e4_T300/"] =  "$N_\mathrm{SD}^\mathrm{(bin)}=10^3$"
data_labels[directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e4Tail_Ens1e3_T300/"] =  "$N_\mathrm{SD}^\mathrm{(bin)}=10^4$"
data_labels[directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e5Tail_Ens100_T300/"] =  "$N_\mathrm{SD}^\mathrm{(bin)}=10^5$"
data_labels[directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_ConstMulti1_Ens10_T300/"] =     "one-to-one"

fig, ax = plt.subplots(2, 2, figsize=(8,8))

plot_DSD(data_labels, data_color, data_ls, data_la, True, fig, ax[0], 0)
plot_DSD(data_labels, data_color, data_ls, data_la, True, fig, ax[1], 300)
ax[0,0].set_xlabel('')
ax[0,1].set_xlabel('')
ax[0,0].set_title('$t=0\ s$')
ax[0,1].set_title('$t=0\ s$')
ax[1,0].set_title('$t=300\ s$')
ax[1,1].set_title('$t=300\ s$')
#single legend for the whole figure
handles, labels = ax[0,0].get_legend_handles_labels()
lgd = fig.legend(handles, labels, handlelength=4, loc='lower center', bbox_to_anchor=(0.475,0))
# a b c d labels
for i,axs in enumerate(ax.flatten()):
  axs.text(0.1, 0.9, labeldict[i], fontsize=10, transform=axs.transAxes)
fig.tight_layout()
fig.subplots_adjust(bottom=0.28, wspace=0.3)#, hspace=0.25)
fig.savefig("/home/piotr/praca/coal_fluctu_dim/LCM_DSD_fluctuations/img/DSD_NSD_t0_t300.pdf")




# diff plots in single file
fig, ax = plt.subplots(2, 2, figsize=(9,9))

# roznica SD wzgledem one-to-one, stale dt=0.1
data_labels = {}
data_labels[directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_ConstMulti1_Ens10_T300/"] =     "one-to-one" # 3.2e8 pairs / s
data_labels[directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e1Tail_Ens1e5_T300/"] =  "$N_\mathrm{SD}^\mathrm{(bin)}=10^1$"
data_labels[directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e2Tail_Ens1e4_T300/"] =  "$N_\mathrm{SD}^\mathrm{(bin)}=10^2$"
data_labels[directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e3Tail_Ens1e4_T300/"]  = "$N_\mathrm{SD}^\mathrm{(bin)}=10^3$"
data_labels[directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e4Tail_Ens1e3_T300/"] =  "$N_\mathrm{SD}^\mathrm{(bin)}=10^4$"
data_labels[directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e5Tail_Ens100_T300/"] =  "$N_\mathrm{SD}^\mathrm{(bin)}=10^5$"
#data_labels[directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal10_Onishi_SdConc1e3Tail_Ens100_T300/"] = "5e4 pairs / s, SD 1e3 tail"
#data_labels[directory_box + "GA17_Np64e6_nx1_dt1_SstpCoal1_Onishi_SdConc1e3Tail_Ens1e3_T300/"]  = "5e2 pairs / s, SD 1e3 tail"
#data_labels[directory_box + "GA17_Np64e6_nx1_dt10_SstpCoal1_Onishi_SdConc1e3Tail_Ens1e3_T300/"]  = "5e1 pairs / s, SD 1e3 tail"
#data_labels[directory_box + "GA17_Np64e6_nx1_dt1_SstpCoal1_Onishi_SdConc1e2Tail_Ens1e4_T300/"] =    "5e1 pairs / s, SD 1e2 tail"
#data_labels[directory_box + "GA17_Np64e6_nx1_dt1_SstpCoal1_Onishi_SdConc1e1Tail_Ens1e5_T300/"] =  "5e0 pairs / s, SD 1e1 tail"
#data_labels[directory_box + "GA17_Np64e6_nx1_dt10_SstpCoal1_Onishi_SdConc1e2Tail_Ens1e4_T300/"] =    "5e0 pairs / s, SD 1e2 tail"
plot_DSD_diff(ax[1,1], data_labels, data_color, data_ls, data_la, "one-to-one", 300, "NSD_dt0.1")

# roznica dt wzgledem one-to-one, stale SD=1e1
data_labels = {}
data_labels[directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_ConstMulti1_Ens10_T300/"] =     "one-to-one" # 3.2e8 pairs / s
data_labels[directory_box + "GA17_Np64e6_nx1_dt0.01_SstpCoal1_Onishi_SdConc1e1Tail_Ens1e5_T300/"] =  "$\Delta t_\mathrm{coal} = 0.01\ \mathrm{s}$"
data_labels[directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e1Tail_Ens1e5_T300/"] =   "$\Delta t_\mathrm{coal} = 0.1\ \mathrm{s}$"
data_labels[directory_box + "GA17_Np64e6_nx1_dt1_SstpCoal1_Onishi_SdConc1e1Tail_Ens1e5_T300/"] =     "$\Delta t_\mathrm{coal} = 1\ \mathrm{s}$"
data_labels[directory_box + "GA17_Np64e6_nx1_dt10_SstpCoal1_Onishi_SdConc1e1Tail_Ens1e5_T300/"] =    "$\Delta t_\mathrm{coal} = 10\ \mathrm{s}$"
plot_DSD_diff(ax[0,0], data_labels, data_color, data_ls, data_la, "one-to-one", 300, "dt_NSD1e1")

# roznica dt wzgledem one-to-one, stale SD=1e2
data_labels = {}
data_labels[directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_ConstMulti1_Ens10_T300/"] =     "one-to-one" # 3.2e8 pairs / s
data_labels[directory_box + "GA17_Np64e6_nx1_dt0.01_SstpCoal1_Onishi_SdConc1e2Tail_Ens1e4_T300/"] =  "$\Delta t_\mathrm{coal} = 0.01\ \mathrm{s}$"
data_labels[directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e2Tail_Ens1e4_T300/"] =   "$\Delta t_\mathrm{coal} = 0.1\ \mathrm{s}$"
data_labels[directory_box + "GA17_Np64e6_nx1_dt1_SstpCoal1_Onishi_SdConc1e2Tail_Ens1e4_T300/"] =     "$\Delta t_\mathrm{coal} = 1\ \mathrm{s}$"
data_labels[directory_box + "GA17_Np64e6_nx1_dt10_SstpCoal1_Onishi_SdConc1e2Tail_Ens1e4_T300/"] =    "$\Delta t_\mathrm{coal} = 10\ \mathrm{s}$"
plot_DSD_diff(ax[0,1], data_labels, data_color, data_ls, data_la, "one-to-one", 300, "dt_NSD1e2")

# roznica dt wzgledem one-to-one, stale SD=1e3
data_labels = {}
data_labels[directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_ConstMulti1_Ens10_T300/"] =     "one-to-one" # 3.2e8 pairs / s
#data_labels[directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal10_Onishi_SdConc1e3Tail_Ens100_T300/"] = "5e4 pairs / s, SD 1e3 tail, $\Delta t_\mathrm{coal} = 0.01\ \mathrm{s}$"
data_labels[directory_box + "GA17_Np64e6_nx1_dt0.01_SstpCoal1_Onishi_SdConc1e3Tail_Ens1e4_T300/"]  = "$\Delta t_\mathrm{coal} = 0.01\ \mathrm{s}$"
data_labels[directory_box + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e3Tail_Ens1e4_T300/"]  =  "$\Delta t_\mathrm{coal} = 0.1\ \mathrm{s}$"
data_labels[directory_box + "GA17_Np64e6_nx1_dt1_SstpCoal1_Onishi_SdConc1e3Tail_Ens1e4_T300/"]  =    "$\Delta t_\mathrm{coal} = 1\ \mathrm{s}$"
data_labels[directory_box + "GA17_Np64e6_nx1_dt10_SstpCoal1_Onishi_SdConc1e3Tail_Ens1e4_T300/"]  =   "$\Delta t_\mathrm{coal} = 10\ \mathrm{s}$"
plot_DSD_diff(ax[1,0], data_labels, data_color, data_ls, data_la, "one-to-one", 300, "dt_NSD1e3")

ax[1,1].set_title("$\Delta t_\mathrm{coal} = 0.1\ \mathrm{s}$")
ax[0,0].set_title("$N_\mathrm{SD}^\mathrm{(bin)}=10^1$")
ax[0,1].set_title("$N_\mathrm{SD}^\mathrm{(bin)}=10^2$")
ax[1,0].set_title("$N_\mathrm{SD}^\mathrm{(bin)}=10^3$")

#ax[0,0].set_ylim(-0.06, 0.05)
#ax[0,1].set_ylim(-0.06, 0.05)
#ax[1,0].set_ylim(-0.06, 0.05)
#ax[1,1].set_ylim(-0.06, 0.05)

ax[0,0].set_xlabel('')
ax[0,1].set_xlabel('')
ax[0,1].set_ylabel('')
ax[1,1].set_ylabel('')
# a b c d labels
for i,axs in enumerate(ax.flatten()):
  axs.text(0.9, 0.9, labeldict[i], fontsize=10, transform=axs.transAxes)
fig.tight_layout()
plt.savefig("/home/piotr/praca/coal_fluctu_dim/LCM_DSD_fluctuations/img/DSD_diff.pdf")


## ----- multi-box plots, has to be at the end because linestyles are overridden ----------

directory_mbox = "/home/piotr/praca/coal_fluctu_dim/LCM_DSD_fluctuations/data/LCM_multi_box/Onishi/"

available_datasets = [
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

#1e10 used to mark variable dt, lame...
data_dt[directory_mbox + "03_07_2023_GA17_Np64e3_nx1e1_eps0.1_dt0.1var_sd1e1_Onishi_HallDavis/"] = 1e10
data_dt[directory_mbox + "03_07_2023_GA17_Np64e6_nx1_eps0.1_dt0.1var_sd1e2_Onishi_HallDavis/"] = 1e10
data_dt[directory_mbox + "04_07_2023_GA17_Np64e3_nx1e1_eps0.1_dt0.1var_sd1_Onishi_HallDavis/"] = 1e10
data_dt[directory_mbox + "GA17_Np1_nx400_eps0.1_dt0.01var_ConstMulti1_Onishi_HallDavis/"] = 1e10
data_dt[directory_mbox + "06_07_2023_GA17_Np64e3_nx1e1_eps0.1_dt0.1var_SdConc1Tail_Spinup30_Onishi_HallDavis/"] = 1e10
data_dt[directory_mbox + "07_07_2023_GA17_Np64e3_nx1e1_eps0.1_dt0.1var_SdConc1Tail_DomainInit_Spinup30_Onishi_HallDavis/"] = 1e10
data_dt[directory_mbox + "07_07_2023_GA17_NpTot64e6_nx2_eps0.1_dt0.1var_SdTot1e4Tail_DomainInit_Spinup30_Onishi_HallDavis/"] = 1e10
data_dt[directory_mbox + "07_07_2023_GA17_NpTot64e6_nx5_eps0.1_dt0.1var_SdTot1e4Tail_DomainInit_Spinup30_Onishi_HallDavis/"] = 1e10
data_dt[directory_mbox + "07_07_2023_GA17_NpTot64e6_nx10_eps0.1_dt0.1var_SdTot1e4Tail_DomainInit_Spinup30_Onishi_HallDavis/"] = 1e10
data_dt[directory_mbox + "07_07_2023_GA17_NpTot64e6_nx20_eps0.1_dt0.1var_SdTot1e4Tail_DomainInit_Spinup30_Onishi_HallDavis/"] = 1e10
data_dt[directory_mbox + "07_07_2023_GA17_NpTot64e6_nx2_eps0.1_dt0.1var_SdTot1e3Tail_DomainInit_Spinup30_Onishi_HallDavis/"] = 1e10
data_dt[directory_mbox + "07_07_2023_GA17_NpTot64e6_nx5_eps0.1_dt0.1var_SdTot1e3Tail_DomainInit_Spinup30_Onishi_HallDavis/"] = 1e10

#override line style and line alpha for multi-box simulations
for ds in available_datasets:
  #if ds in data_nx and data_nx[ds] > 1:
  if ds in data_nx:
    if data_nx[ds] == 1:
      data_ls[ds] = 'solid'
      data_la[ds] = 0.4
    else:
      data_ls[ds] = 'dashdot'
      nsd_conc = float(data_nsd[ds]) / float(data_nx[ds])**3
      data_la[ds] = 0.4 + (1 - 0.4) * (1e4 - nsd_conc) / 1e4



## porownanie symulacji multi-box dla roznej liczby SD dla roznej wielkosci komorek
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
