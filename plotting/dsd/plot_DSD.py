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

directory_base = "/home/piotr/praca/coal_fluctu_dim/LCM_DSD_fluctuations/data/LCM/Onishi/"

#data_nsd = {}
#data_nsd[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_ConstMulti1_Ens10_T300/"] = 64e6
#data_nsd[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e5Tail_Ens100_T300/"] = 1e5
#data_nsd[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e4Tail_Ens100_T300/"] = 1e4
#data_nsd[directory_base + "GA17_Np64e6_nx1_dt0.01_SstpCoal1_Onishi_SdConc1e3Tail_Ens100_T300/"] = 1e3
#data_nsd[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e2Tail_Ens1e4_T300/"] = 1e2
#data_nsd[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e1Tail_Ens1e4_T300/"] = 1e1

# colors fom the default
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']
i = 0
data_color = {}
data_color[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e3Tail_Ens100_T300/"]               = colors[i] 
i = i+1
data_color[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal10_Onishi_SdConc1e3Tail_Ens100_T300/"]              = colors[i]
i = i+1
data_color[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_ConstMulti1_Ens10_T300/"]                  = colors[i]
i = i+1
data_color[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_ConstMulti1_Ens10_T300_RngSeedInit44/"]    = colors[i]
i = i+1
data_color[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e2Tail_Ens1e4_T300/"]               = colors[i]
i = i+1
data_color[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e2Tail_Ens1e4_T300_RngSeedInit44/"] = colors[i]
i = i+1
data_color[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e5Tail_Ens100_T300/"]               = colors[i]
i = i+1
data_color[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e4Tail_Ens100_T300/"]               = colors[i]
i = i+1
data_color[directory_base + "GA17_Np64e6_nx1_dt0.01_SstpCoal1_Onishi_SdConc1e3Tail_Ens100_T300/"]              = colors[i]
i = i+1
data_color[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e1Tail_Ens1e4_T300/"]               = colors[i]
i = i+1
data_color[directory_base + "GA17_Np64e6_nx1_dt1_SstpCoal1_Onishi_SdConc1e2Tail_Ens1e4_T300/"]               = colors[i]


# nie ma roznic miedzy tymi trzema, dt i sstp_coal nie maja wplywu
data_labels = {}
#data_labels[directory_base + "GA17_Np64e6_nx1_dt0.01_SstpCoal1_Onishi_SdConc1e3Tail_Ens100_T300/"] = "DT 0.01 SstpCoal1  Sd1e3Tail"
data_labels[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e3Tail_Ens100_T300/"]  = "const SD, $N_\mathrm{SD}=1e3$, large tail init, $\Delta t = 0.1\ \mathrm{s}$"
data_labels[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal10_Onishi_SdConc1e3Tail_Ens100_T300/"] = "const SD, $N_\mathrm{SD}=1e3$, large tail init, $\Delta t = 0.01\ \mathrm{s}$"
data_labels[directory_base + "GA17_Np64e6_nx1_dt1_SstpCoal1_Onishi_SdConc1e2Tail_Ens1e4_T300/"] = "const SD, $N_\mathrm{SD}=1e2$, large tail init, $\Delta t = 1\ \mathrm{s}$"
data_labels[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e2Tail_Ens1e4_T300/"] = "const SD, $N_\mathrm{SD}=1e2$, large tail init, $\Delta t = 0.1\ \mathrm{s}$"

fig, ax = plt.subplots(1, 2, figsize=(8,4))
plot_DSD(data_labels, data_color, False, fig, ax, 300)
#single legend for the whole figure
handles, labels = ax[0].get_legend_handles_labels()
lgd = fig.legend(handles, labels, handlelength=4, loc='lower center', bbox_to_anchor=(0.475,0))
# a b c d labels
for i,axs in enumerate(ax.flatten()):
  axs.text(0.1, 0.9, labeldict[i], fontsize=10, transform=axs.transAxes)
fig.tight_layout()
fig.subplots_adjust(bottom=0.26, wspace=0.3)#, hspace=0.25)
fig.savefig("/home/piotr/praca/coal_fluctu_dim/LCM_DSD_fluctuations/img/DSD_timestep_t300.pdf")



# porownanie symulacji z i bez losowosci w rozkladzie poczatkowym - nie ma znaczenia, losowosc po 300 s i tak jest taka sama
data_labels = {}
data_labels[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_ConstMulti1_Ens10_T300/"] = "one-to-one"
data_labels[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_ConstMulti1_Ens10_T300_RngSeedInit44/"] = "one-to-one, no randomness in initial DSD"
#data_labels[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e1Tail_Ens1e4_T300/"] =               "DT 0.1  SstpCoal1  Sd1e1Tail"
#data_labels[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e1Tail_Ens100_T300_RngSeedInit44/"] = "DT 0.1  SstpCoal10 Sd1e1Tail RngSeedInit"
data_labels[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e2Tail_Ens1e4_T300/"] =               "const SD, $N_\mathrm{SD}=1e2$, large tail init"
data_labels[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e2Tail_Ens1e4_T300_RngSeedInit44/"] = "const SD, $N_\mathrm{SD}=1e2$, large tail init, no randomness in initial DSD"

fig, ax = plt.subplots(2, 2, figsize=(8,8))
plot_DSD(data_labels, data_color, False, fig, ax[0], 0)
plot_DSD(data_labels, data_color, False, fig, ax[1], 300)
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
data_labels[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_ConstMulti1_Ens10_T300/"] =     "one-to-one"
data_labels[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e5Tail_Ens100_T300/"] =  "const SD, $N_\mathrm{SD}=1e5$, large tail init"
data_labels[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e4Tail_Ens100_T300/"] =  "const SD, $N_\mathrm{SD}=1e4$, large tail init"
data_labels[directory_base + "GA17_Np64e6_nx1_dt0.01_SstpCoal1_Onishi_SdConc1e3Tail_Ens100_T300/"] = "const SD, $N_\mathrm{SD}=1e3$, large tail init"
data_labels[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e2Tail_Ens1e4_T300/"] =  "const SD, $N_\mathrm{SD}=1e2$, large tail init"
data_labels[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e1Tail_Ens1e4_T300/"] =  "const SD, $N_\mathrm{SD}=1e1$, large tail init"


fig, ax = plt.subplots(2, 2, figsize=(8,8))
plot_DSD(data_labels, data_color, True, fig, ax[0], 0)
plot_DSD(data_labels, data_color, True, fig, ax[1], 300)
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

# roznica SD wzgledem one-to-one
data_labels = {}
data_labels[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_ConstMulti1_Ens10_T300/"] =     "one-to-one"
data_labels[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e5Tail_Ens100_T300/"] =  "SD 1e5 tail"
data_labels[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e2Tail_Ens1e4_T300/"] =  "SD 1e2 tail"
data_labels[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e1Tail_Ens1e4_T300/"] =  "SD 1e1 tail"
plot_DSD_diff(data_labels, data_color, "one-to-one", 300, "NSD")
