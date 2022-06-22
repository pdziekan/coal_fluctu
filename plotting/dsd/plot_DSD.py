from plot_DSD_common import plot_DSD, plot_DSD_diff
import matplotlib.pyplot as plt

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


# nie ma roznic miedzy tymi trzema, dt i sstp_coal nie maja wplywu
data_labels = {}
#data_labels[directory_base + "GA17_Np64e6_nx1_dt0.01_SstpCoal1_Onishi_SdConc1e3Tail_Ens100_T300/"] = "DT 0.01 SstpCoal1  Sd1e3Tail"
data_labels[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e3Tail_Ens100_T300/"]  = "dt 0.1  SD 1e3 tail"
data_labels[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal10_Onishi_SdConc1e3Tail_Ens100_T300/"] = "dt 0.01 SD 1e3 tail"

fig, ax = plt.subplots(1, 2, figsize=(24,9))
plot_DSD(data_labels, data_color, fig, ax, 300, "timestep")


# porownanie symulacji z i bez losowosci w rozkladzie poczatkowym - nie ma znaczenia, losowosc po 300 s i tak jest taka sama
data_labels = {}
data_labels[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_ConstMulti1_Ens10_T300/"] = "one-to-one"
data_labels[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_ConstMulti1_Ens10_T300_RngSeedInit44/"] = "one-to-one, no randomness in initial DSD"
#data_labels[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e1Tail_Ens1e4_T300/"] =               "DT 0.1  SstpCoal1  Sd1e1Tail"
#data_labels[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e1Tail_Ens100_T300_RngSeedInit44/"] = "DT 0.1  SstpCoal10 Sd1e1Tail RngSeedInit"
data_labels[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e2Tail_Ens1e4_T300/"] =               "SD 1e2 tail"
data_labels[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e2Tail_Ens1e4_T300_RngSeedInit44/"] = "SD 1e2 tail, no randomness in initial DSD"

fig, ax = plt.subplots(1, 2, figsize=(24,9))
plot_DSD(data_labels, data_color, fig, ax, 300, "initial_randomness")
fig, ax = plt.subplots(1, 2, figsize=(24,9))
plot_DSD(data_labels, data_color, fig, ax, 0, "initial_randomness")


# porownanie symulacji dla roznej liczby SD
data_labels = {}
data_labels[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_ConstMulti1_Ens10_T300/"] =     "one-to-one"
data_labels[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e5Tail_Ens100_T300/"] =  "SD 1e5 tail"
data_labels[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e4Tail_Ens100_T300/"] =  "SD 1e4 tail"
data_labels[directory_base + "GA17_Np64e6_nx1_dt0.01_SstpCoal1_Onishi_SdConc1e3Tail_Ens100_T300/"] = "SD 1e3 tail"
data_labels[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e2Tail_Ens1e4_T300/"] =  "SD 1e2 tail"
data_labels[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e1Tail_Ens1e4_T300/"] =  "SD 1e1 tail"

fig, ax = plt.subplots(1, 2, figsize=(24,9))
plot_DSD(data_labels, data_color, fig, ax, 300, "NSD")
fig, ax = plt.subplots(1, 2, figsize=(24,9))
plot_DSD(data_labels, data_color, fig, ax, 0, "NSD")

plt.clf()
plot_DSD_diff(data_labels, data_color, "one-to-one", 300, "NSD")



