from plot_DSD_common import plot_DSD
import matplotlib.pyplot as plt

data_labels = {}
data_colors = {}
data_nsd = {}

directory_base = "/home/piotr/praca/coal_fluctu_dim/LCM_DSD_fluctuations/data/LCM/Onishi/"

data_labels[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_ConstMulti1_Ens10_T300/"] = "DT 0.1 ConstMulti1"
data_nsd[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_ConstMulti1_Ens10_T300/"] = 64e6
#data_labels[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_ConstMulti1_Ens10_T300_RngSeedInit44/"] = "DT 0.1 ConstMulti1 RngSeedInit"

data_labels[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e5Tail_Ens100_T300/"] = "DT 0.1 SstpCoal1  Sd1e5Tail"
data_nsd[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e5Tail_Ens100_T300/"] = 1e5

data_labels[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e4Tail_Ens100_T300/"] = "DT 0.1 SstpCoal1  Sd1e4Tail"
data_nsd[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e4Tail_Ens100_T300/"] = 1e4

# nie ma roznic miedzy tymi trzema, dt i sstp_coal nie maja wplywu
data_labels[directory_base + "GA17_Np64e6_nx1_dt0.01_SstpCoal1_Onishi_SdConc1e3Tail_Ens100_T300/"] = "DT 0.01 SstpCoal1  Sd1e3Tail"
data_nsd[directory_base + "GA17_Np64e6_nx1_dt0.01_SstpCoal1_Onishi_SdConc1e3Tail_Ens100_T300/"] = 1e3
#data_labels[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e3Tail_Ens100_T300/"]  = "DT 0.1  SstpCoal1  Sd1e3Tail"
#data_labels[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal10_Onishi_SdConc1e3Tail_Ens100_T300/"] = "DT 0.1  SstpCoal10 Sd1e3Tail"

data_labels[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e2Tail_Ens1e4_T300/"] =               "DT 0.1  SstpCoal1  Sd1e2Tail"
data_nsd[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e2Tail_Ens1e4_T300/"] = 1e2
#data_labels[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e2Tail_Ens100_T300_RngSeedInit44/"] = "DT 0.1  SstpCoal10 Sd1e2Tail RngSeedInit"

data_labels[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e1Tail_Ens1e4_T300/"] =               "DT 0.1  SstpCoal1  Sd1e1Tail"
data_nsd[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e1Tail_Ens1e4_T300/"] = 1e1
#data_labels[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e1Tail_Ens100_T300_RngSeedInit44/"] = "DT 0.1  SstpCoal10 Sd1e1Tail RngSeedInit"

fig, ax = plt.subplots(1, 2, figsize=(24,9))

plot_DSD(data_labels, fig, ax, 0)
