import numpy as np
import matplotlib.pyplot as plt

data_labels = {}
data_colors = {}

directory = "/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/data/OnishihalfN/LCM/zero_dim/GA17_Np27e6_nx1_eps0.1_dt1var_OnishihalfN_HallDavis_stopr300"
data_labels[directory+ "/"] = "LCM Np27e6 nx1 stopr300, var dt 1, GA17 libcloudphxx eps 0.1, OnishihalfN, HallDavis"
data_labels[directory+ "/"] = "LCM one large coalescence cell"
data_colors[directory+ "/"] = "red"

directory = "/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/data/OnishihalfN/LCM/zero_dim/GA17_Np27e6_nx1_eps0.1_dt1var_OnishihalfN_HallDavis_stopr300_2"
data_labels[directory+ "/"] = "LCM Np27e6 nx1 stopr300, var dt 1, GA17 libcloudphxx eps 0.1, OnishihalfN, HallDavis RE2"
data_labels[directory+ "/"] = ""
data_colors[directory+ "/"] = "red"

directory = "/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/data/OnishihalfN/LCM/zero_dim/GA17_Np27e6_nx1_eps0.1_dt1var_OnishihalfN_HallDavis_stopr300_3"
data_labels[directory+ "/"] = "LCM Np27e6 nx1 stopr300, var dt 1, GA17 libcloudphxx eps 0.1, OnishihalfN, HallDavis RE3"
data_labels[directory+ "/"] = ""
data_colors[directory+ "/"] = "red"

directory = "/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/data/OnishihalfN/LCM/zero_dim/GA17_Np27e6_nx1_eps0.1_dt1var_OnishihalfN_HallDavis_stopr300_4"
data_labels[directory+ "/"] = "LCM Np27e6 nx1 stopr300, var dt 1, GA17 libcloudphxx eps 0.1, OnishihalfN, HallDavis RE4"
data_labels[directory+ "/"] = ""
data_colors[directory+ "/"] = "red"

directory = "/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/data/OnishihalfN/LCM/zero_dim/GA17_Np27e6_nx1_eps0.1_dt1var_OnishihalfN_HallDavis_stopr300_5"
data_labels[directory+ "/"] = "LCM Np27e6 nx1 stopr300, var dt 1, GA17 libcloudphxx eps 0.1, OnishihalfN, HallDavis RE5"
data_labels[directory+ "/"] = ""
data_colors[directory+ "/"] = "red"



directory = "/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/data/OnishihalfN/LCM/GA17_from_libcloudphxx/GA17_Np1e0_nx300_eps0.1_dt1var_OnishihalfN_HallDavis_stopr300"
data_labels[directory+ "/"] = "LCM Np1e0 nx300 stopr300, var dt 1, GA17 libcloudphxx eps 0.1, OnishihalfN, HallDavis"
data_labels[directory+ "/"] = "LCM multiple coalescence cells, each with 1 droplet on average"
data_colors[directory+ "/"] = "blue"

directory = "/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/data/OnishihalfN/LCM/GA17_from_libcloudphxx/GA17_Np1e0_nx300_eps0.1_dt1var_OnishihalfN_HallDavis_stopr300_2"
data_labels[directory+ "/"] = "LCM Np1e0 nx300 stopr300, var dt 1, GA17 libcloudphxx eps 0.1, OnishihalfN, HallDavis, RE2"
data_labels[directory+ "/"] = ""
data_colors[directory+ "/"] = "blue"

directory = "/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/data/OnishihalfN/LCM/GA17_from_libcloudphxx/GA17_Np1e0_nx300_eps0.1_dt1var_OnishihalfN_HallDavis_stopr300_3"
data_labels[directory+ "/"] = "LCM Np1e0 nx300 stopr300, var dt 1, GA17 libcloudphxx eps 0.1, OnishihalfN, HallDavis, RE3"
data_labels[directory+ "/"] = ""
data_colors[directory+ "/"] = "blue"

directory = "/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/data/OnishihalfN/LCM/GA17_from_libcloudphxx/GA17_Np1e0_nx300_eps0.1_dt1var_OnishihalfN_HallDavis_stopr300_4"
data_labels[directory+ "/"] = "LCM Np1e0 nx300 stopr300, var dt 1, GA17 libcloudphxx eps 0.1, OnishihalfN, HallDavis, RE4"
data_labels[directory+ "/"] = ""
data_colors[directory+ "/"] = "blue"

directory = "/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/data/OnishihalfN/LCM/GA17_from_libcloudphxx/GA17_Np1e0_nx300_eps0.1_dt1var_OnishihalfN_HallDavis_stopr300_5"
data_labels[directory+ "/"] = "LCM Np1e0 nx300 stopr300, var dt 1, GA17 libcloudphxx eps 0.1, OnishihalfN, HallDavis, RE5"
data_labels[directory+ "/"] = ""
data_colors[directory+ "/"] = "blue"




fig, ax = plt.subplots(figsize=(12,9))

aggregated_rad = {} # dictionary of numpy arrays containing all data of given type, where all data with the same color are assumed to be of the same type
aggregated_time = {}

for pre in data_colors:
  aggregated_rad[data_colors[pre]] = np.zeros(0)
  aggregated_time[data_colors[pre]] = np.zeros(0)

for pre in data_labels:
  rad = np.zeros(0)
  time = np.zeros(0)
  fs = open(pre+"rstop.dat","r")
  all_stopr=[x.split() for x in fs.readlines()]
#  print(all_stopr)
  for row in all_stopr:
    rad = np.append(rad,float(row[3]))
    time = np.append(time,float(row[5]))
#  print(time,rad)
  ax.scatter(time,rad, label=data_labels[pre], color=data_colors[pre], s=3, alpha=0.3)
  aggregated_time[data_colors[pre]] = np.append(aggregated_time[data_colors[pre]], time)
  aggregated_rad[data_colors[pre]] = np.append(aggregated_rad[data_colors[pre]], rad)

# plot mean from each type
for data_color in aggregated_time:
  time = aggregated_time[data_color]
  rad = aggregated_rad[data_color]
  ax.errorbar(np.mean(time),np.mean(rad), xerr = np.std(time) / np.sqrt(len(time)), yerr = np.std(rad) / np.sqrt(len(rad)), color=data_color)

#  for idx, (row_t, row_d) in enumerate(zip(all_time, all_data)):
#    time = np.array(row_t, dtype=np.float32)
#    series = np.array(row_d, dtype=np.float32)
#    ax.plot(time, series[0:len(time)], label=data[pre]+'('+str(idx)+')')
#
#
#
plt.title('Time at which largest droplet exceeds 300 um radius and size of the droplet')
ax.set_xlabel('time [s]')
ax.set_ylabel('radius [um]')

plt.legend()
fig.savefig("/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/img/stopr/Onishi_stopr_vs_cell_size_LCM.png")
#plt.show()
#
#
