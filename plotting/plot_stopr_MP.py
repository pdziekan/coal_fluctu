import numpy as np
import matplotlib.pyplot as plt

data_labels = {}
data_colors = {}

directory = "/home/piotr/praca/coal_fluctu_dim/Smolu_min_cell_Marshall_Palmer/data/LCM/zero_dim/GA17_Np64e6_nx1_dt3var_MP_R10_HallDavis_cutoff2.5mm_stopr3mm_0"
data_labels[directory+ "/"] = "LCM Np64e6 nx1 stopr3mm, var dt 3, GA17 libcloudphxx eps 0.1, Marshall-Palmer R10, HallDavis"
data_labels[directory+ "/"] = "LCM one large coalescence cell"
data_colors[directory+ "/"] = "red"



directory = "/home/piotr/praca/coal_fluctu_dim/Smolu_min_cell_Marshall_Palmer/data/LCM/synth_turb_from_libcloudphxx/GA17_Np1_nx400_eps0.1_dt0.01var_MP_R10_HallDavis_cutoff2.5mm_stopr3mm_part"
data_labels[directory+ "/"] = "LCM Np1e0 nx400 stopr3mm, var dt 0.1, GA17 libcloudphxx eps 0.1, Marshall-Palmer R10, HallDavis"
data_labels[directory+ "/"] = "LCM multiple coalescence cells, each with 1 droplet on average"
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
  print(all_stopr)
  for row in all_stopr:
    rad = np.append(rad,float(row[3]))
    time = np.append(time,float(row[6]))
  print(time,rad)
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
plt.title('Time at which largest droplet exceeds 3 mm radius and size of the droplet')
ax.set_xlabel('time [s]')
ax.set_ylabel('radius [um]')

plt.legend()
fig.savefig("/home/piotr/praca/coal_fluctu_dim/lucky_droplets_vs_cell_size/img/lucky_droplets_vs_cell_size_LCM_rain.png")
plt.show()
#
#
