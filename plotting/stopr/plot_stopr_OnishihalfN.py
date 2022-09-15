import numpy as np
import matplotlib.pyplot as plt

data_labels = {}
data_colors = {}

directory = "/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/data/OnishihalfN/LCM/zero_dim/ares/GA17_Np27e6_nx1_eps0.1_dt1var_OnishihalfN_HallDavis_stopr300"
data_labels[directory + "/1/"]   = "nx1 ares float"
data_labels[directory + "/2/"]   = "nx1 ares GPU float, reszta double"
data_labels[directory + "/3/"]  =  "nx1 ares double"

data_colors[directory + "/1/"]   = "pink"
data_colors[directory + "/2/"]   = "pink"
data_colors[directory + "/3/"]   = "pink"

directory = "/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/data/OnishihalfN/LCM/GA17_from_libcloudphxx/ares/nx300/GA17_Np1e0_nx300_eps0.1_dt1var_OnishihalfN_HallDavis_stopr300"
data_labels[directory + "/1/"]   = "nx300 ares float"
data_labels[directory + "/2/"]   = ""
data_labels[directory + "/3/"]   = ""
data_labels[directory + "/4/"]   = ""
data_labels[directory + "/5/"]   = ""
data_labels[directory + "/6/"]   = ""
data_labels[directory + "/7/"]   = ""
data_labels[directory + "/8/"]   = ""
data_labels[directory + "/9/"]   = "nx300 ares GPU float, reszta double"
data_labels[directory + "/10/"]  = ""
data_labels[directory + "/11/"]  = ""
data_labels[directory + "/12/"]  = ""
data_labels[directory + "/13/"]  = ""
data_labels[directory + "/14/"]  = ""
data_labels[directory + "/15/"]  = ""
data_labels[directory + "/16/"]  = ""
data_labels[directory + "/17/"]  = "nx300 ares double"
data_labels[directory + "/18/"]  = ""
data_labels[directory + "/19/"]  = ""
data_labels[directory + "/20/"]  = ""
data_labels[directory + "/21/"]  = ""
data_labels[directory + "/22/"]  = ""
data_labels[directory + "/23/"]  = ""
data_labels[directory + "/24/"]  = ""

data_colors[directory + "/1/"]   = "orange"
data_colors[directory + "/2/"]   = "orange"
data_colors[directory + "/3/"]   = "orange"
data_colors[directory + "/4/"]   = "orange"
data_colors[directory + "/5/"]   = "orange"
data_colors[directory + "/6/"]   = "orange"
data_colors[directory + "/7/"]   = "orange"
data_colors[directory + "/8/"]   = "orange"
data_colors[directory + "/9/"]   = "yellow"
data_colors[directory + "/10/"]  = "yellow"
data_colors[directory + "/11/"]  = "yellow"
data_colors[directory + "/12/"]  = "yellow"
data_colors[directory + "/13/"]  = "yellow"
data_colors[directory + "/14/"]  = "yellow"
data_colors[directory + "/15/"]  = "yellow"
data_colors[directory + "/16/"]  = "yellow"
data_colors[directory + "/17/"]  = "green"
data_colors[directory + "/18/"]  = "green"
data_colors[directory + "/19/"]  = "green"
data_colors[directory + "/20/"]  = "green"
data_colors[directory + "/21/"]  = "green"
data_colors[directory + "/22/"]  = "green"
data_colors[directory + "/23/"]  = "green"
data_colors[directory + "/24/"]  = "green"


directory = "/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/data/OnishihalfN/LCM/zero_dim/cuda_k_4/GA17_Np27e6_nx1_eps0.1_dt1var_OnishihalfN_HallDavis_stopr300"
data_labels[directory+ "/"] = "LCM Np27e6 nx1 stopr300, var dt 1, GA17 libcloudphxx eps 0.1, OnishihalfN, HallDavis"
data_labels[directory+ "/"] = "cuda_k_4 nx1"
data_colors[directory+ "/"] = "red"

directory = "/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/data/OnishihalfN/LCM/zero_dim/cuda_k_4/GA17_Np27e6_nx1_eps0.1_dt1var_OnishihalfN_HallDavis_stopr300_2"
data_labels[directory+ "/"] = "LCM Np27e6 nx1 stopr300, var dt 1, GA17 libcloudphxx eps 0.1, OnishihalfN, HallDavis RE2"
data_labels[directory+ "/"] = ""
data_colors[directory+ "/"] = "red"

directory = "/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/data/OnishihalfN/LCM/zero_dim/cuda_k_4/GA17_Np27e6_nx1_eps0.1_dt1var_OnishihalfN_HallDavis_stopr300_3"
data_labels[directory+ "/"] = "LCM Np27e6 nx1 stopr300, var dt 1, GA17 libcloudphxx eps 0.1, OnishihalfN, HallDavis RE3"
data_labels[directory+ "/"] = ""
data_colors[directory+ "/"] = "red"

directory = "/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/data/OnishihalfN/LCM/zero_dim/cuda_k_4/GA17_Np27e6_nx1_eps0.1_dt1var_OnishihalfN_HallDavis_stopr300_4"
data_labels[directory+ "/"] = "LCM Np27e6 nx1 stopr300, var dt 1, GA17 libcloudphxx eps 0.1, OnishihalfN, HallDavis RE4"
data_labels[directory+ "/"] = ""
data_colors[directory+ "/"] = "red"

directory = "/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/data/OnishihalfN/LCM/zero_dim/cuda_k_4/GA17_Np27e6_nx1_eps0.1_dt1var_OnishihalfN_HallDavis_stopr300_5"
data_labels[directory+ "/"] = "LCM Np27e6 nx1 stopr300, var dt 1, GA17 libcloudphxx eps 0.1, OnishihalfN, HallDavis RE5"
data_labels[directory+ "/"] = ""
data_colors[directory+ "/"] = "red"
#
#
#
directory = "/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/data/OnishihalfN/LCM/GA17_from_libcloudphxx/cuda_k_4/nx300/GA17_Np1e0_nx300_eps0.1_dt1var_OnishihalfN_HallDavis_stopr300"
data_labels[directory+ "/"] = "LCM Np1e0 nx300 stopr300, var dt 1, GA17 libcloudphxx eps 0.1, OnishihalfN, HallDavis"
data_labels[directory+ "/"] = "cuda_k_4 nx300"
data_colors[directory+ "/"] = "blue"

directory = "/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/data/OnishihalfN/LCM/GA17_from_libcloudphxx/cuda_k_4/nx300/GA17_Np1e0_nx300_eps0.1_dt1var_OnishihalfN_HallDavis_stopr300_2"
data_labels[directory+ "/"] = "LCM Np1e0 nx300 stopr300, var dt 1, GA17 libcloudphxx eps 0.1, OnishihalfN, HallDavis, RE2"
data_labels[directory+ "/"] = ""
data_colors[directory+ "/"] = "blue"

directory = "/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/data/OnishihalfN/LCM/GA17_from_libcloudphxx/cuda_k_4/nx300/GA17_Np1e0_nx300_eps0.1_dt1var_OnishihalfN_HallDavis_stopr300_3"
data_labels[directory+ "/"] = "LCM Np1e0 nx300 stopr300, var dt 1, GA17 libcloudphxx eps 0.1, OnishihalfN, HallDavis, RE3"
data_labels[directory+ "/"] = ""
data_colors[directory+ "/"] = "blue"

directory = "/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/data/OnishihalfN/LCM/GA17_from_libcloudphxx/cuda_k_4/nx300/GA17_Np1e0_nx300_eps0.1_dt1var_OnishihalfN_HallDavis_stopr300_4"
data_labels[directory+ "/"] = "LCM Np1e0 nx300 stopr300, var dt 1, GA17 libcloudphxx eps 0.1, OnishihalfN, HallDavis, RE4"
data_labels[directory+ "/"] = ""
data_colors[directory+ "/"] = "blue"

directory = "/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/data/OnishihalfN/LCM/GA17_from_libcloudphxx/cuda_k_4/nx300/GA17_Np1e0_nx300_eps0.1_dt1var_OnishihalfN_HallDavis_stopr300_5"
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
    try:
      time = np.append(time,float(row[5]))
    except:
      time = np.append(time,float(row[6]))
#  print(time,rad)
  ax.scatter(time,rad, label=data_labels[pre], color=data_colors[pre], s=3, alpha=0.15)
  aggregated_time[data_colors[pre]] = np.append(aggregated_time[data_colors[pre]], time)
  aggregated_rad[data_colors[pre]] = np.append(aggregated_rad[data_colors[pre]], rad)

# plot mean from each type
for data_color in aggregated_time:
  time = aggregated_time[data_color]
  rad = aggregated_rad[data_color]
  # mean with mean square error
#  ax.errorbar(np.mean(time),np.mean(rad), xerr = np.std(time) / np.sqrt(len(time)), yerr = np.std(rad) / np.sqrt(len(rad)), color=data_color)
  # mean with 1 std dev 
  ax.errorbar(np.mean(time),np.mean(rad), xerr = np.std(time) , yerr = np.std(rad) , color=data_color)

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
fig.savefig("/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/img/stopr/OnishihalfN_stopr_vs_cell_size_LCM.pdf")
#plt.show()
#
#
