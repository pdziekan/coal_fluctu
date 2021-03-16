import numpy as np
import matplotlib.pyplot as plt

plots = ["tau", "nrain", "rmax"]

data = {
        "/mnt/local/pdziekan/wyniki/coalfluctu_dim/ST_Np1e3_nx20_eps1_dt0.1var_Wang/build/" : "1"
        ,"/mnt/local/pdziekan/wyniki/coalfluctu_dim/ST_Np1e3_nx20_eps1_dt0.1var_Wang/build1/" : "2"
        ,"/mnt/local/pdziekan/wyniki/coalfluctu_dim/ST_Np1e3_nx20_eps1_dt0.1var_Wang/build2/" : "3"
        ,"/mnt/local/pdziekan/wyniki/coalfluctu_dim/ST_Np1e3_nx20_eps1_dt0.1var_Wang/build3/" : "4"
       }


for plot in plots:
  fig, ax = plt.subplots()
  for pre in data:
    time = np.loadtxt(pre+"time.dat")
    series = np.loadtxt(pre+plot+".dat")
    ax.plot(time, series, label=data[pre])
  
  plt.xlim([2000, 5000])
  plt.legend()
  
  fig.savefig("/home/pracownicy/pdziekan/tmp/"+plot+".png")


