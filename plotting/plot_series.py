import numpy as np
import matplotlib.pyplot as plt

data = {
        "/mnt/local/pdziekan/wyniki/coalfluctu_dim/ST_Np1e3_nx20_eps1_dt0.1var_Wang/build/" : "1"
        ,"/mnt/local/pdziekan/wyniki/coalfluctu_dim/ST_Np1e3_nx20_eps1_dt0.1var_Wang/build1/" : "2"
        ,"/mnt/local/pdziekan/wyniki/coalfluctu_dim/ST_Np1e3_nx20_eps1_dt0.1var_Wang/build2/" : "3"
        ,"/mnt/local/pdziekan/wyniki/coalfluctu_dim/ST_Np1e3_nx20_eps1_dt0.1var_Wang/build3/" : "4"
       }

fig, ax = plt.subplots()

for pre in data:
  time = np.loadtxt(pre+"time.dat")
  tau = np.loadtxt(pre+"tau.dat")
  ax.plot(time, tau, label=data[pre])

plt.xlim([2000, 5000])
plt.legend()

fig.savefig("/home/pracownicy/pdziekan/tmp/tau.png")


