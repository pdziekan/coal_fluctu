import numpy as np
import matplotlib.pyplot as plt

four_over_three_pi_rhow = 4./3. * np.pi * 1e3 # [kg/m3]
volume = 0.451 # cell volume [m3]

def plot_DSD(data_labels, data_colors, fig, ax, time, outname):
  for pre in data_labels:
    LCM_r         = np.zeros(0)
    LCM_m_r    = np.zeros(0)
    LCM_m_r_std_dev  = np.zeros(0)
  
    fs = open(pre+"size_spectr_mean.dat","r")
    rows = [x.split() for x in fs.readlines()]
    for row in rows[1:]: # first row is data description
      LCM_r                = np.append(LCM_r,                float(row[0])) # [um]
      if time == 0:
        LCM_m_r          = np.append(LCM_m_r,          float(row[1])) # same
        LCM_m_r_std_dev  = np.append(LCM_m_r_std_dev,  float(row[2])) # same
      else:
        LCM_m_r          = np.append(LCM_m_r,          float(row[3])) # same
        LCM_m_r_std_dev  = np.append(LCM_m_r_std_dev,  float(row[4])) # same
  
    LCM_dlogr = np.log(LCM_r[1]) - np.log(LCM_r[0]) # [log(um)]
    LCM_r *= 1e-6 # [m]
  
    LCM_m_logr = LCM_m_r * LCM_r # [g/m3 / unit(log[um])] 
    LCM_m_logr_std_dev = LCM_m_r_std_dev * LCM_r # [g/m3 / unit(log[um])] 
    LCM_M = LCM_m_logr * LCM_dlogr * volume * 1e-3 # mass of droplets of this size in the cell [kg]
    LCM_N = LCM_M / four_over_three_pi_rhow / np.power(LCM_r,3) # number of droplets of this size in the cell [1]
    print("total number of droplets in the cell in LCM at the end: ", np.sum(LCM_N))
    print("total mass of droplets in the cell in LCM at the end: ", np.sum(LCM_M))
  
  # plot m(log r)
    ax[0].plot(LCM_r * 1e6, LCM_m_logr, color=data_colors[pre], label= data_labels[pre])
    ax[1].plot(LCM_r * 1e6, LCM_m_logr_std_dev, color=data_colors[pre], label= data_labels[pre])
  # std dev scaled as sqrt(N_SD)
  #  ax[1].plot(LCM_r * 1e6, LCM_m_logr_std_dev  * np.sqrt(data_nsd[pre]) / np.sqrt(64e6), label='end ' + data_labels[pre])
  
  # plot m(r)
  #  ax[0].plot(LCM_r * 1e6, LCM_m_r, label='end ' + data_labels[pre])
  #  ax[1].plot(LCM_r * 1e6, LCM_m_r_std_dev, label='end ' + data_labels[pre])
  
  # plot N (number of droplets in the bin)
  #  ax[0].plot(LCM_r * 1e6, LCM_N, label='end ' + data_labels[pre])
    
  
  #  ax[0].plot(LCM_r, LCM_m_r, label='init ' + data_labels[pre]) # LCM data is mass density m(r) with LCM_rius in meters
  #  ax[1].plot(LCM_r, LCM_m_r_std_dev, label='init ' + data_labels[pre])
   # ax[1].plot(LCM_r, LCM_m_r_std_dev * np.sqrt(data_nsd[pre]) / np.sqrt(64e6), label='init ' + data_labels[pre])
   # ax[1].plot(LCM_r, LCM_m_r_std_dev  * np.sqrt(data_nsd[pre]) / np.sqrt(64e6), label='end ' + data_labels[pre])
  
  # EFM results (Smoluchowski)
  efm_data = np.genfromtxt("/home/piotr/praca/coal_fluctu_dim/LCM_DSD_fluctuations/data/EFM/dt 0.1_rq015.0_xmw 2.0_scal10.0_isw3_cut****_rmr****_tmax301._boplot00.out", unpack=True)
  #efm_data = np.genfromtxt("/home/piotr/praca/coal_fluctu_dim/LCM_DSD_fluctuations/data/EFM/dt 0.1_rq015.0_xmw 2.0_scal15.0_isw3_cut****_rmr****_tmax301._boplot00.out", unpack=True)
  # efm_data[0] - radius[um]
  # efm_data[1] - mass density m(ln r) [kg/m3] with radius in microns
  EFM_r = efm_data[0] * 1e-6 # [m]
  EFM_m_logr = efm_data[1] # [kg/m3 / unit(log[um])]
  EFM_dlogr = np.log(efm_data[0][1]) - np.log(efm_data[0][0]) # [log(um)]
  EFM_m_r = EFM_m_logr / EFM_r # mass density m(r) [kg/m3 / m]
  EFM_M = EFM_m_logr * EFM_dlogr * volume # mass of droplets of this size in the cell [kg]
  EFM_N = EFM_M / four_over_three_pi_rhow / np.power(EFM_r,3) # number of droplets of this size in the cell [1]
  EFM_M_std_dev = four_over_three_pi_rhow * np.power(EFM_r,3) * np.sqrt(EFM_N) # std_dev of total mass of droplets of this size in the whole cell [kg]
  # both give the same
  EFM_m_logr_std_dev = EFM_M_std_dev / EFM_dlogr / volume # std_dev of mass density m(log r) [kg/m3]
  #EFM_m_logr_std_dev = np.sqrt(four_over_three_pi_rhow * np.power(EFM_r,3) * EFM_m_logr / EFM_dlogr / volume) # std_dev of mass density m(log r) [kg/m3]
  
  # plot m(log r)
  if time == 0:
    EFM_m_logr = EFM_m_logr[:int(EFM_N.size/2)]
    EFM_m_logr_std_dev = EFM_m_logr_std_dev[:int(EFM_N.size/2)]
  else:
    EFM_m_logr = EFM_m_logr[int(EFM_N.size/2):]
    EFM_m_logr_std_dev = EFM_m_logr_std_dev[int(EFM_N.size/2):]

  ax[0].plot(EFM_r[:int(EFM_N.size/2)] * 1e6, EFM_m_logr * 1e3, ls='--', color='black', label="SCE") # radius in [um], mass density in [g/m3 / unit(log[um])]
  if time > 0: # dont plot SCE initial std dev - it is zero...
    ax[1].plot(EFM_r[:int(EFM_N.size/2)] * 1e6, EFM_m_logr_std_dev * 1e3, ls='--', color='black', label="SCE estimate") # radius in [um], mass density in [g/m3 / unit(log[um])]
  
  # plot N (number of droplets in the bin)
  #ax[0].plot(EFM_r * 1e6, EFM_N, label='end ' + data_labels[pre])
  
  #  ax[0].plot(EFM_r * 1e6, efm_data[1] / efm_data[0] * 1e6 * 1e3, label="EFM") 
  #  ax[1].plot(EFM_r * 1e6, np.sqrt(four_over_three_pi_rhow * np.power(efm_data[0] * 1e-6,3) * EFM_M), label="EFM") # std_dev of total mass of droplets of this size in the whole cell [kg]
  #  ax[0].plot(efm_data[0], efm_data[1] * EFM_dlogr / four_over_three_pi_rhow / np.power(efm_data[0] * 1e-6,3), label="EFM") # EFM data is mass density m(ln r) [kg/m3] with radius in microns; LCM data is m(r) with radius in meters
  
  print("total number of droplets in the cell in EFM at the end: ", np.sum(EFM_N[int(EFM_N.size/2):])) # EFM_N contains init spectrum and final spectrum
  print("total mass of droplets in the cell in EFM at the end: ", np.sum(EFM_M[int(EFM_M.size/2):]))
  
  #plt.title('Time at which largest droplet exceeds 3 mm radius and size of the droplet')
  ax[0].set_xlabel('radius [um]')
  ax[1].set_xlabel('radius [um]')
  #ax[0].set_ylabel('mean mass density m(r) [g/m^3 / m]')
  ax[0].set_ylabel('mean mass density m(log r) [g/m^3 / unit(log[um])')
  ax[1].set_ylabel('standard deviation of mass density m(log r) [g/m^3 / unit(log[um])')
  
  ax[0].set_xscale('log')
  ax[0].set_yscale('log')
  ax[1].set_xscale('log')
  ax[1].set_yscale('log')
  ax[0].set_xlim(1,1e3)
  ax[1].set_xlim(1,1e3)
  ax[0].set_ylim(1e-4,1e1)
  ax[1].set_ylim(1e-6,1e1)
  
  #ax[0].set_xlim(40,)
  #ax[0].set_ylim(0,0.30)
  
  plt.grid(axis='y', color='0.5')
  ax[0].legend()
  ax[1].legend()
  fig.savefig("/home/piotr/praca/coal_fluctu_dim/LCM_DSD_fluctuations/img/DSD_"+str(outname)+"_t"+str(time)+".svg")
  #plt.show()
