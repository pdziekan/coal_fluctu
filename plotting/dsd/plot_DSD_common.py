import numpy as np
import matplotlib.pyplot as plt

SE_scale = 1.96 # used to scale standard error, 1.96 gives 95% confidence interval
four_over_three_pi_rhow = 4./3. * np.pi * 1e3 # [kg/m3]
volume = 0.451 # cell volume [m3]

# terminal fall velocity at sea level according to Beard (1977)
def beard77(r): # r - radius [m]
  # use 3rd degree polynominal for r<20um 
  m_s = [0.105035e2, 0.108750e1, -0.133245, -0.659969e-2]
  # use 7th degree polynominal for r>20um
  m_l = [0.65639e1,    -0.10391e1,    -0.14001e1,    -0.82736e0,    -0.34277e0,    -0.83072e-1,    -0.10583e-1,    -0.54208e-3]

  x = np.log(2.*100.*(r))
  y = 0;
  # calc V0 (sea-level velocity)
  if(r <= 20e-6):
    for i in np.arange(4): 
      y += m_s[i] * pow(x, i);
  else:
    for i in np.arange(8): 
      y += m_l[i] * pow(x, i);

  return np.exp(y) / 100. # [m/s]

# calcualte mean size spectra from results of individual runs
# NOTE: assumnig that all used the same bins!
def mean_DSD(dicts):
  r_all = []
  m_pre_all = []
  m_post_all = []
  
  for pre in dicts:
    r = np.zeros(0)
    m_pre = np.zeros(0)
    m_post = np.zeros(0)

    fs = open(pre+"size_spectr.dat","r")
    rows = [x.split() for x in fs.readlines()]
    for row in rows: 
      if len(row) == 0: # empty row indicates new simulation
        r_all.append(r)
        m_pre_all.append(m_pre)
        m_post_all.append(m_post)

        r = np.zeros(0)
        m_pre = np.zeros(0)
        m_post = np.zeros(0)
      else:
        r      = np.append(r,      float(row[0])) # [um]
        m_pre  = np.append(m_pre,  float(row[1]))
        m_post = np.append(m_post, float(row[2]))
  
  return np.average(r_all, axis=0), np.average(m_pre_all, axis=0), np.std(m_pre_all, axis=0), np.average(m_post_all, axis=0), np.std(m_post_all, axis=0)
    

def read_DSD(pre, time):
  LCM_r      = np.zeros(0)
  LCM_m_r    = np.zeros(0)
  LCM_m_r_std_dev  = np.zeros(0)
  
  fs = open(pre+"size_spectr_mean.dat","r")
  rows = [x.split() for x in fs.readlines()]
  for row in rows[1:]: # first row is data description
    LCM_r              = np.append(LCM_r,            float(row[0])) # [um]
    if time == 0:
      LCM_m_r          = np.append(LCM_m_r,          float(row[1])) # same
      LCM_m_r_std_dev  = np.append(LCM_m_r_std_dev,  float(row[2])) # same
    else:
      LCM_m_r          = np.append(LCM_m_r,          float(row[3])) # same
      LCM_m_r_std_dev  = np.append(LCM_m_r_std_dev,  float(row[4])) # same
  
  LCM_r *= 1e-6 # [m]
  
  LCM_m_logr = LCM_m_r * LCM_r # [g/m3 / unit(log[um])] 
  LCM_m_logr_std_dev = LCM_m_r_std_dev * LCM_r # [g/m3 / unit(log[um])] 

  vterm_r = np.zeros(len(LCM_r))
  for i,r in enumerate(LCM_r):
    vterm_r[i] = beard77(r)
  LCM_m_logr_times_vterm = LCM_m_logr * vterm_r

  # read ensemble size
  LCM_ens = 0
  fs = open(pre+"setup.dat","r")
  rows = [x.split() for x in fs.readlines()]
  for row in rows:
    if row[0] == "n_rep":
      LCM_ens = float(row[2])

  # calculate error of the mean
  LCM_m_logr_err = LCM_m_logr_std_dev / np.sqrt(LCM_ens)

  return(LCM_r, LCM_m_logr, LCM_m_logr_std_dev, LCM_m_logr_times_vterm, LCM_m_logr_err, vterm_r)



def plot_DSD(data_labels, data_colors, data_ls, data_la, EFM_flag, fig, ax, time, STD_flag=True, xlim=[1,1e3], ylim=[1e-4,1e1]):
  if(STD_flag):
    ax0=ax[0]
    ax1=ax[1]
  else:
    ax0=ax
  for pre in data_labels:
    LCM_r, LCM_m_logr, LCM_m_logr_std_dev, LCM_m_logr_times_vterm, LCM_m_logr_err, vterm_r = read_DSD(pre, time)
    LCM_dlogr = np.log(LCM_r[1]) - np.log(LCM_r[0]) # [log(um)]
    LCM_M = LCM_m_logr * LCM_dlogr * volume * 1e-3 # mass of droplets of this size in the cell [kg]
    LCM_N = LCM_M / four_over_three_pi_rhow / np.power(LCM_r,3) # number of droplets of this size in the cell [1]
    print(pre)
    print("total number of droplets in the cell in LCM at the end: ", np.sum(LCM_N))
    print("total mass of droplets in the cell in LCM at the end: ", np.sum(LCM_M))
    print("0-th moment of number density (total droplet number): ", np.sum(LCM_N * np.power(LCM_r,0)))
    print("3-rd moment of number density (total droplet mass): ", np.sum(LCM_N * np.power(four_over_three_pi_rhow * LCM_r,3)))
    print("6-th moment of number density (radar reflectivity): ", np.sum(LCM_N * np.power(four_over_three_pi_rhow * LCM_r,6)))
    print("9-th moment of number density (for comparison with Unterstrasser 2017,2020): ", np.sum(LCM_N * np.power(four_over_three_pi_rhow * LCM_r,9)))
  
  # plot m(log r)
    ax0.plot(LCM_r * 1e6, LCM_m_logr, color=data_colors[pre], ls=data_ls[pre], alpha=data_la[pre], label= data_labels[pre])
    if(STD_flag):
      ax1.plot(LCM_r * 1e6, LCM_m_logr_std_dev, color=data_colors[pre], ls=data_ls[pre], alpha=data_la[pre], label= data_labels[pre])
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
  if EFM_flag:
    efm_data = np.genfromtxt("/home/piotr/praca/coal_fluctu_dim/LCM_DSD_fluctuations/data/EFM/dt 0.1_rq015.0_xmw 2.0_scal10.0_isw3_cut_____rmr_____tmax301._boplot00.out", unpack=True)
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
  
    ax0.plot(EFM_r[:int(EFM_N.size/2)] * 1e6, EFM_m_logr * 1e3, ls='--', color='black', label="SCE") # radius in [um], mass density in [g/m3 / unit(log[um])]
    if time > 0 and STD_flag: # dont plot SCE initial std dev - it is zero...
      ax1.plot(EFM_r[:int(EFM_N.size/2)] * 1e6, EFM_m_logr_std_dev * 1e3, ls='--', color='black', label="SCE") # radius in [um], mass density in [g/m3 / unit(log[um])]
    
    # plot N (number of droplets in the bin)
    #ax[0].plot(EFM_r * 1e6, EFM_N, label='end ' + data_labels[pre])
    
    #  ax[0].plot(EFM_r * 1e6, efm_data[1] / efm_data[0] * 1e6 * 1e3, label="EFM") 
    #  ax[1].plot(EFM_r * 1e6, np.sqrt(four_over_three_pi_rhow * np.power(efm_data[0] * 1e-6,3) * EFM_M), label="EFM") # std_dev of total mass of droplets of this size in the whole cell [kg]
    #  ax[0].plot(efm_data[0], efm_data[1] * EFM_dlogr / four_over_three_pi_rhow / np.power(efm_data[0] * 1e-6,3), label="EFM") # EFM data is mass density m(ln r) [kg/m3] with radius in microns; LCM data is m(r) with radius in meters
    
    print("total number of droplets in the cell in EFM at the end: ", np.sum(EFM_N[int(EFM_N.size/2):])) # EFM_N contains init spectrum and final spectrum
    print("total mass of droplets in the cell in EFM at the end: ", np.sum(EFM_M[int(EFM_M.size/2):]))
  
  #plt.title('Time at which largest droplet exceeds 3 mm radius and size of the droplet')
  ax0.set_xlabel('$r\ [\mathrm{\mu m}]$')
  #ax0.set_ylabel('mean mass density m(r) [g/m^3 / m]')
#  ax0.set_ylabel('mean mass density m(log r) [g/m^3 / unit(log[um])')
  ax0.set_ylabel('$<m>\ [\mathrm{gm^{-3} \ / \ unit(ln(\mu m))}]$')
#  ax1.set_ylabel('standard deviation of mass density m(log r) [g/m^3 / unit(log[um])')
  
  ax0.set_xscale('log')
  ax0.set_yscale('log')
  ax0.set_xlim(xlim)
  ax0.set_ylim(ylim)

  if(STD_flag):
    ax1.set_xlabel('$r\ [\mathrm{\mu m}]$')
    ax1.set_ylabel('$\sigma(m)\ [\mathrm{gm^{-3} \ / \ unit(ln(\mu m))}]$')
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlim(xlim)
    ax1.set_ylim(ylim)
  
  #ax[0].set_xlim(40,)
  #ax[0].set_ylim(0,0.30)
  
#  plt.grid(axis='y', color='0.5')
#  ax[0].legend()
#  ax[1].legend()
#  fig.savefig("/home/piotr/praca/coal_fluctu_dim/LCM_DSD_fluctuations/img/DSD_"+str(outname)+"_t"+str(time)+".pdf")
  #plt.show()


def plot_DSD_diff(ax, data_labels, data_colors, data_ls, data_la, ref_label, time, outname):
#  fig, ax = plt.subplots(figsize=(8,8))

  key_list = list(data_labels.keys())
  val_list = list(data_labels.values())
  pre_ref = key_list[val_list.index(ref_label)]

  # read reference data
  LCM_r, LCM_m_logr_ref, LCM_m_logr_std_dev_ref, LCM_m_logr_times_vterm_ref, LCM_m_logr_err_ref, vterm_r = read_DSD(pre_ref, time)

  # plot diff from reference data
  for pre in data_labels:
    if pre == pre_ref:
      continue

    LCM_r, LCM_m_logr, LCM_m_logr_std_dev, LCM_m_logr_times_vterm, LCM_m_logr_err, vterm_r = read_DSD(pre, time)

    LCM_m_logr_diff_err = np.sqrt(pow(LCM_m_logr_err,2) + pow(LCM_m_logr_err_ref,2))
    LCM_m_logr_times_vterm_diff_err = LCM_m_logr_diff_err * vterm_r # assume vterm has no error 
#    plt.plot(LCM_r * 1e6, LCM_m_logr_times_vterm - LCM_m_logr_times_vterm_ref, color=data_colors[pre], label= data_labels[pre])
    ax.errorbar(LCM_r * 1e6, LCM_m_logr_times_vterm - LCM_m_logr_times_vterm_ref, yerr = SE_scale * LCM_m_logr_times_vterm_diff_err, color=data_colors[pre], ls=data_ls[pre], alpha=data_la[pre], label= data_labels[pre])

#    plt.plot(LCM_r * 1e6, LCM_m_logr - LCM_m_logr_ref, color=data_colors[pre], label= data_labels[pre])
  ax.axhline(y=0., color='grey', linestyle='dotted')

  
  ax.set_xlabel('$r\ [\mathrm{\mu m}]$')
  ax.set_ylabel('$(<m>_\mathrm{SD} - <m>_\mathrm{0})\ v_\mathrm{t}\ \ [\mathrm{g m^{-2} s^{-1} \ / \ unit(ln(um))}]$')
  ax.set_xscale('log')
#  ax.set_yscale('symlog', linthresh=1e-3)
  ax.set_xlim(7,3e2)
  ax.set_ylim()
#  ax.grid(axis='y', color='0.5')
  ax.legend(loc='upper left')
#  plt.savefig("/home/piotr/praca/coal_fluctu_dim/LCM_DSD_fluctuations/img/DSD_diff_"+str(outname)+"_t"+str(time)+".pdf")
