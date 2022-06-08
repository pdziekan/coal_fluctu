import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from math import floor

SE_scale = 1.96 # used to scale standard error, 1.96 gives 95% confidence interval

def plot_coal_series(plot, data, fig, axs):
#  plt.rcParams['text.usetex'] = True # latex labels
#  fig, ax = plt.subplots(figsize=(12,9))

  labels = []
  for pre in data:
    label = data[pre]
    if label not in labels:
      labels.append(label)

  # colors fom the default 
  prop_cycle = plt.rcParams['axes.prop_cycle']
  colors = prop_cycle.by_key()['color']
  color_lab = {}
  for i, lab in enumerate(labels):
    color_lab[lab] = colors[i]

  # y axis labels
  y_lab = {}
  y_lab["tau"] = "\\theta"
  y_lab["nrain"] = "N_r"
  y_lab["rmax"] = "r_\mathrm{max}"

  
  #fig, ax = plt.subplots(figsize=(12,9))
  
  #  time_min = {} 
  time_max = 1e10     # do not plot results for time exceeding the moment fastest simulation reached the end (produced large droplet)
  all_time_lab = {}
  all_data_lab = {}
  ensemble_lab = {}
  bin_mean_lab = {}
  bin_mean_err_lab = {}
  bin_std_dev_lab = {}
  bin_std_dev_err_lab = {}
  bin_mean_time_lab = {}
  
  # plot each simulation
  for pre in data:
    label = data[pre]
    print(plot)
    print(pre)
    print(label)
  
    ft = open(pre+"time.dat","r")
    fd = open(pre+plot+".dat","r")
    all_time=[x.split() for x in ft.readlines()]
    all_data=[x.split() for x in fd.readlines()]
  
    if label not in ensemble_lab:
      ensemble_lab[label] = 0
  
    for idx, (row_t, row_d) in enumerate(zip(all_time, all_data)):
      time = np.array(row_t, dtype=np.float32)
      series = np.array(row_d, dtype=np.float32)
      if plot == "rmax":
        series *= 1e6 # m -> um
  #      ax.plot(time, series[0:len(time)], label=data[pre]+'('+str(idx)+')')
  #      time_min[label] = min(time_min[label], np.min(time)) if label in time_min else np.min(time)
      time_max            = min(time_max, time[-1])
      all_time_lab[label] = np.append(all_time_lab[label], time) if label in all_time_lab else time
      all_data_lab[label] = np.append(all_data_lab[label], series) if label in all_data_lab else series
      ensemble_lab[label] += 1
      #print(all_time_lab)
    print(time_max)
    #print(time_min[label], time_max[label], data_min[label], data_max[label])
  
    #plt.xlim([0, 800])
    
    #plt.legend()
    #fig.savefig(outname+"_series_"+plot+".png")
  plt.close()
  
  
  # plot statistics from ensembles of same simulations - needs binning in time space, because of variable time step
  #time_bin_edges = {}
  #data_binned = {}
  #  print(shape(all_time_lab[label]))

  outfreq = 10 # [s], output frequency in coal fluctu dim

  time_max = floor(time_max)
  time_max -= time_max % outfreq
#    time_max -= 10*outfreq # just to be safe?
  print(time_max)
  #  nbins_ = [1000]
  nbins_ = [int(time_max / outfreq)]
  _bin_edges = np.arange(0,time_max+outfreq,outfreq)
  _bin_edges[-1] -= 1e-10 # last bin is closed on the right, but we dont want t==time_max

    
  print(nbins_)
  for nbins in nbins_:
#      fig, axs = plt.subplots(2,1)
    for label in labels:
    #    label = data[pre]

    #    bin_count,   bin_edges, binnumber = stats.binned_statistic(all_time_lab[label], all_data_lab[label], bins=nbins, range=(0, time_max), statistic='count')
    #    bin_mean   , bin_edges, binnumber = stats.binned_statistic(all_time_lab[label], all_data_lab[label], bins=nbins, range=(0, time_max), statistic='mean')
    #    bin_std_dev, bin_edges, binnumber = stats.binned_statistic(all_time_lab[label], all_data_lab[label], bins=nbins, range=(0, time_max), statistic='std')

      time_till_tmax = all_time_lab[label][all_time_lab[label]<=time_max].copy()
      data_till_tmax = all_data_lab[label][all_time_lab[label]<=time_max].copy()
#      print("max of time till tmax:", np.max(time_till_tmax))
#      print("max of all time lab:", np.max(all_time_lab[label]))
#      print("max of data till tmax:", np.max(data_till_tmax))
#      print("max of all data lab:", np.max(all_data_lab[label]))
#      print("bin edges:", _bin_edges)
#
#      print("time between (tmax-outfreq, tmax>:", time_till_tmax[time_till_tmax>time_max-outfreq])
#      print("number of time points between (tmax-outfreq, tmax>:", np.size(time_till_tmax[time_till_tmax>time_max-outfreq]))
#      print("data between (tmax-outfreq, tmax>:", data_till_tmax[time_till_tmax>time_max-outfreq])
#      print("number of data points between (tmax-outfreq, tmax>:", np.size(data_till_tmax[time_till_tmax>time_max-outfreq]))

      bin_count,   bin_edges, binnumber = stats.binned_statistic(time_till_tmax, data_till_tmax, bins=_bin_edges, statistic='count')
      bin_mean   , bin_edges, binnumber = stats.binned_statistic(time_till_tmax, data_till_tmax, bins=_bin_edges, statistic='mean')
      bin_std_dev, bin_edges, binnumber = stats.binned_statistic(time_till_tmax, data_till_tmax, bins=_bin_edges, statistic='std')

      bin_mean_time,   bin_edges, binnumber = stats.binned_statistic(time_till_tmax, time_till_tmax, bins=_bin_edges, statistic='mean')

#      print("time in the last bin:", time_till_tmax[binnumber==nbins])
#      print("data in the last bin:", data_till_tmax[binnumber==nbins])

      bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.
      bin_width = bin_edges[1] - bin_edges[0]
        #plt.hlines(bin_stats, bin_edges[:-1], bin_edges[1:], colors='g', lw=5, label='binned statistic of data')
    #    plt.plot((bin_edges[:-1] + bin_edges[1:]) / 2., bin_stats, label=label)# colors='g', lw=5, label='binned statistic of data')

#        axs[0].scatter(all_time_lab[label], all_data_lab[label], color=color_lab[label], s=2, alpha=0.15)

      bin_mean_time_lab[label] = bin_mean_time[bin_count>0]
      bin_mean_lab[label] = bin_mean[bin_count>0]
      bin_mean_err_lab[label] = bin_std_dev[bin_count>0] / np.sqrt(bin_count[bin_count>0])
      bin_std_dev_lab[label] = bin_std_dev[bin_count>0]
      bin_std_dev_err_lab[label] = bin_std_dev[bin_count>0] / np.sqrt(2*(bin_count[bin_count>0]-1))
  
      axs[0].errorbar(bin_mean_time_lab[label], bin_mean_lab[label], yerr=SE_scale * bin_mean_err_lab[label],       label=label, color=color_lab[label], ls='', marker='.', markersize=4)
      axs[1].errorbar(bin_mean_time_lab[label], bin_std_dev_lab[label], yerr=SE_scale * bin_std_dev_err_lab[label], label=label, color=color_lab[label], ls='', marker='.', markersize=4)
  
#        axs[0].errorbar(bin_centers[bin_count>0], bin_mean[bin_count>0]   , yerr=1.96 * bin_std_dev[bin_count>0] / np.sqrt(bin_count[bin_count>0]),       label=label, color=color_lab[label], ls='')
#        axs[1].errorbar(bin_centers[bin_count>0], bin_std_dev[bin_count>0], yerr=1.96 * bin_std_dev[bin_count>0] / np.sqrt(2*(bin_count[bin_count>0]-1)), label=label, color=color_lab[label], ls='')
#        axs[0].hlines(bin_mean[bin_count>0],    bin_edges[:-1][bin_count>0], bin_edges[1:][bin_count>0], color=color_lab[label])
#        axs[1].hlines(bin_std_dev[bin_count>0], bin_edges[:-1][bin_count>0], bin_edges[1:][bin_count>0], color=color_lab[label])
  
#      print("time max: " + str(time_max))
#      print("ensemble size: " + str(ensemble_lab[label]))
#      print("numebr of points in bins:")
#      print(bin_count)
#      print("mean value in bins:")
#      print(bin_mean)
#      print("std dev in bins:")
#      print(bin_std_dev)

#      print("mean time in bins:", bin_mean_time)
  
  #  plt.legend()
    axs[0].set_ylabel('$<'+y_lab[plot]+'>$')
    axs[1].set_ylabel('$\sigma ('+y_lab[plot]+')$')
    plt.legend()

#    fig.savefig(outname+"_series_"+plot+"_stats_nbins_"+str(nbins)+".pdf")
 #   fig.savefig(outname+"_series_"+plot+"_stats_nbins_"+str(nbins)+".svg")
    #plt.show()
#    axs[0].clear()
#    axs[1].clear()
#      fig.clear()

  return bin_mean_time_lab, bin_mean_lab, bin_mean_err_lab, bin_std_dev_lab, bin_std_dev_err_lab
      
    #  time_bin_edges[label] = np.linspace(time_min[label], time_max[label], bin_count+1)
    #  data_binned[label] = np.zeros(bin_count)
    #  ft = open(pre+"time.dat","r")
    #  fd = open(pre+plot+".dat","r")
    #  all_time=[x.split() for x in ft.readlines()]
    #  all_data=[x.split() for x in fd.readlines()]
    #  for idx, (row_t, row_d) in enumerate(zip(all_time, all_data)):
    #    time = np.array(row_t, dtype=np.float32)
    #    series = np.array(row_d, dtype=np.float32)
    #    ax.plot(time, series[0:len(time)], label=data[pre]+'('+str(idx)+')')



def plot_coal_series_diff(labelA, labelB, time, mean, mean_err, stddev, stddev_err, fig, axs):
#  print(time[labelA], time[labelB])
  avg_time = (time[labelA] + time[labelB]) / 2.
  #print("avg time:",avg_time)
  avg_time_err = abs(time[labelA] - time[labelB]) / 2.

  mean_diff = mean[labelB] - mean[labelA]
  mean_diff_err = np.sqrt(pow(mean_err[labelA], 2) + pow(mean_err[labelB], 2))
  stddev_diff = stddev[labelB] - stddev[labelA]
  stddev_diff_err = np.sqrt(pow(stddev_err[labelA], 2) + pow(stddev_err[labelB], 2))

#  print("mean diff:", mean_diff)
#  print("stddev diff:", stddev_diff)

  # note that scaling coefficients are not included in error calculation!
  # /100 to get percent
  mean_scale = mean[labelA] / 100.
  stddev_scale = stddev[labelA] / 100.
  mean_scale = 1
#  stddev_scale = 1

  scaled_mean_diff       = np.where(mean_scale != 0, mean_diff / mean_scale, mean_diff)
  scaled_mean_diff_err   = np.where(mean_scale != 0, mean_diff_err / mean_scale, mean_diff_err)
  scaled_stddev_diff     = np.where(stddev_scale != 0, stddev_diff / stddev_scale, stddev_diff)
  scaled_stddev_diff_err = np.where(stddev_scale != 0, stddev_diff_err / stddev_scale, stddev_diff_err)
  
  axs[2].errorbar(avg_time, scaled_mean_diff, xerr = avg_time_err, yerr=SE_scale * scaled_mean_diff_err, ls='', marker='.', markersize=4)
  axs[2].axhline(y=0., color='r', linestyle='--')
  axs[3].errorbar(avg_time, scaled_stddev_diff, xerr = avg_time_err, yerr=SE_scale * scaled_stddev_diff_err, ls='', marker='.', markersize=4)
  axs[3].axhline(y=0., color='r', linestyle='--')
