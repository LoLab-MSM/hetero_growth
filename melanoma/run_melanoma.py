from pysb import *
from pysb.bng import run_ssa
from melanoma import model
import numpy as np
import scipy.stats as ss
from scipy.special import erf
import matplotlib.pyplot as plt
import scipy.stats as ss

def plot_timecourses(label, data, exclude=[]):
    
    plt.figure(label+"_expdata")
    slopes = []
    count = np.array(data[:,0])
    day = data[0:21,1]
    timecourses = np.zeros( ( len(day) , 1+len(count)/len(day) ) )
    timecourses[:,0] = day
#     t = timecourses[8:,0]
    t = timecourses[:,0]
    initial_pops = []
    for i in range(1,len(count)/len(day)):
        if i not in exclude:
            timecourses[:,i] = count[(i-1)*len(day):i*len(day)]
#             tc = np.array([p for p in timecourses[8:,i] if p > 0.])
            tc = np.array([p for p in timecourses[:,i] if p > 0.])
            initial_pops.append(tc[0])
            plt.plot(t[:len(tc)], tc, linewidth=2)
            try:
                slope, intercept, r_value, p_value, std_err = ss.linregress(t[:len(tc)],np.log(tc))
            except RuntimeWarning:
                pass
            else:
                slopes.append(slope)
    plt.yscale('log', basey=2)
    plt.xlabel('Time post drug (days)', fontsize=18)
    plt.ylabel('Cell count', fontsize=18)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    return initial_pops, slopes

def plot_slopes(label, slopes):
    
    plt.figure(label+"_expslopes") #"Distribution of slopes (kdiv=%6.4f, DIP=%6.4f)" % (kdiv, kdiv-kdeg))
        
    mean = np.average(slopes)
    sdev = np.std(slopes)
    skew = ss.skew(slopes)
    print "mean:", mean
    print "sdev:", sdev
    print "skew:", skew
    print
    
    # Histogram
    hist, bin_edges = np.histogram(slopes, bins=10, density=True)
    bin_width = bin_edges[1]-bin_edges[0]
    plt.bar(bin_edges[:-1], hist, width=bin_width)
    
    # Skew-normal fit
    xmin = np.floor(bin_edges[0]*10)/10
    xmax = np.ceil(bin_edges[-1]*10)/10
    xpoints = np.linspace(xmin,xmax,200)
    delta2 = 0.5*np.pi*np.abs(skew)**(2/3) / ( np.abs(skew)**(2/3) + (0.5*(4-np.pi))**(2/3) )
    delta = np.sign(skew)*np.sqrt(delta2)
    alpha = delta/np.sqrt(1-delta2)
    omega = sdev/np.sqrt(1-2*delta2/np.pi)
    xi = mean-omega*delta*np.sqrt(2/np.pi)
    snorm = np.array([1./omega/np.sqrt(2.*np.pi)*np.exp(-0.5*((x-xi)/omega)**2)*(1.+erf(alpha*(x-xi)/np.sqrt(2.)/omega)) for x in xpoints])
    
    plt.plot(xpoints, snorm, 'r', linewidth=5)
    plt.plot([mean,mean],[0,np.max(hist)],'r--', linewidth=5)
    plt.xlim(xmin=min(xpoints),xmax=max(xpoints))
    # plt.ylim(ymax=100)
    plt.xlabel('Slope (ln(cell count)/day)', fontsize=18)
    plt.ylabel('Prob. density', fontsize=18)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)

data = np.genfromtxt('/Users/lopezlab/Documents/workspace/Melanoma/SKMEL5-2uMcFP data.csv', delimiter=",", skip_header=1, usecols=(2,9))
initial_pops_2uM, slopes = plot_timecourses('SKMEL5-2uMcFP', data, exclude=[])#100,93,88,81,63,56,30,18])
plot_slopes('SKMEL5-2uMcFP', slopes)

data = np.genfromtxt('/Users/lopezlab/Documents/workspace/Melanoma/SKMEL5-8uMcFP data.csv', delimiter=",", skip_header=1, usecols=(2,9))
initial_pops_8uM, slopes = plot_timecourses('SKMEL5-8uMcFP', data, exclude=[])#104,102,80,77,52,43,8,2])
plot_slopes('SKMEL5-8uMcFP', slopes)
 
data = np.genfromtxt('/Users/lopezlab/Documents/workspace/Melanoma/SKMEL5-16uMcFP data.csv', delimiter=",", skip_header=1, usecols=(2,9))
initial_pops_16uM, slopes = plot_timecourses('SKMEL5-16uMcFP', data, exclude=[])#81,79,77,76,74,73,65,64,58,46,42,41,40,39,38,36,35,33,32,31,30,29,27,26,25,23,22,21,20,17,12,11,10,6,3,2,1])
plot_slopes('SKMEL5-16uMcFP', slopes)

# plt.show()
# quit()

tspan = 24.*np.array([  7.6926375,    8.92476667,   9.73106667,  10.7534125,   10.99522083,
         11.78024167,  12.73057917,  12.9729875,   13.70632917,  14.875775,
         16.12682917,  16.7051625,   17.764825  ])

plt.figure("SKMEL5-2uMcFP_simulations")
n_sims = len(initial_pops_2uM)
seeds = [500+i for i in range(n_sims)]
slopes_sim = []
for i in range(n_sims):
    print i
    model.parameters['Cell_0'].value = initial_pops_2uM[i]
    x = run_ssa(model, tspan, seed=seeds[i], verbose=False)
    tc = np.array([p for p in x['Cell_total'] if p > 0.])
    try:
        slope, intercept, r_value, p_value, std_err = ss.linregress(x['time'][:len(tc)],np.log(tc))
    except RuntimeWarning:
        pass
    else:
        slopes_sim.append(slope)
    plt.plot(x['time']/24., x['Cell_total'], linewidth=2)
    plt.yscale('log', basey=2)
    plt.xlabel('Time (days)', fontsize=18)
    plt.ylabel('Cell count', fontsize=18)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)

plt.figure("SKMEL5-2uMcFP_simslopes") 
    
mean = np.average(slopes_sim)
sdev = np.std(slopes_sim)
skew = ss.skew(slopes_sim)
print "mean:", mean
print "sdev:", sdev
print "skew:", skew

# Histogram
hist, bin_edges = np.histogram(slopes_sim, bins=10, density=True)
bin_width = bin_edges[1]-bin_edges[0]
plt.bar(bin_edges[:-1], hist, width=bin_width)

plt.show()



