import numpy as np
import scipy.stats as ss
import snorm
import matplotlib.pyplot as plt

def rebound(cellline, erl=True, abt=False, sb=False, an=False, il=False, tnf=False, birb=False):

    # PC9.par+erl: Wells D02 and D03 not valid (don't know why)
    # HCC4006+erl
    # H1118+erl: Transient ends at day 5 (rather than day 3)
    # PC9.par+erl+ABT
    # PC9.par+erl+TNF

    data = np.genfromtxt('/Users/lopezlab/git/r-code-peter_nature/Data used for figures (pulled by R)/Rebound combinations (try).csv', 
                         delimiter=',', dtype=None, names=True)
    
    exclude_wells = []
    if cellline == 'PC9.par' and erl and not any([abt,sb,an,il,tnf,birb]):
        exclude_wells = ['D02','D03']        

    if cellline == 'H1118': minday = 5
    else: minday = 3
    
    timeday   = [d['Timeday'] for d in data 
                if d['Cellline'] == cellline 
                and d['Timeday'] >= minday 
                and (d['Erl']  > 0) == erl 
                and (d['ABT']  > 0) == abt 
                and (d['SB']   > 0) == sb 
                and (d['An']   > 0) == an 
                and (d['IL']   > 0) == il 
                and (d['TNF']  > 0) == tnf 
                and (d['BIRB'] > 0) == birb
                and d['Well'] not in exclude_wells]
    
    cellcount = [float(d['Cellcount']) for d in data 
                if d['Cellline'] == cellline 
                and d['Timeday'] >= minday 
                and (d['Erl']  > 0) == erl 
                and (d['ABT']  > 0) == abt 
                and (d['SB']   > 0) == sb 
                and (d['An']   > 0) == an 
                and (d['IL']   > 0) == il 
                and (d['TNF']  > 0) == tnf 
                and (d['BIRB'] > 0) == birb
                and d['Well'] not in exclude_wells]
    
    exp_times = [[timeday[0]]]
    exp_pops = [[cellcount[0]]]
    for i in range(1,len(timeday)):
        if timeday[i] > timeday[i-1]:
            exp_times[-1].append(timeday[i])
            exp_pops[-1].append(cellcount[i])
        else:
            exp_times.append([timeday[i]])
            exp_pops.append([cellcount[i]])

    timecourses = {}
    for i in range(len(exp_times)):
        timecourses[i] = {}
        timecourses[i]['Time'] = np.array(exp_times[i])
        timecourses[i]['CellCount'] = np.array(exp_pops[i])
    
    return timecourses

def colony_tracking(cellline, erl=True, abt=False, sb=False, an=False, il=False, tnf=False, chx=False, fsk=False, trm=False, _17a=False):

    if   cellline == 'DS9': cellline = 'PC9.C03'
    elif cellline == 'DS6': cellline = 'PC9.DS6'
    elif cellline == 'DS8': cellline = 'PC9.B03'

    if cellline == 'H1118': minday = 5
    else: minday = 3
    
    if erl and not any([abt,sb,an,il,tnf,chx,fsk,trm,_17a]):
        data = np.genfromtxt('/Users/lopezlab/git/r-code-peter_nature/Data used for figures (pulled by R)/DCGA colony tracking_3.csv', 
                             delimiter=',', dtype=None, names=True)
        
        timeday     =  [d['Time_Day'] for d in data 
                       if d['Cellline'] == cellline
                       and d['Time_Day'] >= minday 
                       and d['Condition'] == 'Erlotinib']
    
        cellcount   =  [float(d['CellCount']) for d in data 
                       if d['Cellline'] == cellline
                       and d['Time_Day'] >= minday 
                       and d['Condition'] == 'Erlotinib']
    else:
        data = np.genfromtxt('/Users/lopezlab/git/r-code-peter_nature/Data used for figures (pulled by R)/AllCombos.csv', 
                             delimiter=',', dtype=None, names=True)
        
        timeday   = [d['Time_Day'] for d in data 
                    if d['Cellline'] == cellline 
                    and d['Time_Day'] >= minday 
                    and (d['Erl']  > 0) == erl 
                    and (d['ABT']  > 0) == abt 
                    and (d['SB']   > 0) == sb 
                    and (d['An']   > 0) == an 
                    and (d['IL']   > 0) == il 
                    and (d['TNF']  > 0) == tnf 
                    and (d['CHX']  > 0) == chx
                    and (d['FSK']  > 0) == fsk
                    and (d['TRM']  > 0) == trm
                    and (d['17A']  > 0) == _17a]

        cellcount = [float(d['CellCount']) for d in data 
                    if d['Cellline'] == cellline 
                    and d['Time_Day'] >= minday 
                    and (d['Erl']  > 0) == erl 
                    and (d['ABT']  > 0) == abt 
                    and (d['SB']   > 0) == sb 
                    and (d['An']   > 0) == an 
                    and (d['IL']   > 0) == il 
                    and (d['TNF']  > 0) == tnf 
                    and (d['CHX']  > 0) == chx
                    and (d['FSK']  > 0) == fsk
                    and (d['TRM']  > 0) == trm
                    and (d['17A']  > 0) == _17a]
    
    exp_times = [[timeday[0]]]
    exp_pops = [[cellcount[0]]]
    for i in range(1,len(timeday)):
        if timeday[i] > timeday[i-1]:
            exp_times[-1].append(timeday[i])
            exp_pops[-1].append(cellcount[i])
        else:
            exp_times.append([timeday[i]])
            exp_pops.append([cellcount[i]])

    timecourses = {}
    for i in range(len(exp_times)):
        timecourses[i] = {}
        timecourses[i]['Time'] = np.array(exp_times[i])
        timecourses[i]['CellCount'] = np.array(exp_pops[i])
        
    exp_slopes = []
    for i in range(len(timecourses)):
        tc = np.array([p for p in timecourses[i]['CellCount'] if p > 0.])
        t = timecourses[i]['Time'][:len(tc)]
        try:
            slope, intercept, r_value, p_value, std_err = ss.linregress(t*24., np.log(tc))
        except RuntimeWarning:
            pass
        else:
            exp_slopes.append(slope)
    
    return timecourses, exp_slopes

def plot_timecourses(label, timecourses, xmin=None, xmax=None, ymin=None, ymax=None):
    
    plt.figure(label)
    
    for tc in timecourses:
        cellcount = [p for p in timecourses[tc]['CellCount'] if p > 0.]
        time = timecourses[tc]['Time'][:len(cellcount)]
        plt.plot(time-time[0], np.log2(cellcount/cellcount[0]), '0.5')
    
    plt.xlabel('Time post drug (days)', fontsize=18)
    plt.ylabel('log2(Normalized cell count)', fontsize=18)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.xlim(xmax=xmax, xmin=xmin)
    plt.ylim(ymin=ymin, ymax=ymax)

def plot_init_pops(label, timecourses, bins, xmin=None, xmax=None):
    
    plt.figure(label)
    init_pops = np.array([timecourses[i]['CellCount'][0] for i in timecourses])
    plot_hist(np.log10(init_pops), bins, xmin, xmax, density=False)
    plt.xlabel('log10(Initial cell count)', fontsize=18)
    plt.ylabel('Number', fontsize=18)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    
def plot_slopes(label, slopes, bins, xmin=None, xmax=None, skewnorm=True, fit=True, textcolor='k'):
    
    plt.figure(label)
    plot_hist(slopes, bins, xmin, xmax, density=True)
    plt.xlabel('Slope (ln(cell count)/h)', fontsize=18)
    plt.ylabel('Prob. density', fontsize=18)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    
    # Fit
    if fit:
        mean = np.average(slopes)
        sdev = np.std(slopes)
        skew = ss.skew(slopes)
        xpoints = np.linspace(xmin,xmax,10*bins)
        if skewnorm:
            pdf = snorm.pdf(xpoints,mean,sdev,skew)
            plt.annotate("skew: %.3g" % skew, (0.05, 0.74), xycoords='axes fraction', fontsize=18, color=textcolor)
        else:
            pdf = ss.norm.pdf(xpoints,mean,sdev)
        plt.plot(xpoints, pdf, 'r', linewidth=5)
        plt.annotate("mean: %.3g" % mean, (0.05, 0.90), xycoords='axes fraction', fontsize=18, color=textcolor)
        plt.annotate("sdev: %.3g" % sdev, (0.05, 0.82), xycoords='axes fraction', fontsize=18, color=textcolor)
        
def plot_hist(data, bins, xmin=None, xmax=None, density=False):
    
    if xmin: rmin = xmin
    else: rmin = min(data)
    if xmax: rmax = xmax
    else: rmax = max(data)
    hist, bin_edges = np.histogram(data, bins=bins, density=density, range=(rmin,rmax))
    bin_width = bin_edges[1]-bin_edges[0]
    plt.bar(bin_edges[:-1], hist, width=bin_width, color='0.5')
    plt.xlim(xmin, xmax)

#########################################################

# cellline = 'PC9.par'
# condition = 'erl+TNF'
# tc = rebound(cellline, abt=False, tnf=True)
# plot_timecourses('Timecourses: %s+%s (rebound)' % (cellline, condition), tc, xmax=35, ymin=-3, ymax=1)
# plt.plot([0,35],[0,0],'r',linewidth=3)
# plot_init_pops('Initial populations: %s+%s (rebound)' % (cellline, condition), tc, bins=20, xmin=3, xmax=5)
#  
# cellline = 'PC9.par'
# condition = 'erl'
# tc, slopes = colony_tracking(cellline, abt=False, tnf=False)
# plot_timecourses('Timecourses: %s+%s (cFP)' % (cellline, condition), tc)
# plot_init_pops('Initial populations: %s+%s (cFP)' % (cellline, condition), tc, bins=30, xmin=0, xmax=3)
# plot_slopes('DIP rate distribution: %s+%s (cFP)' % (cellline, condition), slopes, bins=20, xmin=-0.05, xmax=0.05, skewnorm=True)
# plt.ylim(ymax=150)
#  
# plt.show()

'''
# Cell counts at 72 hr (initial point for linear fit)
plt.figure()
init_pops = np.log10(timecourses[3,1:]) #np.log10(np.loadtxt('/Users/lopezlab/Documents/workspace/HeteroGrowth/expt_data/DS9_erl_init_pops.txt'))
hist, bin_edges = np.histogram(init_pops, bins=30, density=True, range=(0,3))
bin_width = bin_edges[1]-bin_edges[0]
plt.bar(bin_edges[:-1], hist, width=bin_width)

# lognormal fit
mean = np.average(init_pops)
sdev = np.std(init_pops)
xpoints = np.linspace(0,3,300)
lognorm = ss.norm.pdf(xpoints,mean,sdev) 
#####
# skew = ss.skew(init_pops)
# delta2 = 0.5*np.pi*np.abs(skew)**(2/3) / ( np.abs(skew)**(2/3) + (0.5*(4-np.pi))**(2/3) )
# delta = np.sign(skew)*np.sqrt(delta2)
# alpha = delta/np.sqrt(1-delta2)
# omega = sdev/np.sqrt(1-2*delta2/np.pi)
# xi = mean-omega*delta*np.sqrt(2/np.pi)
# slognorm = np.array([1./omega/np.sqrt(2.*np.pi)*np.exp(-0.5*((x-xi)/omega)**2)*(1.+erf(alpha*(x-xi)/np.sqrt(2.)/omega)) for x in xpoints])
#####
plt.plot(xpoints, lognorm, 'r', linewidth=5)

plt.ylabel("Prob. density", fontsize=18)
plt.xlabel("log10(Cell count at 72 h post drug)", fontsize=18)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)

# Distribution of linear best-fit slopes
plt.figure()
slopes = np.loadtxt('/Users/lopezlab/Documents/workspace/HeteroGrowth/expt_data/PC9par_erl_slopes.txt')
hist, bin_edges = np.histogram(slopes, range=(-0.05,0.05), bins=40, density=True)
bin_width = bin_edges[1]-bin_edges[0]
plt.bar(bin_edges[:-1], hist, width=bin_width)

print "mean:", np.average(slopes)
print "sdev:", np.std(slopes)
print "skew:", ss.skew(slopes)

# Skew-normal fit
mean = np.average(slopes)
sdev = np.std(slopes)
skew = ss.skew(slopes)
xpoints = np.linspace(-0.05,0.05,200)
delta2 = 0.5*np.pi*np.abs(skew)**(2/3) / ( np.abs(skew)**(2/3) + (0.5*(4-np.pi))**(2/3) )
delta = np.sign(skew)*np.sqrt(delta2)
alpha = delta/np.sqrt(1-delta2)
omega = sdev/np.sqrt(1-2*delta2/np.pi)
xi = mean-omega*delta*np.sqrt(2/np.pi)
snorm = np.array([1./omega/np.sqrt(2.*np.pi)*np.exp(-0.5*((x-xi)/omega)**2)*(1.+erf(alpha*(x-xi)/np.sqrt(2.)/omega)) for x in xpoints])

plt.plot(xpoints, snorm, 'r', linewidth=5)
plt.plot([mean,mean],[0,np.max(hist)],'r--', linewidth=5)
plt.xlim(xmin=-0.05,xmax=0.02)
plt.ylim(ymax=100)
plt.xlabel('Slope (ln(cell count)/h)', fontsize=18)
plt.ylabel('Prob. density', fontsize=18)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)

# max_index = 0
# for i in range(1,len(snorm)):
#     if snorm[i] > snorm[max_index]:
#         max_index = i
# print xpoints[max_index], snorm[max_index]

plt.show()
'''