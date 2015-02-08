import hetero_growth as hg
from hetero_growth import model
from pysb.integrate import odesolve
from pysb.bng import run_ssa
from pysb.bng import set_bng_path
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as ss
import snorm
import read_expt_data as exp
import os

set_bng_path('/Users/lopezlab/git/bionetgen/bng2')
OUTDIR = '/Users/lopezlab/Documents/workspace/HeteroGrowth/RESULTS'
n_replicates = 1
envelope = False
run_ode = False
output_to_file = False

##################
condition = {'cellline':'PC9.par', 'abt':False, 'tnf':False}
xpop  = 1. #0.3 
xmean = 0.7 #0.42
xsdev = 0.8 #0.45 
xskew = 1.1 #-0.5
kdiv_base = 0.02 #0.5
nbins = 40
annot_mean = (0.65, 0.25)
annot_sdev = (0.65, 0.17)
annot_skew = (0.65, 0.09)
annot_DIPn = (0.05, 0.20)
annot_DIPp = (0.05, 0.12)
# FOR FIGS
# ---------
# n_replicates = 2
# seed = 150
##################
# condition = {'cellline':'HCC4006', 'abt':False, 'tnf':False}
# xpop  = 1.
# xmean = 1.6
# xsdev = 2.05 #2.
# xskew = 1.2
# kdiv_base = 0.05
# nbins = 40
# annot_mean = (0.65, 0.93)
# annot_sdev = (0.65, 0.85)
# annot_skew = (0.65, 0.77)
# annot_DIPn = (0.05, 0.90)
# annot_DIPp = (0.05, 0.83)
# FOR FIGS
# ---------
# n_replicates = 4
# seed = 150
##################
# condition = {'cellline':'H1118', 'abt':False, 'tnf':False}
# xpop  = 1.
# xmean = 0.25
# xsdev = 0.5
# xskew = 1.2
# kdiv_base = 0.05
# nbins = 40
# annot_mean = (0.65, 0.25)
# annot_sdev = (0.65, 0.17)
# annot_skew = (0.65, 0.09)
# annot_DIPn = (0.05, 0.20)
# annot_DIPp = (0.05, 0.12)
# FOR FIGS
# ---------
# n_replicates = 5
# seed = 150
##################
# condition = {'cellline':'PC9.par', 'abt':True, 'tnf':False}
# xpop  = 1.
# xmean = 0.9
# xsdev = 1.1
# xskew = 1.
# kdiv_base = 0.05
# nbins = 40
# annot_mean = (0.05, 0.92)
# annot_sdev = (0.05, 0.84)
# annot_skew = (0.05, 0.76)
# annot_DIPn = (0.02, 0.10)
# annot_DIPp = (0.02, 0.02)
# FOR FIGS
# ---------
# n_replicates = 5
# seed = 150
##################
# condition = {'cellline':'PC9.par', 'abt':False, 'tnf':True}
# xpop  = 0.5 #1.
# xmean = 0.9
# xsdev = 1.
# xskew = 0.8
# kdiv_base = 0.05
# nbins = 40
# annot_mean = (0.05, 0.92)
# annot_sdev = (0.05, 0.84)
# annot_skew = (0.05, 0.76)
# annot_DIPn = (0.05, 0.10)
# annot_DIPp = (0.05, 0.02)
# FOR FIGS
# ---------
# n_replicates = 6
# seed = 150
##################

OUTFILES = {}
if output_to_file:
    prefix = condition['cellline'].replace('.','') + '_erl'
    if condition['abt']: prefix += '_abt'
    if condition['tnf']: prefix += '_tnf'
    if envelope: OUTFILES['envelope']    = open(os.path.join(OUTDIR, prefix+'_envelope.csv'), 'w')
    else:        OUTFILES['timecourses'] = open(os.path.join(OUTDIR, prefix+'_timecourses.csv'), 'w')

def run_ssa_sims(n_sims, t, init_pops, seed=None, stop_mult=2, full_output=[], TRACK=None, verbose=False, **simargs):
    
    ref = {}
    for p in model.parameters:
        ref[p.name] = p.value

    if seed == None:
        seeds = np.random.randint(low=1,size=n_sims)
    else:
        seeds = [seed+i for i in range(n_sims)]
    timecourses = {}
    for i in range(n_sims):
        timecourses[i] = {}
        if TRACK:
            TRACK.write(str(i+1)+'\n')
        else:
            print i+1
        for j in range(len(init_pops[i])):
            if init_pops[i][j] > 0.:
                model.parameters['Cell%d_Init' % j].value = init_pops[i][j]
                model.parameters['kdiv_%d' % j].value = ref['kdiv_%d' % j]
                model.parameters['kdeg_%d' % j].value = ref['kdeg_%d' % j]
            else:
                model.parameters['Cell%d_Init' % j].value = 0.
                model.parameters['kdiv_%d' % j].value = 0.
                model.parameters['kdeg_%d' % j].value = 0.
        if stop_mult > 0:
            x = run_ssa(model, t, seed=seeds[i], stop_if='Cell_total>=%d*%d' % (stop_mult, sum(init_pops[i])), verbose=verbose, **simargs)
        else:
            x = run_ssa(model, t, seed=seeds[i], verbose=verbose, **simargs)
        timecourses[i]['Time'] = x['time']
        timecourses[i]['CellCount'] = x['Cell_total']
        if i in full_output:
            for key in x.dtype.names:
                if key not in ['time','Cell_total']:
                    timecourses[i][key] = x[key]

    # Restore kdiv and kdeg values
    for i in range(len(init_pops[0])):
        model.parameters['kdiv_%d' % i].value = ref['kdiv_%d' % i]
        model.parameters['kdeg_%d' % i].value = ref['kdeg_%d' % i]

    return timecourses
####################################################################################################

rb_tc = exp.rebound(**condition)
exp.plot_timecourses('Rebound', rb_tc, xmin=0, xmax=35, ymin=-3, ymax=1)
# plt.show()
# quit()
total_cells = np.array([rb_tc[i]['CellCount'][0] for i in rb_tc])
total_cells = np.round(xpop*total_cells)

##### For HCC4006
# total_cells = [i for i in total_cells if i > 5000.]
# print total_cells
# quit()
#####

cfp_tc, cfp_slopes = exp.colony_tracking(**condition)
exp.plot_slopes('DIP rate distribution', cfp_slopes, bins=nbins, xmin=-0.05, xmax=0.05, skewnorm=True, fit=True, textcolor='r')
mean = np.mean(cfp_slopes)
sdev = np.std(cfp_slopes) 
skew = ss.skew(cfp_slopes)
#####
mean *= xmean
sdev *= xsdev
skew *= xskew
plt.plot(np.linspace(-0.05,0.05,3*40), snorm.pdf(np.linspace(-0.05,0.05,3*40),mean,sdev,skew), 'b', linewidth=5)
plt.annotate("mean: %.3g" % mean, (0.65, 0.90), xycoords='axes fraction', fontsize=18, color='b')
plt.annotate("sdev: %.3g" % sdev, (0.65, 0.82), xycoords='axes fraction', fontsize=18, color='b')
plt.annotate("skew: %.3g" % skew, (0.65, 0.74), xycoords='axes fraction', fontsize=18, color='b')
# plt.show()
# quit()
#####
DIP = np.linspace(-0.05,0.05,nbins+1)
kdiv = np.array([kdiv_base]*len(DIP))
# kdiv = np.array([kdiv_base for d in DIP if d <= 0.] + [10.*kdiv_base for d in DIP if d > 0.])
# for i in range(len(DIP)):
#     print 's%d:' % i, DIP[i]
# quit()
kdeg = np.array([kdiv[i]-DIP[i] for i in range(len(DIP)) if kdiv[i]-DIP[i] > 0.])
kdiv = kdiv[:len(kdeg)]
DIP = DIP[:len(kdeg)]
# probs = ss.norm.pdf(DIP,mean,sdev)*(DIP[1]-DIP[0])
probs = snorm.pdf(DIP,mean,sdev,skew)*(DIP[1]-DIP[0])
probs /= sum(probs)
init_cells = []
np.random.seed(1)
for i in range(n_replicates):
    init_cells.append([np.random.multinomial(n, probs, size=1)[0] for n in total_cells])
init_cells = np.reshape(init_cells, (n_replicates*len(total_cells),len(DIP)))
init_cells = np.round(init_cells)

hg.define_model(kdiv, kdeg)
t = np.linspace(0, 35*24, 36)

# SSA simulations
n_sims = len(init_cells)
full_output = [] #[33,78]
tc = run_ssa_sims(n_sims, t, init_cells, seed=150, stop_mult=0, full_output=full_output, verbose=False)
plt.figure('Rebound')
cell_count = [np.log2(tc[i]['CellCount']/tc[i]['CellCount'][0]) for i in range(len(tc))]
if not envelope:
    for i in range(len(tc)):
        plt.plot(tc[i]['Time']/24., np.log2(tc[i]['CellCount']/tc[i]['CellCount'][0]), 'b', linewidth=2)
else:
    ssa_low = np.percentile(cell_count, 2.5, axis=0)
    ssa_high = np.percentile(cell_count, 97.5, axis=0)
    plt.fill_between(t/24., ssa_low, ssa_high, facecolor='b', alpha=0.25, edgecolor='b')

# ODE simulation
if run_ode:
    for i in range(len(probs)):
        model.parameters['Cell%d_Init' % i].value = 1e4*probs[i]
    x = odesolve(model, t, verbose=False)
    plt.figure('Rebound')
    tc_ode = np.log2(x['Cell_total']/x['Cell_total'][0])
    plt.plot(t/24., tc_ode, 'r', linewidth=5)

plt.annotate("mean: %6.4f" % mean, annot_mean, xycoords='axes fraction', fontsize=18)
plt.annotate("sdev: %6.4f" % sdev, annot_sdev, xycoords='axes fraction', fontsize=18)
plt.annotate("skew: %6.4f" % skew, annot_skew, xycoords='axes fraction', fontsize=18)
#####
plt.annotate(r"DIP$\leq$0: kdiv = %g" % kdiv[0], annot_DIPn, xycoords='axes fraction', fontsize=18)
plt.annotate(r"DIP$>$0: kdiv = %g" % kdiv[-1], annot_DIPp, xycoords='axes fraction', fontsize=18)
#####
# fast_rebounder = 0
# slow_rebounder = 0
# for i in range(1,len(cell_count)):
#     if min(cell_count[i]) < min(cell_count[slow_rebounder]):
#         slow_rebounder = i
#     if cell_count[i][21] > cell_count[fast_rebounder][21]:
#         fast_rebounder = i
# print "slow: %d" % slow_rebounder
# print "fast: %d" % fast_rebounder
# plt.figure('Divergent Responders')
# plt.plot(tc[slow_rebounder]['Time']/24., cell_count[slow_rebounder], 'r', linewidth=3, label='SSA_%d' % slow_rebounder)
# plt.plot(tc[fast_rebounder]['Time']/24., cell_count[fast_rebounder], 'g', linewidth=3, label='SSA_%d' % fast_rebounder)
# plt.ylim(ymin=-3.,ymax=1.)
# plt.legend()
# plt.show()
# quit()
#####
# Clonal time courses
if len(full_output) > 0:
    n = 1
    fig, axs = plt.subplots(1, len(full_output)+1) #, sharey=True)
    fig.canvas.set_window_title('SSA time courses')
for i in full_output:
#     plt.figure('Total cell count v. time')
#     plt.plot(tc[i]['Time']/24., np.log2(tc[i]['CellCount']/tc[i]['CellCount'][0]), linewidth=2, label='SSA_%d' % i)
#     plt.ylim(ymin=-3,ymax=1)
#     plt.xlabel('Time post drug (days)', fontsize=18)
#     plt.ylabel('log2(Normalized cell count)', fontsize=18)
#     plt.xticks(fontsize=18)
#     plt.yticks(fontsize=18)
#     plt.legend(loc=0)
    #####
#     axs[0].set_title('Divergent rebounders')
    axs[0].plot(tc[i]['Time']/24., np.log2(tc[i]['CellCount']/tc[i]['CellCount'][0]), linewidth=2, label='SSA_%d' % i)
    axs[0].set_ylim(ymin=-3,ymax=1)
    axs[0].set_xlabel('Time post drug (days)', fontsize=18)
    axs[0].set_ylabel('log2(Normalized cell count)', fontsize=18)
    axs[0].tick_params(labelsize=18)
    axs[0].legend(loc='upper left')
    #####
#     plt.figure('SSA_%d' % i) #'SSA time courses')
#     axs[n].set_title('SSA_%d' % i)
    axs[n].annotate("SSA_%d" % i, (0.70, 0.95), xycoords='axes fraction', fontsize=18, fontstyle='italic', fontweight='bold')
    # Shink current axis
#     ax = plt.subplot(111)
#     box = ax.get_position()
#     ax.set_position([box.x0, box.y0, box.width*0.9, box.height])
    ###
    print '------'
    print 'SSA_%d' % i
    print '------'
    for j in range(len(DIP)):
        if tc[i]['__s%d' % j][0] > 0.:
            print 's%d:' % j, tc[i]['__s%d' % j][0], '(DIP=%g)' % DIP[j]
            if tc[i]['__s%d' % j][-1] > 0.:
                axs[n].plot(tc[i]['Time']/24., tc[i]['__s%d' % j], linewidth=2, label='s%d' % j)
            elif DIP[j] > 0.:
                axs[n].plot(tc[i]['Time']/24., tc[i]['__s%d' % j], '--', linewidth=3, label='s%d' % j)
            else:
                axs[n].plot(tc[i]['Time']/24., tc[i]['__s%d' % j], '--', linewidth=1, label='s%d' % j)
    axs[n].set_yscale('log', basey=2)
    axs[n].set_ylim(ymax=2**15)
#     plt.legend(loc=0)
    axs[n].legend(loc='upper left') #, bbox_to_anchor=(1, 1), ncol=1)
    axs[n].set_xlabel('Time post drug (days)', fontsize=18)
    axs[n].set_ylabel('Cell count', fontsize=18)
    axs[n].tick_params(labelsize=18)
    n += 1

# Output to file
if OUTFILES.get('envelope'):
    with OUTFILES['envelope'] as o:
        o.write('t,ssa_low,ssa_high')
        if run_ode: o.write(',ode')
        o.write('\n')
        for i in range(len(t)):
            o.write('%g,%g,%g' % (t[i]/24., ssa_low[i], ssa_high[i]))
            if run_ode: o.write('%g' % tc_ode[i])
            o.write('\n')
        o.close()

if OUTFILES.get('timecourses'):
    with OUTFILES['timecourses'] as o:
        o.write('t')
        for i in range(len(tc)):
            o.write(',SSA_%d' % i)
        o.write('\n')
        for i in range(len(t)):
            o.write('%g' % (t[i]/24.))
            for j in range(len(tc)):
                o.write(',%g' % (np.log2(tc[j]['CellCount'][i]/tc[j]['CellCount'][0])))
            o.write('\n')

# fig.savefig(os.path.join(OUTDIR,'SSA_timecourses.pdf'), format='pdf')

plt.show()

# STATS = open(os.path.join(OUTDIR,'stats.txt'), 'a')

# init_cells = []
# for i in range(10):
#     init_cells.append(exp_pops)
# init_cells = np.reshape(init_cells, 10*len(exp_pops))

# plt.savefig(os.path.join(OUTDIR,'tot_cells_v_time_ALL.pdf'), format='pdf')
# plt.savefig(os.path.join(OUTDIR,'tot_cells_v_time_ODE.pdf'), format='pdf')
# plt.savefig(os.path.join(OUTDIR,'cell_types_v_time_ODE_clipped.pdf'), format='pdf')
'''
for tc in timecourses.keys():
    plt.figure("Cell type time courses (SSA_"+str(tc)+")")
    # Shink current axis by 20%
    ax = plt.subplot(111)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width*0.75, box.height])
    ###
    x = timecourses[tc]
    final_pops = [x["__s"+str(i)][-1] for i in range(len(model.species)-2)]
    n_largest = []
    for i in range(5):
        largest_index = numpy.argmax(final_pops)
        n_largest.append(largest_index)
        final_pops[largest_index] = 0.0
    for i in range(len(model.species)):
        if re.match("Cell",str(model.species[i])) and sum(x['__s'+str(i)]) > 0.0:
            if i in n_largest:
                plt.plot(x['time'], x['__s'+str(i)], label="s"+str(i), linewidth=2)
            elif x["__s"+str(i)][-1] == 0:
                plt.plot(x['time'], x['__s'+str(i)], ':', label="s"+str(i), linewidth=1)
            else:
                plt.plot(x['time'], x['__s'+str(i)], label="s"+str(i), linewidth=0.5)
    plt.yscale('log', basey=2)
    plt.xlabel('time')
    plt.ylabel('population')
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1), ncol=3, shadow=True, prop={'size':8})
'''