import hetero_growth as hg
from hetero_growth import model
from pysb.integrate import odesolve
from pysb.bng import run_ssa
from pysb.bng import set_bng_path
import matplotlib.pyplot as plt
import numpy as np
import snorm
import os

set_bng_path('/Users/lopezlab/git/bionetgen/bng2')
OUTDIR = '/Users/lopezlab/Documents/workspace/HeteroGrowth/GraphicalAbstract'

n_sims = 20
envelope = False
run_ode = True
output_to_file = True

total_cells = [3000]*n_sims
mean = -0.00581
sdev = 0.00571
skew = -0.566
kdiv_base = 0.05 
nbins = 40

MEAN = [mean, 2.*mean, 0.5*mean]
SDEV = [sdev, 2.*sdev, 0.5*sdev]
SKEW = [skew]*3

colors = ['b','g','r']

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

DIP = np.linspace(-0.05,0.05,nbins+1)
kdiv = np.array([kdiv_base]*len(DIP))
# kdiv = np.array([kdiv_base for d in DIP if d <= 0.] + [10.*kdiv_base for d in DIP if d > 0.])
kdeg = np.array([kdiv[i]-DIP[i] for i in range(len(DIP)) if kdiv[i]-DIP[i] > 0.])
kdiv = kdiv[:len(kdeg)]
DIP = DIP[:len(kdeg)]
hg.define_model(kdiv, kdeg)
    
for nn in range(len(MEAN)):
    mean = MEAN[nn]
    sdev = SDEV[nn]
    skew = SKEW[nn]
    plt.figure('DIP rates')
    plt.plot(np.linspace(-0.05,0.05,3*40), snorm.pdf(np.linspace(-0.05,0.05,3*40),mean,sdev,skew), linewidth=5) #, color='g')
#     plt.annotate("mean: %.3g" % mean, (0.65, 0.90), xycoords='axes fraction', fontsize=18, color='g')
#     plt.annotate("sdev: %.3g" % sdev, (0.65, 0.82), xycoords='axes fraction', fontsize=18, color='g')
#     plt.annotate("skew: %.3g" % skew, (0.65, 0.74), xycoords='axes fraction', fontsize=18, color='g')
#     plt.show()
#     quit()
    
    # probs = ss.norm.pdf(DIP,mean,sdev)*(DIP[1]-DIP[0])
    probs = snorm.pdf(DIP,mean,sdev,skew)*(DIP[1]-DIP[0])
    probs /= sum(probs)
    init_cells = []
    np.random.seed(1)
    init_cells.append([np.random.multinomial(n, probs, size=1)[0] for n in total_cells])
    init_cells = np.reshape(init_cells, (len(total_cells),len(DIP)))
    
    n_days = 50
    t = np.linspace(0, n_days*24, n_days+1)
    
    # SSA simulations
    tc = run_ssa_sims(n_sims, t, init_cells, seed=150, stop_mult=32, verbose=False)
    plt.figure('Rebound')
    cell_count = [np.log2(tc[i]['CellCount']/tc[i]['CellCount'][0]) for i in range(len(tc))]
    if envelope:
        ssa_low = np.percentile(cell_count, 2.5, axis=0)
        ssa_high = np.percentile(cell_count, 97.5, axis=0)
        plt.fill_between(t/24., ssa_low, ssa_high, facecolor='b', alpha=0.25, edgecolor='b')
    else:
        for i in range(len(tc)):
            plt.plot(tc[i]['Time']/24., np.log2(tc[i]['CellCount']/tc[i]['CellCount'][0]), lw=1, ls='--', color=colors[nn])
    
    # ODE simulation
    if run_ode:
        for i in range(len(probs)):
            model.parameters['Cell%d_Init' % i].value = 1e4*probs[i]
        x = odesolve(model, t, verbose=False)
        plt.figure('Rebound')
        plt.ylim(ymin=-2,ymax=5)
        tc_ode = np.log2(x['Cell_total']/x['Cell_total'][0])
        plt.plot(t/24., tc_ode, linewidth=5) #, color='r')

    # Output to file
    if output_to_file:
        # distribution
#         np.linspace(-0.05,0.05,3*40), snorm.pdf(np.linspace(-0.05,0.05,3*40),mean,sdev,skew)
#         output_dist = open(os.path.join(OUTDIR, 'distribution_%d.csv' % nn), 'w')
        x = np.linspace(-0.05,0.05,3*40+1)
        x = np.reshape(x,(len(x),1))
        y = np.reshape(snorm.pdf(x,mean,sdev,skew),(len(x),1))
        np.savetxt(os.path.join(OUTDIR,'distribution_%d.csv' % nn), np.hstack((x,y)), delimiter=',', header='dip,density', comments='')
#         output_dist.close()
        # time courses
        output_tc = open(os.path.join(OUTDIR, 'timecourses_%d.csv' % nn), 'w')
        if envelope:
            output_tc.write('t,ssa_low,ssa_high')
            if run_ode: output_tc.write(',ode')
            output_tc.write('\n')
            for i in range(len(t)):
                output_tc.write('%g,%g,%g' % (t[i]/24., ssa_low[i], ssa_high[i]))
                if run_ode: output_tc.write(',%g' % tc_ode[i])
                output_tc.write('\n')
        else:
            output_tc.write('t')
            for i in range(len(tc)):
                output_tc.write(',SSA_%d' % i)
            if run_ode: output_tc.write(',ode')
            output_tc.write('\n')
            for i in range(len(t)):
                output_tc.write('%g' % (t[i]/24.))
                for j in range(len(tc)):
                    if i < len(tc[j]['CellCount']):
                        output_tc.write(',%g' % (np.log2(tc[j]['CellCount'][i]/tc[j]['CellCount'][0])))
                    else:
                        output_tc.write(',\t')
                if run_ode: output_tc.write(',%g' % tc_ode[i])
                output_tc.write('\n')
        output_tc.close()

# fig.savefig(os.path.join(OUTDIR,'SSA_timecourses.pdf'), format='pdf')
plt.show()
