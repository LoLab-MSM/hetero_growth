import numpy as np
import scipy.stats as ss
import matplotlib.pyplot as plt
import matplotlib.colors
import matplotlib.cm

# colors = matplotlib.colors.cnames.keys()
# colors = ['blue', 'green', 'red', 'cyan', 'magenta', 'DarkOrange', 'indigo', 'maroon', 'black']
# colors = ['0.0', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8']
colors = [matplotlib.cm.spectral(i) for i in np.linspace(0.9, 0.2, 9)]
index = 0

data = np.genfromtxt('/Users/lopezlab/git/LungCancer/AvgDrugTitration/Jameson/DS8-DS20_fulldataset_Jan2015.csv', 
                         delimiter=',', dtype=None, names=True)

print data.dtype.names

conc = []
for d in data:
    if d['Erlotinib'] not in conc:
        conc.append(d['Erlotinib'])

expts = []
for c in np.sort(conc):
    expts += [(d['Column'], d['Row'], d['plateID'], d['Subline'], d['Erlotinib']) for d in data if d['Timept'] == 1 and d['Erlotinib'] == c]

lcolors = {}
record = {}

for i in range(len(expts)):
    
    print i
    
    col = expts[i][0]
    row = expts[i][1]
    plateID = expts[i][2]
    subline = expts[i][3]
    erl = expts[i][4] if expts[i][4] != np.min(conc) else 0.0
    
    plt.figure(subline)
    
    # Get time points
    t = [d['Time'] for d in data
        if  d['Column'] == col
        and d['Row'] == row
        and d['plateID'] == plateID][1:]
    
    # Get cell counts
    c = [float(d['CellNucleus']) for d in data
         if  d['Column'] == col
         and d['Row'] == row
         and d['plateID'] == plateID][1:]

    if erl not in lcolors.keys():
        lcolors[erl] = colors[index]
        index += 1
    
    if subline not in record.keys():
        record[subline] = [erl]
        plt.plot(t, np.log2(np.array(c)/c[0]), lw=2, color=lcolors[erl], label=str(erl))
        plt.plot([0]+t+[plt.xticks()[0][-1]], [0]*(len(t)+2), '--', lw=4, color='black')
    else:
        if erl not in record[subline]:
            record[subline] += [erl]
            plt.plot(t, np.log2(np.array(c)/c[0]), lw=2, color=lcolors[erl], label=str(erl))
        else:
            plt.plot(t, np.log2(np.array(c)/c[0]), lw=2, color=lcolors[erl])
            
    plt.legend(loc=0)
    plt.ylabel(r'log$_2$(Cell count)')
    plt.xlabel('Time (h)')
    plt.ylim(ymin=-3, ymax=7)
    plt.yticks(range(-3,8))
    plt.title(subline)

plt.show()