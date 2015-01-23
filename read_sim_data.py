import numpy as np
import matplotlib.pyplot as plt
import re

kpts = np.linspace(0.01, 0.04, 151) # kdiv
dpts = np.linspace(-0.011, -0.013, 21) # DIP
data = np.zeros((len(kpts),len(dpts)))

f = open('/Users/lopezlab/Documents/workspace/HeteroGrowth/stats.txt', 'r')
kdiv = []
DIP  = []
pval = []
for line in f:
    x = re.match(r'kdiv:\s+(.+)', line)
    y = re.match(r'DIP:\s+(.+)', line)
    z = re.match(r'pval:\s+(.+)', line)
    if   x: kdiv.append(float(x.group(1)))
    elif y: DIP.append(float(y.group(1)))
    elif z: pval.append(float(z.group(1)))

count = 0
for i in range(len(kpts)):
    for k in range(len(kdiv)):
        if int(round(kpts[i]*1e4)) == int(round(kdiv[k]*1e4)):
            for j in range(len(dpts)):
                if int(round(dpts[j]*1e4)) == int(round(DIP[k]*1e4)):
                    if pval[k] > 0.1: 
                        data[i][j] = pval[k]
                        count += 1 
                    else: 
                        data[i][j] = 0.0

max_pval = max(max(data[i]) for i in range(len(data)))
max_index = 0
for i in range(len(pval)):
    if pval[i] == max_pval:
        max_index = i
print max_pval
print max_index
print "kdiv:", kdiv[max_index]
print "DIP:", DIP[max_index]
print "pval:", pval[max_index]

plt.pcolormesh(dpts, kpts, data, vmin=0, vmax=1) 
plt.axis([dpts.min(), dpts.max(), kpts.min(), kpts.max()])
plt.colorbar(label = 'p-value')
plt.xlabel('DIP rate (kdiv-kdeg)')
plt.ylabel('kdiv (1/h)')

plt.show()
