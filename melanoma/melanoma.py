from pysb import *
from pysb.macros import *

# Cell[i] -> Cell[i] + Cell[i]  k_div[i]
# Cell[i] -> 0  k_deg[i]
# dCell[i]/dt = (k_div[i]-k_deg[i])*Cell[i]

model = Model('Melanoma')

Na = 6.022e23 # molec/mol
vol = 2e-4 # L
DIP_rate = 0.005 #-0.013305153 # k_div - k_deg

Parameter('Cell_0', 100)
Parameter('drug_0', 2e-6*Na*vol)
Parameter('k_div', 0.05)
Parameter('k_deg', (k_div.value-DIP_rate)/drug_0.value)

Monomer('Cell')
Monomer('drug')

Initial(Cell(), Cell_0)
Initial(drug(), drug_0)

Rule('Deg_drug', Cell() + drug() >> drug(), k_deg)
Rule('Division', Cell() >> Cell() + Cell(), k_div)

Observable("Cell_total", Cell())

