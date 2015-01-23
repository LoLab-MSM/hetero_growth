from pysb import *
import pysb.core
from pysb.macros import *
import numpy as np

model = Model('PolyGrowth')

def define_model(kdiv, kdeg):

    # Cell[i] -> Cell[i] + Cell[i]  k_div[i]
    # Cell[i] -> 0  k_deg[i]
    
    if len(kdiv) != len(kdeg):
        raise Exception('kdiv and kdeg must have equal lengths.')
    
    Monomer('Cell', ['type'], {'type': [str(i) for i in range(len(kdiv))] + ['Q']})
    [Initial(Cell(type=str(i)), Parameter('Cell%d_Init' % i)) for i in range(len(kdiv))]
    [Rule('divide_Cell%d' % i, Cell(type=str(i)) >> Cell(type=str(i)) + Cell(type=str(i)), Parameter('kdiv_%d' % i, kdiv[i])) for i in range(len(kdiv))]
    [degrade(Cell(type=str(i)), Parameter('kdeg_%d' % i, kdeg[i])) for i in range(len(kdiv))]
    Observable("Cell_total", Cell())

#     Parameter('kqp', 0.05)
#     Parameter('kqm', 0.05)
#     [Rule('quiesce_Cell%d' %i, Cell(type=str(i)) <> Cell(type='Q'), kqp, kqm) for i in range(len(kdiv))]

#     model.initial_conditions = np.zeros(n_cell_types, dtype=(pysb.core.ComplexPattern, Parameter)).tolist()
#     model.initial_conditions[counter] = (pysb.core.as_complex_pattern(Cell(type=str(counter))), Parameter('Init_%d' % counter, init_pops[counter]))
