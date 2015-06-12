import json
from collections import OrderedDict

import numpy as np
import dolfin as df

from savingdata import Create, Read

mesh = df.UnitSquareMesh(10, 10)
filename = 'file_vector'
functionspace = df.VectorFunctionSpace(mesh, 'CG', 1, 3)
f = df.Function(functionspace)
t_array = np.linspace(0, 1e-9, 5)

# Load data.
ld = Read(filename)
for t in t_array:
    f.assign(df.Constant((1+t, 2, 3)))
    f_loaded = ld.load_field('f', t)

    assert df.assemble(f[0]*df.dx) == df.assemble(f_loaded[0]*df.dx)
ld.close()
