import numpy as np
import dolfin as df

from dolfinh5tools import Create, Read

mesh = df.UnitSquareMesh(10, 10)
filename = 'file_mpi'
functionspace = df.VectorFunctionSpace(mesh, 'CG', 1, 3)
f = df.Function(functionspace)
t_array = np.linspace(0, 1, 5)

# Save data.
sd = Create(filename, functionspace)
sd.save_mesh()
for t in t_array:
    f.assign(df.Constant((1 + t, 2, 3)))
    sd.save_field(f, 'f', t)
sd.close()
