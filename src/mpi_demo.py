import numpy as np
import dolfin as df

from savingdata import Create, Read

mesh = df.UnitSquareMesh(10, 10)
filename = 'file_vector'
functionspace = df.VectorFunctionSpace(mesh, 'CG', 1, 3)
f = df.Function(functionspace)
t_array = np.linspace(0, 1e-9, 5)

# Save data.
sd = Create(filename, functionspace)
sd.save_mesh()
for i in range(len(t_array)):
    f.assign(df.Constant((1, 2, 3)))
    sd.save_field(f, 'f', t_array[i])
sd.close()

# Load data.
ld = Read(filename)
for t in t_array:
    f_loaded = ld.load_field('f', t)
ld.close()
