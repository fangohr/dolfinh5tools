import os
import json
import subprocess
from collections import OrderedDict

import numpy as np
import dolfin as df

from savingdata import Create, Read

mesh = df.UnitSquareMesh(10, 10)
t_array = np.linspace(0, 1e-9, 5)

def test_save_vector_field_data():
    filename = 'file_vector'
    functionspace = df.VectorFunctionSpace(mesh, 'CG', 1, 3)
    f = df.Function(functionspace)

    sd = Create(filename, functionspace)

    sd.save_mesh()

    for i in range(len(t_array)):
        f.assign(df.Constant((t_array[i], 0, 0)))
    
        sd.save_field(f, 'f', t_array[i])

    sd.close()

def test_save_scalar_field_data():
    filename = 'file_scalar'
    functionspace = df.FunctionSpace(mesh, 'CG', 1)
    f = df.Function(functionspace)

    sd = Create(filename, functionspace)

    sd.save_mesh()

    for i in range(len(t_array)):
        f.assign(df.Constant(t_array[i]))
    
        sd.save_field(f, 'f', t_array[i])

    sd.close()


def test_load_vector_data():
    filename = 'file_vector'
    functionspace = df.VectorFunctionSpace(mesh, 'CG', 1, 3)
    f = df.Function(functionspace)

    ld = Read(filename)

    for t in t_array:
        f_loaded = ld.load_field('f', t)

        f.assign(df.Constant((t, 0, 0)))
        assert np.all(f.vector().array() == f_loaded.vector().array())

    ld.close()

    os.system('rm file_vector.h5')
    os.system('rm file_vector.json')


def test_load_scalar_data():
    filename = 'file_scalar'
    functionspace = df.FunctionSpace(mesh, 'CG', 1)
    f = df.Function(functionspace)

    ld = Read(filename)

    for t in t_array:
        f_loaded = ld.load_field('f', t)

        f.assign(df.Constant(t))
        assert np.all(f.vector().array() == f_loaded.vector().array())

    ld.close()

    os.system('rm file_scalar.h5')


def test_saved():
    jsonfilename = 'file_scalar.json'
    with open(jsonfilename) as jsonfile:
        jsonData = json.load(jsonfile, object_pairs_hook=OrderedDict)
    jsonfile.close()

    # assert times from keys
    for i, t in enumerate(t_array):
        name = "f{}".format(i)
        jsonTime = jsonData['f']['data'][name]
        assert(jsonTime == t)

    # assert times by iterating through values (as they should be ordered)
    index = 0
    for jsonTime in jsonData['f']['data'].itervalues():
         assert(jsonTime == t_array[index])
         index += 1
    os.system('rm file_scalar.json')

def test_mpi():
    for i in range(1, 4):
        for j in range(1, 4):
            print 'Writing with ' + str(i) + ' and reading with ' + \
                str(j) + ' cores.'
            os.system('mpirun -np ' + str(i) + ' python mpi_demo.py')
            subprocess.check_call('mpirun -np ' + str(j) + ' python mpi_check.py',
                                  shell=True)
    os.system('rm file_mpi.h5')
    os.system('rm file_mpi.json')
