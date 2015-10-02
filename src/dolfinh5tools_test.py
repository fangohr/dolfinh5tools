import os
import json
import subprocess
from collections import OrderedDict

import numpy as np
import dolfin as df

from dolfinh5tools import openh5

mesh = df.UnitSquareMesh(10, 10)
t_array = np.linspace(0, 1e-9, 5)


def test_writing_usage():
    filename = 'file_usage'
    functionspace = df.VectorFunctionSpace(mesh, 'CG', 1, 3)
    f1 = df.Function(functionspace)
    f2 = df.Function(functionspace)

    h5file = openh5(filename, functionspace, mode='w')

    h5file.save_mesh()

    for t in t_array:
        f1.assign(df.Constant((t, 0, 0)))
        f2.assign(df.Constant((t, 1, 0)))
        
        h5file.write(f1, 'f1', t)
        h5file.write(f2, 'f2', t)


def test_field2_has_more_timesteps_than_field1():
    filename = 'file_usage_diff'
    functionspace = df.VectorFunctionSpace(mesh, 'CG', 1, 3)
    f1 = df.Function(functionspace)
    f2 = df.Function(functionspace)

    h5file = openh5(filename, functionspace, mode='w')

    h5file.save_mesh()

    for t in t_array:
        f1.assign(df.Constant((t, 0, 0)))
        f2.assign(df.Constant((t, 1, 0)))
        
        h5file.write(f1, 'f1', t)
        h5file.write(f2, 'f2', t)
    f2.assign(df.Constant((50, 1, 0)))
    h5file.write(f2, 'f2', t)

    h5file = openh5(filename, mode='r')
    print h5file.fields()
    for field in h5file.fields():
        for t in h5file.times(field):
            h5file.read(t_array[1], field)

    os.system('rm file_usage_diff.h5')
    os.system('rm file_usage_diff.json')

def test_reading_usage():
    filename = 'file_usage'
    h5file = openh5(filename, mode='r')
    for field in h5file.fields():
        for t in h5file.times(field):
            h5file.read(t_array[1], field)


def test_reading_fixed_field_usage():
    filename = 'file_usage'
    h5file = openh5(filename, field_name='f1', mode='r')
    for t in h5file.times():
        _ = h5file.read(t)

    os.system('rm file_usage.h5')
    os.system('rm file_usage.json')

def test_save_vector_field_data():
    filename = 'file_vector'
    functionspace = df.VectorFunctionSpace(mesh, 'CG', 1, 3)
    f = df.Function(functionspace)

    sd = openh5(filename, functionspace, mode='w')

    sd.save_mesh()

    for i in range(len(t_array)):
        f.assign(df.Constant((t_array[i], 0, 0)))
    
        sd.write(f, 'f', t_array[i])

    sd.close()


def test_save_scalar_field_data():
    filename = 'file_scalar'
    functionspace = df.FunctionSpace(mesh, 'CG', 1)
    f = df.Function(functionspace)

    sd = openh5(filename, functionspace, mode='w')

    sd.save_mesh()

    for i in range(len(t_array)):
        f.assign(df.Constant(t_array[i]))
    
        sd.write(f, 'f', t_array[i])

    sd.close()


def test_load_vector_data():
    filename = 'file_vector'
    functionspace = df.VectorFunctionSpace(mesh, 'CG', 1, 3)
    f = df.Function(functionspace)

    ld = openh5(filename, mode='r')

    for t in t_array:
        f_loaded = ld.read(t, 'f')

        f.assign(df.Constant((t, 0, 0)))
        assert np.all(f.vector().array() == f_loaded.vector().array())

    ld.close()

    os.system('rm file_vector.h5')
    os.system('rm file_vector.json')


def test_load_scalar_data():
    filename = 'file_scalar'
    functionspace = df.FunctionSpace(mesh, 'CG', 1)
    f = df.Function(functionspace)

    ld = openh5(filename, mode='r')

    for t in t_array:
        f_loaded = ld.read(t, 'f')

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


