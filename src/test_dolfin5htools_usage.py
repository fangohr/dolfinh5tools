def test_field2_has_more_timesteps_than_field1():
    pass


def test_reading():

    h5file = h5tools.open(filename)
    for field in h5file.fields:
        for t in h5file.times(field):
            _ = h5file.read(field, t)


def test_reading_fixed_field():

    h5file = h5tools.open(filename, field='m')
    for t in h5file.times:
        _ = h5file.read(t)

def test_writing():

    h5file = h5tools.open(filename, 'w')
    for t in simulation_time_steps:
        # compute field1
        # compute field2
        for field in [field1, field2]:
            h5file.write(field, time, fieldname)




















