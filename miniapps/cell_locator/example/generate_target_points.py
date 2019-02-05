import netCDF4
import numpy
import argparse

"""
Generate random target points
"""

parser = argparse.ArgumentParser(description='Generate random target points')
parser.add_argument('-n', metavar='n', type=int, default=11, help='Number of points')
parser.add_argument('-seed', metavar='seed', type=int, default=123, help='Random number seed')
parser.add_argument('-zlo', metavar='zlo', type=float, default=0., help='Low vertical level')
parser.add_argument('-zhi', metavar='zhi', type=float, default=1., help='High vertical level')
parser.add_argument('-f', metavar='f', type=str, default='target_points.nc', help='File name')


args = parser.parse_args()

npts = args.n
numpy.random.seed(args.seed)
lons = 0. + 2*numpy.pi * numpy.random.random((npts,))
lats = -numpy.pi/2. + numpy.pi * numpy.random.random((npts,))
levs = args.zlo + (args.zhi - args.zlo) * numpy.random.random((npts,))
points = numpy.zeros((npts, 3), numpy.float64)
points[:, 0] = lons
points[:, 1] = lats
points[:, 2] = levs

print('Saving the points in file {}'.format(args.f))
nc = netCDF4.Dataset(args.f, 'w')
# create dimensions
three = nc.createDimension('three', 3)
npts = nc.createDimension('npts', npts)
# create variable
target_points = nc.createVariable('target_points', numpy.float64, ('npts', 'three'))
# write
target_points[:] = points

nc.close()

