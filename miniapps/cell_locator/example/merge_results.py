import netCDF4
import numpy
import argparse
import glob

"""
Merge per processor interpolation results
"""

parser = argparse.ArgumentParser(description='Merge results stored in different NetCDF files')
parser.add_argument('-i', metavar='i', type=str, default='cell_ids_[0-9]*.nc', 
	                help='Result file names (cell_ids_[0-9]*.nc)')
parser.add_argument('-o', metavar='o', type=str, default='cell_ids_merged.nc', 
	                help='Output file name after merging (cell_ids_merged.nc)')

args = parser.parse_args()

npts = 0
point_ids_0 = []
cell_ids_0 = []
dist_error_square = []
pcoords = []
for filename in glob.glob(args.i):

    print('reading {}...'.format(filename))
    nc = netCDF4.Dataset(filename)

    local_point_ids_0 = nc.variables['point_ids_0'][:]
    local_cell_ids_0 = nc.variables['cell_ids_0'][:]
    local_dist_error_square = nc.variables['dist_error_square'][:]
    local_pcoords = nc.variables['pcoords'][:]
    nc.close()

    npts += len(local_point_ids_0)

    n = len(local_point_ids_0)
    for i in range(n):
        point_ids_0.append(local_point_ids_0[i])
        cell_ids_0.append(local_cell_ids_0[i])
        dist_error_square.append(local_dist_error_square[i])
        pcoords.append(local_pcoords[i, :])


# write the data
print('writing the combined results to file {}'.format(args.o))
nc = netCDF4.Dataset(args.o, 'w')
# create dimensions
three = nc.createDimension('three', 3)
npts = nc.createDimension('npts', npts)
# create variables
point_ids_0_var = nc.createVariable('point_ids_0', numpy.int, ('npts',))
cell_ids_0_var = nc.createVariable('cell_ids_0', numpy.int64, ('npts',))
dist_error_square_var = nc.createVariable('dist_error_square', numpy.float64, ('npts',))
pcoords_var = nc.createVariable('pcoords', numpy.float64, ('npts', 'three'))
# write
point_ids_0_var[:] = point_ids_0
cell_ids_0_var[:] = cell_ids_0
dist_error_square_var[:] = dist_error_square
pcoords_var[:] = pcoords


