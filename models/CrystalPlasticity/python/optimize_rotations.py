import sys
import numpy as np
import scipy.optimize as opt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import quaternion as quat

import time

global_param = {'energy_exponent': 20, 'eps_num_der': 1e-8}
plot_result = False
write_result = True

v9order = np.array([[1,1], [2,2], [3,3], [1,2], [2,3], [3,1], [1,3], [2,1], [3,2]])-1
m3x3order= np.array([[1,4,6], [8,2,7], [6,9,3]])-1


if plot_result:
    import vpython as vp


# Optimization of csys rotations, specified by a vector whose length is the magnitude of rotation about it's own axis
def main(argv):
    num_points = np.arange(2, 25)
    num_tests = 10
    Emin = []
    Emean= []
    Estd = []
    Eall = []
    for num in num_points:
        Etmp = []
        for t in range(num_tests):
            Etmp.append(get_minimum_energy(num, t))
        Emin.append(np.min(Etmp))
        Emean.append(np.average(Etmp))
        Estd.append(np.std(Etmp))
        Eall.append(Etmp)
        print('{:3.0f} {:.10f} {:.10f} {:.10f}'.format(num, Emin[-1], Emean[-1], Estd[-1]))

    plt.figure()
    plt.errorbar(num_points, Emean, Estd)

    plt.figure()
    plt.semilogx(num_points, Emin, c='r')

    plt.figure()
    plt.plot(num_points, Emin, c='r', label='Minimum')
    plt.errorbar(num_points, Emean, Estd, linestyle='--', label='Mean and std. deviation')
    plt.legend()

    np.savetxt('all_errors.txt', np.transpose(np.array(Eall)))

    plt.show()


def get_minimum_energy(num_points, filenum):
    num_optim = num_points - 1
    r0 = get_initial_guess(num_optim)

    # Transform into variable array
    x0 = np.reshape(r0, (num_optim * 3))

    res = opt.minimize(objective_function, x0, method='Powell',
                       options={'xtol': 1e-5, 'disp': False, 'maxfev': 10000})

    if write_result:
        write_result_file(num_points, res, filenum)


    if plot_result:
        rot_vectors = np.reshape(res.x, (3, num_optim))
        coords = get_coords(rot_vectors)
        vp.canvas(title='Coordinate system rotations', width=1800, height=1200)
        draw_initial(coords)

    return res.fun


# Result write routines
def write_result_file(num_points, res, filenum):
    fname = 'N'+str(num_points)+'_csys_rot_' + str(filenum) + '.txt'
    rot_vectors = np.reshape(res.x, (3, num_points-1))
    with open(fname, 'w') as fid:
        fid.write('! Results for N = ' + str(num_points) + ' coordinate systems\n')
        fid.write('! Energy measure = ' + str(res.fun) +
                  '(energy exponent = ' + str(global_param['energy_exponent']) + ')\n')
        fid.write('! Rotation matrices:\n')
        fid.write('{:<12s}'.format('!'))
        for ind in v9order:
            fid.write('{:>19s}{:1.0f}{:1.0f}'.format('Q', ind[0]+1, ind[1]+1))
        fid.write('\n')
        rownum = 1
        for rv in np.transpose(rot_vectors):
            print_rot_matrix(fid, rv, rownum)
            rownum = rownum + 1

def print_rot_matrix(fid, rot_vector, rownum):
    rot_matrix = quat.as_rotation_matrix(quat.from_rotation_vector(rot_vector))
    rot_mat_v9 = [rot_matrix[ind[0],ind[1]] for ind in v9order]
    fid.write('Q(:, {:0=3.0f}) = '.format(rownum))
    fid.write(('(/' + '{:19.16f}, '*8 + '{:19.16f}/)' + '\n').format(*rot_mat_v9))

# Plotting routines
def draw_initial(coords):
    color = [1,0,0]
    csys = []
    for cor in coords:
        csys.append(plot_csys(cor[0], cor[1], cor[2], color))
        color = [color[i] for i in [1, 2, 0]]

    return csys


def plot_csys(ex, ey, ez, color):
    csys = []
    csys.append(plot_arrow(ex, ey, color))
    csys.append(plot_arrow(ey, ez, color))
    csys.append(plot_arrow(ez, ex, color))
    return csys


def plot_arrow(vec, up, color):
    vp.arrow(axis=vp.vec(-vec[0], -vec[1], -vec[2]), up=vp.vec(-up[0], -up[1], -up[2]),
             color=vp.vec(color[0], color[1], color[2]))
    arrow = vp.arrow(axis=vp.vec(vec[0], vec[1], vec[2]), up=vp.vec(up[0], up[1], up[2]),
                     color=vp.vec(color[0], color[1], color[2]))

    return arrow


def objective_function(x):
    rot_vectors = np.reshape(x, (3, int(x.size/3)))
    coords = get_coords(rot_vectors)

    return get_total_energy(coords)


def obj_num_deriv(x):
    e0 = objective_function(x)

    de_dx = np.zeros((x.size))
    for i in range(x.size):
        xtmp = x
        xtmp[i] = xtmp[i] + global_param['eps_num_der']
        de_dx[i] = (objective_function(xtmp)-e0)/global_param['eps_num_der']

    return de_dx


def get_coords(rot_vectors):
    coords = [[np.array([1, 0, 0]), np.array([0, 1, 0]), np.array([0, 0, 1])]]
    for rv in np.transpose(rot_vectors):
        coords.append(rotate_coord(coords[0], rv))

    return coords


def rotate_coord(coord0, rotation_vector):
    coord = []
    for c0 in coord0:
        coord.append(quat.rotate_vectors(quat.from_rotation_vector(rotation_vector), c0))

    return coord


def get_initial_guess(num_directions):
    alpha = np.random.rand(num_directions)*np.pi            # Initial rotation about z axis
    beta = np.arccos(2*np.random.rand(num_directions)-1)  # Rotation vector angle to z-axis
    gamma = np.random.rand(num_directions)*np.pi*2          # Rotation about rotation vector

    x0 = np.array([np.sin(beta)*np.cos(alpha), np.sin(beta)*np.sin(alpha), np.cos(beta)]) @ np.diag(gamma)

    return x0


def get_total_energy(coords):
    energy = 0
    for c0 in coords:
        for c in coords:
            energy = energy + get_coupling_energy(c0, c)
        energy = energy - get_coupling_energy(c0, c0)
    return (energy/(len(coords)**2))**(1.0/global_param['energy_exponent'])


def get_coupling_energy(coords0, coords1):
    energy = 0
    for v0 in coords0:
        for v1 in coords1:
            energy = energy + 0.5*(np.dot(v0, v1))**global_param['energy_exponent']

    return energy


def set_axes_equal(ax):
    '''Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc..  This is one possible solution to Matplotlib's
    ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

    Input
      ax: a matplotlib axis, e.g., as output from plt.gca().
    '''

    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5*max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])


if __name__=='__main__':
    main(sys.argv)