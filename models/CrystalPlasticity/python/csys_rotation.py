import sys
import numpy as np
import matplotlib.pyplot as plt
import quaternion as quat
import vpython as vp
import time

param = {'kangle': np.pi*1, 'kred': 1.01,
         'energy_exponent': 20}


def main(argv):
    num_points = 10
    num_steps = 500
    draw_steps = 10
    coords = get_initial_directions(num_points)

    vp.canvas(title='Coordinate system rotations', width=1800, height=1200)

    csys = draw_initial(coords)

    energy_hist = []
    energy_hist.append(get_total_energy(coords))
    print('{:4.0f}: {:f} ({:f})'.format(0, energy_hist[-1], param['kangle']))
    for step in range(num_steps):
        coords = update_positions(coords)
        energy_hist.append(get_total_energy(coords))

        print('{:4.0f}: {:f} ({:f})'.format(step, energy_hist[-1], param['kangle']))
        param['kangle'] = param['kangle']/param['kred']
        if (step + 1) % draw_steps == 0:
            # time.sleep(1/(1+step/10))
            time.sleep(0.001)
            update_drawing(csys, coords)

    plt.semilogy(energy_hist)
    print(np.min(energy_hist))
    plt.show()


def get_total_energy(coords):
    energy = 0
    for c0 in coords:
        for c in coords:
            energy = energy + get_coupling_energy(c0, c)
        energy = energy - get_coupling_energy(c0, c0)
    return (energy/(len(coords)**2))**(1.0/param['energy_exponent'])


def get_coupling_energy(coords0, coords1):
    energy = 0
    for v0 in coords0:
        for v1 in coords1:
            energy = energy + 0.5*(np.dot(v0, v1))**param['energy_exponent']

    return energy


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
    arrow = vp.arrow(axis=vp.vec(vec[0], vec[1], vec[2]), up=vp.vec(up[0], up[1], up[2]),
                     color=vp.vec(color[0], color[1], color[2]), make_trail=True)

    return arrow


def update_drawing(csys, coords):
    for cs, coord in zip(csys, coords):
        ctmp = [coord[1], coord[2], coord[0]]
        for a, v, u in zip(cs, coord, ctmp):
            update_arrow(a, v, u)


def update_arrow(a, vec, up):
    a.axis = vp.vec(vec[0], vec[1], vec[2])
    a.up = vp.vec(up[0], up[1], up[2])


# Initial setup
def get_initial_directions(num_points):
    alpha = 2*np.pi*np.random.rand(num_points)
    beta = np.arcsin(2*np.random.rand(num_points)-1)
    gamma = np.pi*np.random.rand(num_points)

    alpha[0] = 0
    beta[0] = 0
    gamma[0] = 0
    coords = []
    for a, b, g in zip(alpha, beta, gamma):
        u0 = np.array([+np.cos(a), np.sin(a), 0])
        v0 = np.array([-np.sin(a), np.cos(a), 0])
        ex = rotate_vector(u0, v0, b)
        ey = rotate_vector(v0, ex, g)
        ez = np.cross(ex, ey)
        if len(coords) == 1:
            x1 = 1
            x2 = (np.sqrt(1/x1**2 + 2)-1/x1)/2
            x3 = -(1/x1 + x2)
            ex = np.array([x1, x1, 1])
            ey = np.array([x2, x3, 1])
            ez = np.array([x3, x2, 1])
            ex = ex / np.linalg.norm(ex)
            ey = ey/np.linalg.norm(ey)
            ez = ez/np.linalg.norm(ez)

        coords.append([ex, ey, ez])

    return coords


# Modification of state
def rotate_vector(vector, axis, angle):
    return quat.rotate_vectors(quat.from_rotation_vector(axis * angle), vector)


def update_positions(coords):

    E0 = get_total_energy(coords)
    sfac = 1
    eincrease = True
    cnt = 0
    while eincrease and cnt<1:
        # Calculate torque vectors for each coordinate set
        coords_new = []
        coords_new.append(coords[0])
        for c0 in coords[1:]:
            t0 = np.zeros(3)
            for c1 in coords:
                t0 = t0 + get_torque(c0, c1)

            coords_new.append(rotate_coords(c0, sfac*t0/len(coords)))
        sfac = sfac/2
        Enew = get_total_energy(coords_new)
        eincrease = Enew > E0
        cnt = cnt + 1

    return coords_new


def get_torque(c0, c1):
    torque = np.zeros(3)
    for v0 in c0:
        for v1 in c1:
            torque = torque - get_normal_direction(v0, v1) * np.dot(v0, v1)

    return torque

def get_normal_direction(vec1, vec2):
    direction = np.cross(vec1, vec2)
    norm = np.linalg.norm(direction)

    if norm < 1e-12:
        direction = np.cross(np.random.rand(3), vec1)
        norm = np.linalg.norm(direction)

    return direction/(norm+1e-12)


def rotate_coords(coords, torque):
    angle = np.linalg.norm(torque)*param['kangle']
    axis = torque/np.linalg.norm(torque)
    new_coords = []
    for v in coords:
        new_coords.append(rotate_vector(v, axis, angle))

    return new_coords




if __name__=='__main__':
    main(sys.argv)