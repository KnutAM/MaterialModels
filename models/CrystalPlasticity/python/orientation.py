import sys
import numpy as np
import matplotlib.pyplot as plt
import quaternion as quat
from mpl_toolkits.mplot3d import Axes3D
import time
import vpython as vp


def main(argv):
    num_points = 4
    num_steps = 10
    plotstep = 2000
    param = {'k': 1, 'max_rot': np.pi/8, 'kred': 0e-3, 'moment_decrease': 1.5}
    u, v = get_initial_directions(num_points)

    if plotstep < num_steps:
        fig1 = plt.figure(1)
        ax1 = fig1.add_subplot(111, projection='3d')
        ax1.set_xlim([0, 1])
        ax1.set_ylim([0, 1])
        ax1.set_zlim([0, 1])
        cstr = 'r'
        fig0 = plt.figure(0)
        ax0 = fig0.add_subplot(111)


    energy_history = []
    for step in range(num_steps):
        u, v = calculate_positions(u, v, param)
        energy = get_energy(u, v)
        energy_history.append(energy)
        print('Step ' + str(step) + ' of ' + str(num_steps))
        param['k'] = param['k']/(1+param['kred'])
        if ((step+1) % plotstep) == 0:
            ax0.plot(energy_history)
            plt.pause(0.5)
            plot_step(ax1, u, v, cstr)
            fig1.show()
            plt.pause(0.5)
            cstr = 'b'


    plot_boxes(u, v)
    print_rotation_matrices(u, v)
    print(energy_history[0])
    print(energy_history[-1])
    print(param['k'])
    if plotstep > num_steps:
        plt.plot(energy_history)
        plt.show(block=False)
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111, projection='3d')
        ax1.set_xlim([0, 1])
        ax1.set_ylim([0, 1])
        ax1.set_zlim([0, 1])

    plot_step(ax1, u, v, 'g')
    fig1.show()
    plt.show()


def plot_boxes(u, v):
    vp.canvas()
    for uu, vv in zip(u, v):
        vp.box(pos=vp.vec(0, 0, 0), axis=vp.vec(uu[0], uu[1], uu[2]), length=1, width=1, height=1,
               up=vp.vec(vv[0], vv[1], vv[2]),
               color=vp.vec(np.random.rand(1)[0], np.random.rand(1)[0], np.random.rand(1)[0]))

def plot_step(ax, u, v, cstr):

    for uu, vv in zip(u, v):
        soa = np.array([np.append(np.zeros(3), uu), np.append(np.zeros(3), vv)])
        X, Y, Z, U, V, W = zip(*soa)
        ax.quiver(X, Y, Z, U, V, W, color=cstr)


def print_rotation_matrices(u, v):
    nr = 1
    for uu, vv in zip(u, v):
        ww = np.cross(uu, vv)
        row = 1
        for vec in [uu, vv, ww]:
            print('Q(:, ' + str(row) + ', ' + str(nr) + ') = (/', end='')
            print('{:+.15f}'.format(vec[0]), end='')
            for ind in vec[1:]:
                print(', {:+.15f}'.format(ind), end='')
            print('/)')
            row = row + 1

        print('')
        nr = nr + 1


def print_euler_angles(u, v):
    for uu, vv in zip(u, v):
        eulang = get_euler_angle(u, v)


def get_euler_angle(u, v):
        w = np.cross(u, v)
        R = np.stack((u, v, w), axis=-1)

        sy = np.linalg.norm(R[0:1, 0])
        if sy > 1e-8:
            x = np.arctan2(R[2, 1], R[2, 2])
            y = np.arctan2(-R[2, 0], sy)
            z = np.arctan2(R[1, 0], R[0, 0])
        else:
            x = np.arctan2(-R[1, 2], R[1, 1])
            y = np.arctan2(-R[2, 0], sy)
            z = 0

        return np.array([x, y, z])


def eulerAnglesToRotationMatrix(theta):
    R_x = np.array([[1, 0, 0],
                    [0, np.cos(theta[0]), -np.sin(theta[0])],
                    [0, np.sin(theta[0]), np.cos(theta[0])]
                    ])

    R_y = np.array([[np.cos(theta[1]), 0, np.sin(theta[1])],
                    [0, 1, 0],
                    [-np.sin(theta[1]), 0, np.cos(theta[1])]
                    ])

    R_z = np.array([[np.cos(theta[2]), -np.sin(theta[2]), 0],
                    [np.sin(theta[2]), np.cos(theta[2]), 0],
                    [0, 0, 1]
                    ])

    R = np.dot(R_z, np.dot(R_y, R_x))

    return R


def get_energy(u, v):
    energy = 0
    for u1, v1 in zip(u, v):
        for u2, v2 in zip(u, v):
            energy = energy + np.dot(u1, u2)**2 + np.dot(u1, v2)**2 + np.dot(v1, u2)**2 + np.dot(v1, v2)**2

    return energy


def calculate_positions(u, v, param):
    moment = get_torque_sum(u, v) + get_torque_sum(v, u)
    eincrease = True
    E0 = get_energy(u, v)
    while eincrease:
        u_new = update_positions(u, moment, param)
        v_new = update_positions(v, moment, param)
        E1 = get_energy(u_new, v_new)
        eincrease = E1 > E0
        new_moment = []
        for m in moment:
            new_moment.append(m/param['moment_decrease'])
        moment = new_moment

    return u_new, v_new


def update_positions(vec, moment, param):
    new_vec = []
    for v, m in zip(vec, moment):
        mag = param['max_rot'] * (1 - np.exp(-np.linalg.norm(m) * param['k']))
        new_vec.append(rotate_vector(v, m, mag))

    return new_vec


def get_torque_sum(vec, others):
    M = []
    for v in vec:
        M.append(np.zeros(3))
        for u2, v2 in zip(vec, others):
            M[-1] = M[-1] + get_torque(v, u2) + get_torque(v, v2)
    return M


def get_torque(vec1, vec2):
    direction = get_normal_direction(vec1, vec2)
    return - direction*np.dot(vec1, vec2)


def get_normal_direction(vec1, vec2):
    direction = np.cross(vec1, vec2)
    norm = np.linalg.norm(direction)

    if norm < 1e-12:
        direction = np.cross(np.random.rand(3), vec1)
        norm = np.linalg.norm(direction)

    return direction/(norm+1e-12)


def get_initial_directions(num_points):
    alpha = np.random.rand(num_points)
    beta = np.random.rand(num_points)
    gamma = np.random.rand(num_points)
#    alpha[0] = 0
 #   beta[0] = 0
  #  gamma[0] = 0
    u = []
    v = []
    for a,b,g in zip(alpha, beta, gamma):
        u0 = np.array([+np.cos(a), np.sin(a), 0])
        v0 = np.array([-np.sin(a), np.cos(a), 0])
        u.append(rotate_vector(u0, v0, b))
        v.append(rotate_vector(v0, u0, g))

    return u, v


def rotate_vector(vector, axis, angle):
    return quat.rotate_vectors(quat.from_rotation_vector(axis * angle), vector)


if __name__=='__main__':
  #  test()
    main(sys.argv)