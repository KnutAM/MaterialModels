import numpy as np
import matplotlib.pyplot as plt
import sys


def main(argv):

    miller_planes, miller_directions = get_crystal_info('bcc')
    planes, directions = get_combined_systems(miller_planes, miller_directions)

    plot_yield_functions(planes, directions)

    # print_system_to_mathematica(planes, directions)
    # print_system_to_matlab(planes, directions)


def plot_yield_functions(planes, directions):
    num_dir = 100
    angle = np.linspace(0, 2 * np.pi, num_dir)
    sigma = np.cos(angle)
    tau = np.sin(angle) / np.sqrt(3)
    plt.figure()


    tau_alpha_max = np.zeros(num_dir)
    for p, d in zip(planes, directions):
        tau_alpha = get_resolved_stress(sigma, tau, p, d)
        plt.plot(angle, tau_alpha,
                 label='p=[' + str(p[0]) + str(p[1]) + str(p[2]) + '], d=[' + str(d[0]) + str(d[1]) + str(d[2]) + ']')
        tau_alpha_max = np.maximum(tau_alpha_max, np.abs(tau_alpha))

    plt.plot(angle, tau_alpha_max, color='black', linestyle='--')
    plt.legend(loc='upper right')
    plt.show()



def get_resolved_stress(sigma, tau, plane, direction):
    tau_alpha = []
    for s, t in zip(sigma, tau):
        smat = np.array([[s, 0, 0], [0, 0, t], [0, t, 0]])
        tau_alpha.append(np.dot(np.matmul(plane, smat), direction) / (
                    np.linalg.norm(plane) * np.linalg.norm(direction)))

    return tau_alpha


def print_system_to_matlab(planes, directions):
    print('Slip planes')
    print_mathematica_vectors(planes)
    print('Slip directions')
    print_mathematica_vectors(directions)
    print('End')


def print_matlab_vectors(vectors):
    print('[', end='')
    print_mathematica_vector_norm(vectors[0])
    for v in vectors[1:]:
        print(', ', end='')
        print_mathematica_vector_norm(v)
    print(']')


def print_matlab_vector(vector):
    print('[', end='')
    print(str(vector[0]), end='')
    for v in vector[1:]:
        print('; ', str(v), end='')
    print(']', end='')


def print_matlab_vector_norm(vector):
    print('[', end='')
    print(str(vector[0]), end='')
    for v in vector[1:]:
        print('; ', str(v), end='')
    print(']/sqrt(', end='')
    print(sum(np.power(vector, 2)), end='')
    print(')', end='')



def print_system_to_mathematica(planes, directions):
    print('{')
    print_mathematica_vectors(planes)
    print(',')
    print_mathematica_vectors(directions)
    print('}')


def print_mathematica_vectors(vectors):
    print('{', end='')
    print_mathematica_vector_norm(vectors[0])
    for v in vectors[1:]:
        print(', ', end='')
        print_mathematica_vector_norm(v)
    print('}')


def print_mathematica_vector(vector):
    print('{', end='')
    print(str(vector[0]), end='')
    for v in vector[1:]:
        print(', ', str(v), end='')
    print('}', end='')


def print_mathematica_vector_norm(vector):
    print('{', end='')
    print(str(vector[0]), end='')
    for v in vector[1:]:
        print(', ', str(v), end='')
    print('}/Sqrt[', end='')
    print(sum(np.power(vector, 2)), end='')
    print(']', end='')


def get_crystal_info(crystal):
    if crystal == 'fcc':
        miller_planes = [[1, 1, 1]]
        miller_directions = [[1, 1, 0]]
    elif crystal == 'bcc':
        miller_planes = [[1, 1, 0], [1, 1, 2], [1, 2, 3]]
        miller_directions = [[1, 1, 1], [1, 1, 1], [1, 1, 1]]
    elif crystal == 'bcc12':
        miller_planes = [[1, 1, 0]]
        miller_directions = [[1, 1, 1]]
    else:
        print('Crystal ' + crystal + ' is not supported')
        miller_planes = None
        miller_directions = None

    return miller_planes, miller_directions


def normalize_vectors(vectors):
    vnew = []
    for v in vectors:
        vnew.append(v/np.linalg.norm(v))

    return vnew


def get_combined_systems(miller_planes, miller_directions):
    planes = []
    directions = []
    for mp, md in zip(miller_planes, miller_directions):
        p = get_permutations(mp)
        p, d = get_systems(p, md)
        for ip, id in zip(p, d):
            planes.append(ip)
            directions.append(id)

    return planes, directions


def get_systems(planes, miller):
    all_directions = get_permutations(miller)
    planes_new = []
    directions_new = []
    for p in planes:
        dirs = get_slip_directions(p, all_directions)
        for d in dirs:
            planes_new.append(p)
            directions_new.append(d)

    return planes_new, directions_new


def get_slip_directions(plane, directions):
    dirs = []
    for d in directions:
        if np.dot(d, plane) == 0:
            dirs.append(d)
    return dirs


def get_permutations(miller):
    # Given a set of miller indicies, flip these to return the set of unique indicies
    perm = []
    ind_perm = ([0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0])
    sgn_perm = np.array(([1, 1, 1], [1, 1, -1], [1, -1, 1], [1, -1, -1]))

    for inds in ind_perm:
        mperm = np.array([miller[ind] for ind in inds])
        for sgn in sgn_perm:
            perm.append(mperm*sgn)

    perm = remove_duplicates(perm)

    return perm


def remove_duplicates(perm):
    pnew = []
    pnew.append(perm[0])
    for p in perm[1:]:
        if len(pnew) == 1:
            if not ((pnew == p).all() or (pnew == (-p)).all()):
                pnew.append(p)
        else:
            if not (any((pnew[:] == p).all(1)) or any((pnew[:] == -p).all(1))):
                pnew.append(p)

    return np.array(pnew)


if __name__=='__main__':
    main(sys.argv)
