# TODO: Documentation

def get_rand_vector_in_sphere():

    import numpy as np

    theta = 2 * np.pi * np.random.rand()
    phi   = np.arccos(1 - 2 * np.random.rand())
    x     = np.sin(phi) * np.cos(theta)
    y     = np.sin(phi) * np.sin(theta)
    z     = np.cos(phi)
    return np.array([x, y, z])


def get_angle_from_two_vectors(a, b):
    """
    In radians.
    """

    import numpy as np

    assert len(a) == 3, "Parameter 'a' should be an [X, Y, Z] vector."
    assert len(b) == 3, "Parameter 'b' should be an [X, Y, Z] vector."
    return np.arccos(np.dot(a, b)/(np.linalg.norm(a) * np.linalg.norm(b)))


def get_rotation_matrix_from_axis_angle(axis, angle):

    import numpy as np

    q0 = np.cos(0.5 * angle)
    q1, q2, q3 = np.sin(0.5 * angle) * axis / np.linalg.norm(axis)
    return np.array([[1-2*q2*q2-2*q3*q3, 2*q1*q2-2*q0*q3, 2*q1*q3+2*q0*q2],
        [2*q2*q1+2*q0*q3, 1-2*q3*q3-2*q1*q1, 2*q2*q3-2*q0*q1],
        [2*q3*q1-2*q0*q2, 2*q3*q2+2*q0*q1, 1-2*q1*q1-2*q2*q2]])


def rotate_coords_from_rotation_matrix(coords, rot_matrix, pivot):

    import numpy as np

    #1. Center coords on pivot
    coords1 = coords - pivot
    #2. Apply rotation matrix
    coords2 = np.matmul(rot_matrix, coords1.transpose())
    #3. Recenter coords in original position
    return coords2.transpose() + pivot 