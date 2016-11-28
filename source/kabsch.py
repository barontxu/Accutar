"""

https://github.com/charnley/rmsd
https://en.wikipedia.org/wiki/Kabsch_algorithm

"""

import numpy as np
import re

def kabsched_Q(Q_origin, P_origin):
    """
    Rotate matrix Q unto P and return Q
    Q and P are numpy array with shape of (N * 3)
    """

    P = np.copy(P_origin)
    Q = np.copy(Q_origin)

    Pc = centroid(P)
    Qc = centroid(Q)

    P -= Pc
    Q -= Qc

    U = kabsch(Q, P)

    # Rotate Q
    Q = np.dot(Q, U)

    Q += Pc

    return Q


def kabsch_rmsd(Q_origin, P_origin):
    """
    Rotate matrix Q unto P and calculate the RMSD
    """
    P = np.copy(P_origin)
    Q = np.copy(Q_origin)

    Pc = centroid(P)
    Qc = centroid(Q)

    P -= Pc
    Q -= Qc

    U = kabsch(Q, P)

    # Rotate Q
    Q = np.dot(Q, U)

    return rmsd(Q, P)


def rotate(Q, P):
    """
    Rotate matrix Q unto matrix P using Kabsch algorithm
    """
    U = kabsch(Q, P)

    # Rotate P
    Q = np.dot(Q, U)
    return Q


def kabsch(Q, P):
    """
    The optimal rotation matrix U is calculated and then used to rotate matrix
    Q unto matrix P so the minimum root-mean-square deviation (RMSD) can be
    calculated.

    Using the Kabsch algorithm with two sets of paired point Q and P,
    centered around the center-of-mass.
    Each vector set is represented as an NxD matrix, where D is the
    the dimension of the space.

    The algorithm works in three steps:
    - a translation of Q and P
    - the computation of a covariance matrix C
    - computation of the optimal rotation matrix U

    http://en.wikipedia.org/wiki/Kabsch_algorithm

    Parameters:
    Q -- (N, number of points)x(D, dimension) matrix
    P -- (N, number of points)x(D, dimension) matrix

    Returns:
    U -- Rotation matrix

    """

    # Computation of the covariance matrix
    C = np.dot(np.transpose(Q), P)

    # Computation of the optimal rotation matrix
    # This can be done using singular value decomposition (SVD)
    # Getting the sign of the det(V)*(W) to decide
    # whether we need to correct our rotation matrix to ensure a
    # right-handed coordinate system.
    # And finally calculating the optimal rotation matrix U
    # see http://en.wikipedia.org/wiki/Kabsch_algorithm
    V, S, W = np.linalg.svd(C)
    d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0

    if d:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]

    # Create Rotation matrix U
    U = np.dot(V, W)

    return U


def centroid(X):
    """
    Calculate the centroid from a vectorset X
    """
    C = sum(X)/len(X)
    return C


def rmsd(V, W):
    """
    Calculate Root-mean-square deviation from two sets of vectors V and W.
    """
    D = len(V[0])
    N = len(V)
    rmsd = 0.0
    for v, w in zip(V, W):
        rmsd += sum([(v[i]-w[i])**2.0 for i in range(D)])
    return np.sqrt(rmsd/N)
