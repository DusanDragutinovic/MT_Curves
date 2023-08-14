#!/usr/bin/sage
import os
import pickle
import time
import numpy as np

from multiprocessing import Pool, Manager

import sage.all
from sage.all import *

with open(os.path.join('data', 'prepared_data.pkl'), 'rb') as f:
    _, _, OP, Stab, _, matrices_to_P = pickle.load(f)

print("Started")
def GF2coef(P):
    return (GF(2)(P[0]), GF(2)(P[1]), GF(2)(P[2]), GF(2)(P[3]), GF(2)(P[4]), GF(2)(P[5]), GF(2)(P[6]), GF(2)(P[7]), GF(2)(P[8]), GF(2)(P[9]), GF(2)(P[10]), GF(2)(P[11]), GF(2)(P[12]), GF(2)(P[13]), GF(2)(P[14]))


def switch_using_m(P, Q, R, m):
    var('x, y, z, t, u')
    poly_ring = GF(2)['X,Y,Z,T,U']
    X,Y,Z,T,U = poly_ring.gens()
    def fP(x, y, z, t, u):
        return P[0]*x**2 + P[1]*x*y + P[2]*x*z + P[3]*x*t + P[4]*x*u + P[5]*y**2 + P[6]*y*z + P[7]*y*t + P[8]*y*u + P[9]*z**2 + P[10]*z*t + P[11]*z*u + P[12]*t**2 + P[13]*t*u + P[14]*u**2
    def fQ(x, y, z, t, u):
        return Q[0]*x**2 + Q[1]*x*y + Q[2]*x*z + Q[3]*x*t + Q[4]*x*u + Q[5]*y**2 + Q[6]*y*z + Q[7]*y*t + Q[8]*y*u + Q[9]*z**2 + Q[10]*z*t + Q[11]*z*u + Q[12]*t**2 + Q[13]*t*u + Q[14]*u**2
    def fR(x, y, z, t, u):
        return R[0]*x**2 + R[1]*x*y + R[2]*x*z + R[3]*x*t + R[4]*x*u + R[5]*y**2 + R[6]*y*z + R[7]*y*t + R[8]*y*u + R[9]*z**2 + R[10]*z*t + R[11]*z*u + R[12]*t**2 + R[13]*t*u + R[14]*u**2
    def gP(x, y, z, t, u):
        return fP(m[0, 0]*x + m[0, 1]*y + m[0, 2]*z + m[0, 3]*t + m[0, 4]*u, m[1, 0]*x + m[1, 1]*y + m[1, 2]*z + m[1, 3]*t + m[1, 4]*u, m[2, 0]*x + m[2, 1]*y + m[2, 2]*z + m[2, 3]*t + m[2, 4]*u, m[3, 0]*x + m[3, 1]*y + m[3, 2]*z + m[3, 3]*t + m[3, 4]*u, m[4, 0]*x + m[4, 1]*y + m[4, 2]*z + m[4, 3]*t + m[4, 4]*u)
    def gQ(x, y, z, t, u):
        return fQ(m[0, 0]*x + m[0, 1]*y + m[0, 2]*z + m[0, 3]*t + m[0, 4]*u, m[1, 0]*x + m[1, 1]*y + m[1, 2]*z + m[1, 3]*t + m[1, 4]*u, m[2, 0]*x + m[2, 1]*y + m[2, 2]*z + m[2, 3]*t + m[2, 4]*u, m[3, 0]*x + m[3, 1]*y + m[3, 2]*z + m[3, 3]*t + m[3, 4]*u, m[4, 0]*x + m[4, 1]*y + m[4, 2]*z + m[4, 3]*t + m[4, 4]*u)
    def gR(x, y, z, t, u):
        return fR(m[0, 0]*x + m[0, 1]*y + m[0, 2]*z + m[0, 3]*t + m[0, 4]*u, m[1, 0]*x + m[1, 1]*y + m[1, 2]*z + m[1, 3]*t + m[1, 4]*u, m[2, 0]*x + m[2, 1]*y + m[2, 2]*z + m[2, 3]*t + m[2, 4]*u, m[3, 0]*x + m[3, 1]*y + m[3, 2]*z + m[3, 3]*t + m[3, 4]*u, m[4, 0]*x + m[4, 1]*y + m[4, 2]*z + m[4, 3]*t + m[4, 4]*u)
    hP = poly_ring(gP(X, Y, Z, T, U))
    nP = (hP.coefficient(X**2), hP.coefficient(X*Y), hP.coefficient(X*Z), hP.coefficient(X*T), hP.coefficient(X*U), hP.coefficient(Y**2), hP.coefficient(Y*Z), hP.coefficient(Y*T), hP.coefficient(Y*U), hP.coefficient(Z**2), hP.coefficient(Z*T), hP.coefficient(Z*U), hP.coefficient(T**2), hP.coefficient(T*U), hP.coefficient(U**2))
    hQ = poly_ring(gQ(X, Y, Z, T, U))
    nQ = (hQ.coefficient(X**2), hQ.coefficient(X*Y), hQ.coefficient(X*Z), hQ.coefficient(X*T), hQ.coefficient(X*U), hQ.coefficient(Y**2), hQ.coefficient(Y*Z), hQ.coefficient(Y*T), hQ.coefficient(Y*U), hQ.coefficient(Z**2), hQ.coefficient(Z*T), hQ.coefficient(Z*U), hQ.coefficient(T**2), hQ.coefficient(T*U), hQ.coefficient(U**2))
    hR = poly_ring(gR(X, Y, Z, T, U))
    nR = (hR.coefficient(X**2), hR.coefficient(X*Y), hR.coefficient(X*Z), hR.coefficient(X*T), hR.coefficient(X*U), hR.coefficient(Y**2), hR.coefficient(Y*Z), hR.coefficient(Y*T), hR.coefficient(Y*U), hR.coefficient(Z**2), hR.coefficient(Z*T), hR.coefficient(Z*U), hR.coefficient(T**2), hR.coefficient(T*U), hR.coefficient(U**2))
    el = [nP, nQ, nR]
    return el


P1 = (1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0) 
P3 = (0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0) 
P4 = (0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0) 
P6 = (0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0)


P4_list = [...]
P4_automorphisms = {}

for i in range(0, len(P4_list)):
    P4_automorphisms[i] = []

V = VectorSpace(GF(2),15)
for i in range(0, len(P4_list)):
    P, Q, R = P4_list[i]
    W = V.subspace([P, Q, R])
    for w in W: 
        w_with_int_entries = (int(w[0]), int(w[1]), int(w[2]), int(w[3]), int(w[4]), int(w[5]), int(w[6]), int(w[7]), int(w[8]), int(w[9]), int(w[10]), int(w[11]), int(w[12]), int(w[13]), int(w[14]))
        if w_with_int_entries in OP[4]:
            ind = OP[4].index(w_with_int_entries)
            #n maps w to P4
            n = (matrices_to_P[4])[ind]
            P_after_n, Q_after_n, R_after_n = switch_using_m(P, Q, R, n)
            for m in Stab[4]:
                P_after_nm, Q_after_nm, R_after_nm = switch_using_m(P_after_n, Q_after_n, R_after_n, m)
                W_after_nm = V.subspace([GF2coef(P_after_nm), GF2coef(Q_after_nm), GF2coef(R_after_nm)])
                if W_after_nm == W:
                   P4_automorphisms[i].append(n*m)
    print(i, len(P4_automorphisms[i]))

weighted_sum = 0
lengths_of_autogrps = []

for i in range(0, len(P4_list)):
    auto_size = len(P4_automorphisms[i])
    lengths_of_autogrps.append(auto_size)
    weighted_sum = weighted_sum + QQ(1/auto_size)
print(weighted_sum)

    
print("done")

path = os.path.join("automorphisms_data")
os.makedirs(path, exist_ok = True)
with open(os.path.join(path, "automorphisms_P4.txt"), "w") as f:
    f.write(str(P4_automorphisms))
with open(os.path.join(path, "automorphisms_P4_lengths.txt"), "w") as f:
    f.write(str(lengths_of_autogrps))
    