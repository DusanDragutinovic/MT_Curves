#!/usr/bin/sage
import os
import pickle
import time

from multiprocessing import Pool, Manager

import sage.all
from sage.all import *

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

"""
In this code, we want to exclude the curves that are isomorphic via a Stab[1] matrix
"""

"""
potential_list_of_third_quadrics, second_quadric_list, OP, 
                  Stab, StabP_cap_StabQ, matrices_to_P
"""


with open(os.path.join('data', 'prepared_data.pkl'), 'rb') as f:
    _, second_quadric_list, OP, Stab, _, matrices_to_P = pickle.load(f)
var('x, y, z, t, u')
poly_ring = GF(2)['X,Y,Z,T,U']
X,Y,Z,T,U = poly_ring.gens()

StabP = Stab[4].copy()

P1 = (1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0) 
P3 = (0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0) 
P4 = (0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0) 
P6 = (0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0)

"""
Initializing the data;


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!change only P4 -> Pi and P[4] to P[i]!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
"""
#initial_P4_list is the one obtained after main_calculation.py
initial_P4_list = [...]

P4_list = {}
P4_list_vspaces = {}
for i in range(1, 8):
    P4_list[i] = []
    P4_list_vspaces[i] = []

V = VectorSpace(GF(2),15)
for triple_of_quadrics in initial_P4_list:
    P, Q, R = triple_of_quadrics
    W = V.subspace([P, Q, R])
    how_many_in_OP4 = 0
    for w in W:
        w_with_int_entries = (int(w[0]), int(w[1]), int(w[2]), int(w[3]), int(w[4]), int(w[5]), int(w[6]), int(w[7]), int(w[8]), int(w[9]), int(w[10]), int(w[11]), int(w[12]), int(w[13]), int(w[14]))
        if w_with_int_entries in OP[4]:
            how_many_in_OP4 = how_many_in_OP4 + 1
    P4_list[how_many_in_OP4].append([P, Q, R])
    P4_list_vspaces[how_many_in_OP4].append(W)

isomph_classes_reprs = []

#Excluding the duplicates and finding the representatives!

for i in range(1, 8):
    print(i)
    while P4_list[i] != []:
        print(len(P4_list[i]))
        P, Q, R = (P4_list[i])[0]
        W = V.subspace([P, Q, R])
        isomph_classes_reprs.append([P, Q, R])
        
        #first remove all appearances of W
        while W in P4_list_vspaces[i]:
            ind_of_W = (P4_list_vspaces[i]).index(W)
            (P4_list[i]).remove((P4_list[i])[ind_of_W])
            (P4_list_vspaces[i]).remove((P4_list_vspaces[i])[ind_of_W])
        
        #now removing the isomorphic images
        for w in W:
            w_with_int_entries = (int(w[0]), int(w[1]), int(w[2]), int(w[3]), int(w[4]), int(w[5]), int(w[6]), int(w[7]), int(w[8]), int(w[9]), int(w[10]), int(w[11]), int(w[12]), int(w[13]), int(w[14]))
            
            if w_with_int_entries == P4:
                for m in StabP:
                    P_after_m, Q_after_m, R_after_m = switch_using_m(P, Q, R, m)
                    W_after_m = V.subspace([GF2coef(P_after_m), GF2coef(Q_after_m), GF2coef(R_after_m)])
                    while W_after_m in P4_list_vspaces[i]:
                        ind_of_W_a_m = (P4_list_vspaces[i]).index(W_after_m)
                        (P4_list[i]).remove((P4_list[i])[ind_of_W_a_m])
                        (P4_list_vspaces[i]).remove((P4_list_vspaces[i])[ind_of_W_a_m])
            elif w_with_int_entries in OP[4]:
                ind_in_OP4 = (OP[4]).index(w_with_int_entries)
                m = (matrices_to_P[4])[ind_in_OP4]
                P_mapping_w_to_P4, Q_mapping_w_to_P4, R_mapping_w_to_P4 = switch_using_m(P, Q, R, m)
                for n in StabP:
                    P_mn, Q_mn, R_mn = switch_using_m(P_mapping_w_to_P4, Q_mapping_w_to_P4, R_mapping_w_to_P4, n)
                    W_mn = V.subspace([GF2coef(P_mn), GF2coef(Q_mn), GF2coef(R_mn)])
                    while W_mn in P4_list_vspaces[i]:
                        ind_of_W_mn = (P4_list_vspaces[i]).index(W_mn)
                        (P4_list[i]).remove((P4_list[i])[ind_of_W_mn])
                        (P4_list_vspaces[i]).remove((P4_list_vspaces[i])[ind_of_W_mn])

print(len(isomph_classes_reprs))

path = os.path.join("isomorphism_classes")
os.makedirs(path, exist_ok = True)
with open(os.path.join(path, "P4_list.txt"), "w") as f:
    f.write(str(isomph_classes_reprs))
print("[INFO] Done")