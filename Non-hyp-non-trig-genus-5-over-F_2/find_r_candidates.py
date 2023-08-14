#!/usr/bin/sage
import os
import pickle
import time

from multiprocessing import Pool, Manager

import sage.all
from sage.all import *

StabP_cap_StabQ = []
potential_list_of_third_quadrics = []
manager = Manager()
R_list_PQi = manager.dict()

PARALLEL = True


def iteration(i):
    var('x, y, z, t, u')
    poly_ring = GF(2)['X,Y,Z,T,U']
    X,Y,Z,T,U = poly_ring.gens()    
    
    G = StabP_cap_StabQ[i]
    L = potential_list_of_third_quadrics.copy()
    temp_R_list_PQi = []
    
    start = time.time()
    j = 0
    
    P1 = (1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0) 
    P3 = (0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0) 
    P4 = (0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0) 
    P6 = (0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0)
    
    """
    P4 computation
    """
    P = P4 
    Q = second_quadric_list[i]
       
    while L != []:
        R = L[0]
        temp_R_list_PQi.append(R)
        exclude = [R]
        """
        a fixed R defines the same ideal as  R + P, R + Q, R + Q + P, 
        and it defines the isomorphic curve when we go through all of the 
        G = StabP_cap_StabQ[i] = Stab[1]\cap Stab(second_quadric_list[i])-conjugates 
        of R, R + P, R + Q, R + Q + P, 
        """
        def fR(x, y, z, t, u):
            return R[0]*x**2 + R[1]*x*y + R[2]*x*z + R[3]*x*t + R[4]*x*u + R[5]*y**2 + R[6]*y*z + R[7]*y*t + R[8]*y*u + R[9]*z**2 + R[10]*z*t + R[11]*z*u + R[12]*t**2 + R[13]*t*u + R[14]*u**2

        def fR_plus_P(x, y, z, t, u):
            return (R[0] + P[0])*x**2 + (R[1] + P[1])*x*y + (R[2] + P[2])*x*z + (R[3] + P[3])*x*t + (R[4] + P[4])*x*u + (R[5] + P[5])*y**2 + (R[6] + P[6])*y*z + (R[7] + P[7])*y*t + (R[8] + P[8])*y*u + (R[9] + P[9])*z**2 + (R[10] + P[10])*z*t + (R[11] + P[11])*z*u + (R[12] + P[12])*t**2 + (R[13] + P[13])*t*u + (R[14] + P[14])*u**2

        def fR_plus_Q(x, y, z, t, u):
            return (R[0] + Q[0])*x**2 + (R[1] + Q[1])*x*y + (R[2] + Q[2])*x*z + (R[3] + Q[3])*x*t + (R[4] + Q[4])*x*u + (R[5] + Q[5])*y**2 + (R[6] + Q[6])*y*z + (R[7] + Q[7])*y*t + (R[8] + Q[8])*y*u + (R[9] + Q[9])*z**2 + (R[10] + Q[10])*z*t + (R[11] + Q[11])*z*u + (R[12] + Q[12])*t**2 + (R[13] + Q[13])*t*u + (R[14] + Q[14])*u**2

        def fR_plus_P_plus_Q(x, y, z, t, u):
            return (R[0] + Q[0] + P[0])*x**2 + (R[1] + Q[1] + P[1])*x*y + (R[2] + Q[2] + P[2])*x*z + (R[3] + Q[3] + P[3])*x*t + (R[4] + Q[4] + P[4])*x*u + (R[5] + Q[5] + P[5])*y**2 + (R[6] + Q[6] + P[6])*y*z + (R[7] + Q[7] + P[7])*y*t + (R[8] + Q[8] + P[8])*y*u + (R[9] + Q[9] + P[9])*z**2 + (R[10] + Q[10] + P[10])*z*t + (R[11] + Q[11] + P[11])*z*u + (R[12] + Q[12] + P[12])*t**2 + (R[13] + Q[13] + P[13])*t*u + (R[14] + Q[14] + P[14])*u**2
        
        #definitions above are correct
        for m in G:
            def gR(x, y, z, t, u):
                return fR(m[0, 0]*x + m[0, 1]*y + m[0, 2]*z + m[0, 3]*t + m[0, 4]*u, m[1, 0]*x + m[1, 1]*y + m[1, 2]*z + m[1, 3]*t + m[1, 4]*u, m[2, 0]*x + m[2, 1]*y + m[2, 2]*z + m[2, 3]*t + m[2, 4]*u, m[3, 0]*x + m[3, 1]*y + m[3, 2]*z + m[3, 3]*t + m[3, 4]*u, m[4, 0]*x + m[4, 1]*y + m[4, 2]*z + m[4, 3]*t + m[4, 4]*u)    
            poly_gR = poly_ring(gR(X, Y, Z, T, U))
            R_conj = (poly_gR.coefficient(X**2), poly_gR.coefficient(X*Y), poly_gR.coefficient(X*Z), poly_gR.coefficient(X*T), poly_gR.coefficient(X*U), poly_gR.coefficient(Y**2), poly_gR.coefficient(Y*Z), poly_gR.coefficient(Y*T), poly_gR.coefficient(Y*U), poly_gR.coefficient(Z**2), poly_gR.coefficient(Z*T), poly_gR.coefficient(Z*U), poly_gR.coefficient(T**2), poly_gR.coefficient(T*U), poly_gR.coefficient(U**2))
            exclude.append(R_conj)
            
            def gR_plusP(x, y, z, t, u):
                return fR_plus_P(m[0, 0]*x + m[0, 1]*y + m[0, 2]*z + m[0, 3]*t + m[0, 4]*u, m[1, 0]*x + m[1, 1]*y + m[1, 2]*z + m[1, 3]*t + m[1, 4]*u, m[2, 0]*x + m[2, 1]*y + m[2, 2]*z + m[2, 3]*t + m[2, 4]*u, m[3, 0]*x + m[3, 1]*y + m[3, 2]*z + m[3, 3]*t + m[3, 4]*u, m[4, 0]*x + m[4, 1]*y + m[4, 2]*z + m[4, 3]*t + m[4, 4]*u)    
            poly_gR_plusP = poly_ring(gR_plusP(X, Y, Z, T, U))
            R_plusP_conj = (poly_gR_plusP.coefficient(X**2), poly_gR_plusP.coefficient(X*Y), poly_gR_plusP.coefficient(X*Z), poly_gR_plusP.coefficient(X*T), poly_gR_plusP.coefficient(X*U), poly_gR_plusP.coefficient(Y**2), poly_gR_plusP.coefficient(Y*Z), poly_gR_plusP.coefficient(Y*T), poly_gR_plusP.coefficient(Y*U), poly_gR_plusP.coefficient(Z**2), poly_gR_plusP.coefficient(Z*T), poly_gR_plusP.coefficient(Z*U), poly_gR_plusP.coefficient(T**2), poly_gR_plusP.coefficient(T*U), poly_gR_plusP.coefficient(U**2))
            exclude.append(R_plusP_conj)
            
            def gR_plusQ(x, y, z, t, u):
                return fR_plus_Q(m[0, 0]*x + m[0, 1]*y + m[0, 2]*z + m[0, 3]*t + m[0, 4]*u, m[1, 0]*x + m[1, 1]*y + m[1, 2]*z + m[1, 3]*t + m[1, 4]*u, m[2, 0]*x + m[2, 1]*y + m[2, 2]*z + m[2, 3]*t + m[2, 4]*u, m[3, 0]*x + m[3, 1]*y + m[3, 2]*z + m[3, 3]*t + m[3, 4]*u, m[4, 0]*x + m[4, 1]*y + m[4, 2]*z + m[4, 3]*t + m[4, 4]*u)    
            poly_gR_plusQ = poly_ring(gR_plusQ(X, Y, Z, T, U))
            R_plusQ_conj = (poly_gR_plusQ.coefficient(X**2), poly_gR_plusQ.coefficient(X*Y), poly_gR_plusQ.coefficient(X*Z), poly_gR_plusQ.coefficient(X*T), poly_gR_plusQ.coefficient(X*U), poly_gR_plusQ.coefficient(Y**2), poly_gR_plusQ.coefficient(Y*Z), poly_gR_plusQ.coefficient(Y*T), poly_gR_plusQ.coefficient(Y*U), poly_gR_plusQ.coefficient(Z**2), poly_gR_plusQ.coefficient(Z*T), poly_gR_plusQ.coefficient(Z*U), poly_gR_plusQ.coefficient(T**2), poly_gR_plusQ.coefficient(T*U), poly_gR_plusQ.coefficient(U**2))
            exclude.append(R_plusQ_conj)
            
            def gR_plus_P_plus_Q(x, y, z, t, u):
                return fR(m[0, 0]*x + m[0, 1]*y + m[0, 2]*z + m[0, 3]*t + m[0, 4]*u, m[1, 0]*x + m[1, 1]*y + m[1, 2]*z + m[1, 3]*t + m[1, 4]*u, m[2, 0]*x + m[2, 1]*y + m[2, 2]*z + m[2, 3]*t + m[2, 4]*u, m[3, 0]*x + m[3, 1]*y + m[3, 2]*z + m[3, 3]*t + m[3, 4]*u, m[4, 0]*x + m[4, 1]*y + m[4, 2]*z + m[4, 3]*t + m[4, 4]*u)    
            poly_gR_plus_P_plus_Q = poly_ring(gR_plus_P_plus_Q(X, Y, Z, T, U))
            R_plus_P_plus_Q_conj = (poly_gR_plus_P_plus_Q.coefficient(X**2), poly_gR_plus_P_plus_Q.coefficient(X*Y), poly_gR_plus_P_plus_Q.coefficient(X*Z), poly_gR_plus_P_plus_Q.coefficient(X*T), poly_gR_plus_P_plus_Q.coefficient(X*U), poly_gR_plus_P_plus_Q.coefficient(Y**2), poly_gR_plus_P_plus_Q.coefficient(Y*Z), poly_gR_plus_P_plus_Q.coefficient(Y*T), poly_gR_plus_P_plus_Q.coefficient(Y*U), poly_gR_plus_P_plus_Q.coefficient(Z**2), poly_gR_plus_P_plus_Q.coefficient(Z*T), poly_gR_plus_P_plus_Q.coefficient(Z*U), poly_gR_plus_P_plus_Q.coefficient(T**2), poly_gR_plus_P_plus_Q.coefficient(T*U), poly_gR_plus_P_plus_Q.coefficient(U**2))
            exclude.append(R_plus_P_plus_Q_conj)    
        j+=1
        #checking if this is okay: 
        L = list(set(L) - set(exclude))
        print("For i={}, done iteration {}/{}".format(i, j, len(L)))
    end = time.time()
    print("[TIMER] i={} done in {}".format(i, end - start))
    R_list_PQi[i] = temp_R_list_PQi.copy()

"""
potential_list_of_third_quadrics, second_quadric_list, OP, 
                  Stab, StabP_cap_StabQ
"""


if __name__ == "__main__":
    try:
        with open(os.path.join('data', 'prepared_data.pkl'), 'rb') as open_doc:
            potential_list_of_third_quadrics, second_quadric_list, _, _, StabP_cap_StabQ, _ = pickle.load(open_doc)
    except Exception as e:
        print(e)
        print("You need to run `prepare_initial_data.py` first!")
        exit(1)

    var('x, y, z, t, u')
    poly_ring = GF(2)['X,Y,Z,T,U']
    X,Y,Z,T,U = poly_ring.gens()

    start = time.time()
    if PARALLEL:
        with Pool() as pool:
            pool.map(iteration, [i for i in range(len(second_quadric_list))])
    else:
        for i in range(len(second_quadric_list)):
            iteration((g, i, j))
        print("[PROGRESS] Done {}%".format(100 * (i+1) / len(second_quadric_list)))
   
    with open(os.path.join("data", "Rs.pkl"), "wb") as f:
        pickle.dump(dict(R_list_PQi), f)

    print("[INFO] Done")
