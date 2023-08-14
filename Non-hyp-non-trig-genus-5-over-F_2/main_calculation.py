#!/usr/bin/sage

import os
import pickle
import time
import sys

from multiprocessing import Pool

import sage.all
from sage.all import *

PARALLEL = True

R_list_PQi = None
Falg = GF(2).algebraic_closure()
Pr2 = ProjectiveSpace(GF(4), 4)
Pr3 = ProjectiveSpace(GF(8), 4)
Pr = ProjectiveSpace(Falg, 4)    

def iteration(data):
    poly_g, i, j = data
    if j % 100 == 0:
        print("i={}, j={} done".format(i, j))
    R = R_list_PQi[i][j]
    poly_ring = GF(2)['X,Y,Z,T,U']
    X,Y,Z,T,U = poly_ring.gens()
    
    poly_h = R[0]*X ** 2 + R[1]*X*Y + R[2]*X*Z + R[3]*X*T + R[4]*X*U + R[5]*Y ** 2 + R[6]*Y*Z + R[7]*Y*T + R[8]*Y*U + R[9]*Z ** 2 + R[10]*Z*T + R[11]*Z*U + R[12]*T ** 2 + R[13]*T*U + R[14]*U ** 2
    
    #checking that all nonzero elements in the vector space
    #<q_P, q_Q, q_R> are irreducible
    ind = 0
    for coef_f in GF(2):
        for coef_g in GF(2):
            for coef_h in GF(2):
                test_poly = poly_ring(coef_f*poly_f + coef_g*poly_g + coef_h*poly_h)
                if test_poly.is_prime() == False:
                    ind = ind + 1
    if ind == 1:
        h1 = poly_f + poly_h
        n1 = (h1.coefficient(X**2), h1.coefficient(X*Y), h1.coefficient(X*Z), h1.coefficient(X*T), h1.coefficient(X*U), h1.coefficient(Y**2), h1.coefficient(Y*Z), h1.coefficient(Y*T), h1.coefficient(Y*U), h1.coefficient(Z**2), h1.coefficient(Z*T), h1.coefficient(Z*U), h1.coefficient(T**2), h1.coefficient(T*U), h1.coefficient(U**2))
        
        h2 = poly_g + poly_h
        n2 = (h2.coefficient(X**2), h2.coefficient(X*Y), h2.coefficient(X*Z), h2.coefficient(X*T), h2.coefficient(X*U), h2.coefficient(Y**2), h2.coefficient(Y*Z), h2.coefficient(Y*T), h2.coefficient(Y*U), h2.coefficient(Z**2), h2.coefficient(Z*T), h2.coefficient(Z*U), h2.coefficient(T**2), h2.coefficient(T*U), h2.coefficient(U**2))
        
        h3 = poly_f + poly_g + poly_h
        n3 = (h3.coefficient(X**2), h3.coefficient(X*Y), h3.coefficient(X*Z), h3.coefficient(X*T), h3.coefficient(X*U), h3.coefficient(Y**2), h3.coefficient(Y*Z), h3.coefficient(Y*T), h3.coefficient(Y*U), h3.coefficient(Z**2), h3.coefficient(Z*T), h3.coefficient(Z*U), h3.coefficient(T**2), h3.coefficient(T*U), h3.coefficient(U**2))
        
        if (n1 in potential_list_of_third_quadrics) and (n2 in potential_list_of_third_quadrics) and (n3 in potential_list_of_third_quadrics): 
            J = poly_ring.ideal(poly_f, poly_g, poly_h)
            Cur2 = Pr2.subscheme(J)
            if Cur2.is_smooth():
                Cur3 = Pr3.subscheme(J)
                if Cur3.is_smooth():    
                    Cur = Pr.subscheme(J)
                    if Cur.is_smooth():
                        path = os.path.join("data", "smooth_curves")
                        os.makedirs(path, exist_ok = True)
                        with open(os.path.join(path, "P4_list.txt"), "a+") as f:
                            f.write("[{}, {} , {}], ".format(P, Q, R))

# data/calc/{i}/{j}.txt << {R}

"""
potential_list_of_third_quadrics, second_quadric_list, OP, 
                  Stab, StabP_cap_StabQ, matrices_to_P
"""

if __name__ == "__main__":
    try:
        with open(os.path.join('data', 'prepared_data.pkl'), 'rb') as f:
                potential_list_of_third_quadrics, second_quadric_list, _, _, _, _= pickle.load(f)
    except Exception as e:
        print(e)
        print("You need to run `prepare_initial_data.py` first!")
        exit(1)
    
    try:
        with open(os.path.join('data', 'Rs.pkl'), 'rb') as f:
                R_list_PQi = pickle.load(f)
    except Exception as e:
        print(e)
        print("You need to run `find_r_candidates.py` first!")
        exit(1)
    
        
    P1 = (1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0) 
    P3 = (0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0) 
    P4 = (0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0) 
    P6 = (0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0)

    
    P = P4
    var('x, y, z, t, u')
    poly_ring = GF(2)['X,Y,Z,T,U']
    X,Y,Z,T,U = poly_ring.gens()

    def f(x, y, z, t, u):
        return P[0]*x**2 + P[1]*x*y + P[2]*x*z + P[3]*x*t + P[4]*x*u + P[5]*y**2 + P[6]*y*z + P[7]*y*t + P[8]*y*u + P[9]*z**2 + P[10]*z*t + P[11]*z*u + P[12]*t**2 + P[13]*t*u + P[14]*u**2
    poly_f = poly_ring(f(X, Y, Z, T, U))
    
    UPPER_BOUND = len(second_quadric_list)
    for i in range(0, UPPER_BOUND):
        Q = second_quadric_list[i]
        poly_g = Q[0]*X ** 2 + Q[1]*X*Y + Q[2]*X*Z + Q[3]*X*T + Q[4]*X*U + Q[5]*Y ** 2 + Q[6]*Y*Z + Q[7]*Y*T + Q[8]*Y*U + Q[9]*Z ** 2 + Q[10]*Z*T + Q[11]*Z*U + Q[12]*T ** 2 + Q[13]*T*U + Q[14]*U ** 2
        start = time.time()

        if PARALLEL:
            with Pool() as pool:
                pool.map(iteration, [(poly_g, i, j) for j in range(len(R_list_PQi[i]))])
        else:
            for j in range(len(R_list_PQi[i])):
                iteration((poly_g, i, j))
        print("[PROGRESS] Done {}%".format(100 * (i+1) / len(second_quadric_list)))

        end = time.time()
        print("[TIMER] Last iteration (for i={}) done in {}".format(i, end - start))
    print("[INFO] Done")
