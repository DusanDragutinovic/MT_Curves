#!/usr/bin/sage
import os
import pickle
import time

from multiprocessing import Pool, Manager

import sage.all
from sage.all import *

candidates = []
lisPr = []
manager = Manager()
R_i_cands = manager.dict()

PARALLEL = True
LOWER_BOUND = 0
UPPER_BOUND = 100

def iteration(i):
    G = candidates[i]
    L = lisPr.copy()
    R = []
    start = time.time()
    j = 0
    
    while L != []:
        P = L[0]
        R.append(P)
        exclude = [P]
        def f(x, y, z, t, u):
            return P[0]*x**2 + P[1]*x*y + P[2]*x*z + P[3]*x*t + P[4]*x*u + P[5]*y**2 + P[6]*y*z + P[7]*y*t + P[8]*y*u + P[9]*z**2 + P[10]*z*t + P[11]*z*u + P[12]*t**2 + P[13]*t*u + P[14]*u**2
        for m in G:
            def g(x, y, z, t, u):
                return f(m[0, 0]*x + m[0, 1]*y + m[0, 2]*z + m[0, 3]*t + m[0, 4]*u, m[1, 0]*x + m[1, 1]*y + m[1, 2]*z + m[1, 3]*t + m[1, 4]*u, m[2, 0]*x + m[2, 1]*y + m[2, 2]*z + m[2, 3]*t + m[2, 4]*u, m[3, 0]*x + m[3, 1]*y + m[3, 2]*z + m[3, 3]*t + m[3, 4]*u, m[4, 0]*x + m[4, 1]*y + m[4, 2]*z + m[4, 3]*t + m[4, 4]*u)
            h = g(X, Y, Z, T, U)
            Q = (h.coefficient(X**2), h.coefficient(X*Y), h.coefficient(X*Z), h.coefficient(X*T), h.coefficient(X*U), h.coefficient(Y**2), h.coefficient(Y*Z), h.coefficient(Y*T), h.coefficient(Y*U), h.coefficient(Z**2), h.coefficient(Z*T), h.coefficient(Z*U), h.coefficient(T**2), h.coefficient(T*U), h.coefficient(U**2))
            exclude.append(Q)
        j+=1
        L = list(set(L) - set(exclude))
        print("For i={}, done iteration {}/{}".format(i, j, len(L)))
    
    end = time.time()
    print("[TIMER] i={} done in {}".format(i, end - start))
    R_i_cands[i] = R

if __name__ == "__main__":
    try:
        with open(os.path.join('data', 'prepared_data.pkl'), 'rb') as f:
            lisPr, Sec, Stab, candidates = pickle.load(f)
    except Exception as e:
        print(e)
        print("You need to run `prepare_initial_data.py` first!")
        exit(1)

    var('x, y, z, t, u')
    teri = GF(2)['x,y,z,t,u']
    X,Y,Z,T,U = teri.gens()

    start = time.time()
    if PARALLEL:
        with Pool() as pool:
            pool.map(iteration, [i for i in range(len(Sec))])
    else:
        for i in range(len(Sec)):
            iteration((g, i, j))
        print("[PROGRESS] Done {}%".format(100 * (i+1) / len(Sec)))
   
    with open(os.path.join("data", "Rs.pkl"), "wb") as f:
        pickle.dump(dict(R_i_cands), f)

    print("[INFO] Done")
