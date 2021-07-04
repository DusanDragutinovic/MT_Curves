#!/usr/bin/sage

import os
import pickle
import time
import sys

from multiprocessing import Pool

import sage.all
from sage.all import *

PARALLEL = True
LOWER_BOUND = 110
UPPER_BOUND = 134
R_i_cands = None

def iteration(data):
    g, i, j = data
    R = R_i_cands[i][j]
    h = R[0]*X ** 2 + R[1]*X*Y + R[2]*X*Z + R[3]*X*T + R[4]*X*U + R[5]*Y ** 2 + R[6]*Y*Z + R[7]*Y*T + R[8]*Y*U + R[9]*Z ** 2 + R[10]*Z*T + R[11]*Z*U + R[12]*T ** 2 + R[13]*T*U + R[14]*U ** 2
    I = Ideal(f1, g, h)
    # prethodna linija koda naivno definise ideal I, dok sledeca omogucava njegovo razumevanje
    J = I.radical()
    # ispod je eliminacioni uslov
    if j % 100 == 0:
        print("i={}, j={} done".format(i, j))
    if I == J:
        if J.ngens() == 3:
            F2 = GF(2).algebraic_closure()
            Pr = ProjectiveSpace(F2, 4)
            Cur = Pr.subscheme(J)
            if Cur.is_smooth():
                path = os.path.join("data", "results", "{}".format(i))
                os.makedirs(path, exist_ok = True)
                with open(os.path.join(path, "{}.txt".format(j)), "w") as f:
                    f.write("[{} , {}]".format(Q, R))

# data/calc/{i}/{j}.txt << {R}


if __name__ == "__main__":
    try:
        with open(os.path.join('data', 'prepared_data.pkl'), 'rb') as f:
                _, Sec, _, _= pickle.load(f)
    except Exception as e:
        print(e)
        print("You need to run `prepare_initial_data.py` first!")
        exit(1)
    
    try:
        with open(os.path.join('data', 'Rs.pkl'), 'rb') as f:
                R_i_cands = pickle.load(f)
    except Exception as e:
        print(e)
        print("You need to run `find_r_candidates.py` first!")
        exit(1)

    P = (0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0)
    var('x, y, z, t, u')
    teri = GF(2)['x,y,z,t,u']
    X,Y,Z,T,U = teri.gens()

    lista = []
    def f(x, y, z, t, u):
            return P[0]*x ** 2 + P[1]*x*y + P[2]*x*z + P[3]*x*t + P[4]*x*u + P[5]*y ** 2 + P[6]*y*z + P[7]*y*t + P[8]*y*u + P[9]*z ** 2 + P[10]*z*t + P[11]*z*u + P[12]*t ** 2 + P[13]*t*u + P[14]*u ** 2
    f1 = f(X, Y, Z, T, U)
    #prolazimo kroz sve elemente liste Sec koja predstavlja F_2 u trojci (F_1, F_2, F_3) koja potpuno odredjuje krivu
    #prolazimo kroz prvih 100 elemenata liste Thi koji odredjuju F_3: u praksi, ja pustim da se prolazi kroz 2000 elemenata liste Thi, 
    #zapamtim rezultate koje sam dobio i pokrenem novi kod sa naredne 2000
    upper_limit = min(UPPER_BOUND, len(Sec))
    for i in range(LOWER_BOUND, upper_limit):
        Q = Sec[i]
        g = Q[0]*X ** 2 + Q[1]*X*Y + Q[2]*X*Z + Q[3]*X*T + Q[4]*X*U + Q[5]*Y ** 2 + Q[6]*Y*Z + Q[7]*Y*T + Q[8]*Y*U + Q[9]*Z ** 2 + Q[10]*Z*T + Q[11]*Z*U + Q[12]*T ** 2 + Q[13]*T*U + Q[14]*U ** 2
        start = time.time()

        if PARALLEL:
            with Pool() as pool:
                pool.map(iteration, [(g, i, j) for j in range(len(R_i_cands[i]))])
        else:
            for j in range(len(R_i_cands[i])):
                iteration((g, i, j))
        print("[PROGRESS] Done {}%".format(100 * (i+1) / len(Sec)))

        end = time.time()
        print("[TIMER] Last iteration (for i={}) done in {}".format(i, end - start))
    print("[INFO] Done")
