#!/usr/bin/sage

import os
import pickle
import time
from multiprocessing import Manager, Pool

import sage.all
from sage.all import *

#zameni py notaciju


def switch_using_m(P, Q, R, m):
    teri = GF(2)['x,y,z,t,u']
    X,Y,Z,T,U = teri.gens()
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
    hP = gP(X, Y, Z, T, U)
    nP = (hP.coefficient(X**2), hP.coefficient(X*Y), hP.coefficient(X*Z), hP.coefficient(X*T), hP.coefficient(X*U), hP.coefficient(Y**2), hP.coefficient(Y*Z), hP.coefficient(Y*T), hP.coefficient(Y*U), hP.coefficient(Z**2), hP.coefficient(Z*T), hP.coefficient(Z*U), hP.coefficient(T**2), hP.coefficient(T*U), hP.coefficient(U**2))
    hQ = gQ(X, Y, Z, T, U)
    nQ = (hQ.coefficient(X**2), hQ.coefficient(X*Y), hQ.coefficient(X*Z), hQ.coefficient(X*T), hQ.coefficient(X*U), hQ.coefficient(Y**2), hQ.coefficient(Y*Z), hQ.coefficient(Y*T), hQ.coefficient(Y*U), hQ.coefficient(Z**2), hQ.coefficient(Z*T), hQ.coefficient(Z*U), hQ.coefficient(T**2), hQ.coefficient(T*U), hQ.coefficient(U**2))
    hR = gR(X, Y, Z, T, U)
    nR = (hR.coefficient(X**2), hR.coefficient(X*Y), hR.coefficient(X*Z), hR.coefficient(X*T), hR.coefficient(X*U), hR.coefficient(Y**2), hR.coefficient(Y*Z), hR.coefficient(Y*T), hR.coefficient(Y*U), hR.coefficient(Z**2), hR.coefficient(Z*T), hR.coefficient(Z*U), hR.coefficient(T**2), hR.coefficient(T*U), hR.coefficient(U**2))
    el = [nP, nQ, nR]
    return el

def indicator(P, Q, R):
    if P in OP[3] and Q not in OP[3] and R not in OP[3]:
        indikator = 3
        ind = OP[3].index(P)
    elif P not in OP[3] and Q in OP[3] and R not in OP[3]:
        indikator = 3
        ind = OP[3].index(Q)
    elif P not in OP[3] and Q not in OP[3] and R in OP[3]:
        indikator = 3
        ind = OP[3].index(R)
    elif P in OP[4] and Q not in OP[4] and R not in OP[4]:
        indikator = 4
        ind = OP[4].index(P)
    elif P not in OP[4] and Q in OP[4] and R not in OP[4]:
        indikator = 4
        ind = OP[4].index(Q)
    elif P not in OP[4] and Q not in OP[4] and R in OP[4]:
        indikator = 4
        ind = OP[4].index(R)
    elif P in OP[6] and Q not in OP[6] and R not in OP[6]:
        indikator = 6
        ind = OP[6].index(P)
    elif P not in OP[6] and Q in OP[6] and R not in OP[6]:
        indikator = 6
        ind = OP[6].index(Q)
    elif P not in OP[6] and Q not in OP[6] and R in OP[6]:
        indikator = 6
        ind = OP[6].index(R)
    elif P in OP[1] and Q not in OP[1] and R not in OP[1]:
        indikator = 1
        ind = OP[1].index(P)
    elif P not in OP[1] and Q in OP[1] and R not in OP[1]:
        indikator = 1
        ind = OP[1].index(Q)
    elif P not in OP[1] and Q not in OP[1] and R in OP[1]:
        indikator = 1
        ind = OP[1].index(R)
    elif P in OP[3] and Q in OP[3] and R in OP[3]:
        indikator = 333
        ind = -1
    elif P in OP[4] and Q in OP[4] and R in OP[4]:
        indikator = 444
        ind = -1
    elif P in OP[6] and Q in OP[6] and R in OP[6]:
        indikator = 666
        ind = -1
    elif P in OP[1] and Q in OP[1] and R in OP[1]:
        indikator = 111
        ind = -1
    else:
        indikator = -1
        ind = -1
    el = [indikator, ind]
    return el

def do_job_iii(indikator, P1, Q1, R1, description):
    for indeks in [OP[indikator%10].index(P1), OP[indikator%10].index(Q1), OP[indikator%10].index(R1)]:
        do_job_i(indikator, indeks, P1, Q1, R1, description)

manager = Manager()
shared_list = manager.list()

def handleMat(args):
    nP, nQ, nR, raw_indikator, index, length, mat, description = args
    #print("Progress: {} Iteration {}/{}.".format(description, index, length))
    P, Q, R = switch_using_m(nP, nQ, nR, mat)
    fP = P[0]*X**2 + P[1]*X*Y + P[2]*X*Z + P[3]*X*T + P[4]*X*U + P[5]*Y**2 + P[6]*Y*Z + P[7]*Y*T + P[8]*Y*U + P[9]*Z**2 + P[10]*Z*T + P[11]*Z*U + P[12]*T**2 + P[13]*T*U + P[14]*U**2
    fQ = Q[0]*X**2 + Q[1]*X*Y + Q[2]*X*Z + Q[3]*X*T + Q[4]*X*U + Q[5]*Y**2 + Q[6]*Y*Z + Q[7]*Y*T + Q[8]*Y*U + Q[9]*Z**2 + Q[10]*Z*T + Q[11]*Z*U + Q[12]*T**2 + Q[13]*T*U + Q[14]*U**2
    fR = R[0]*X**2 + R[1]*X*Y + R[2]*X*Z + R[3]*X*T + R[4]*X*U + R[5]*Y**2 + R[6]*Y*Z + R[7]*Y*T + R[8]*Y*U + R[9]*Z**2 + R[10]*Z*T + R[11]*Z*U + R[12]*T**2 + R[13]*T*U + R[14]*U**2
    I = Ideal(fP, fQ, fR)
    if raw_indikator < 10:
        if I in ideals[raw_indikator]:
            shared_list.append(ideals[raw_indikator].index(I))
    else:
        for i in range(1, 4):
            if I in ideals[raw_indikator][i]:
                shared_list.append(ideals[raw_indikator][i].index(I))


def do_job_i(raw_indikator, index, P1, Q1, R1, description):
    indikator = raw_indikator %10
    n = matricesP[indikator][index]
    m = n.inverse()
    nP, nQ, nR = switch_using_m(P1, Q1, R1, m)

    shared_list[:] = []
    pool = Pool()
    data = [[nP, nQ, nR, raw_indikator, index, len(Stab[indikator]), Stab[indikator][index], description] for index in range(len(Stab[indikator]))]
    pool.map(handleMat, data)
    pool.close()

    blacklist = list(set(list(shared_list)))
    if blacklist != []:
        if raw_indikator in [1, 3, 4, 6]:
            P_blacklist = [P_list[indikator][e] for e in blacklist]
            ideals_blacklist = [ideals[indikator][e] for e in blacklist]
            P_list[indikator] = [e for e in P_list[indikator] if e not in P_blacklist]
            ideals[indikator] = [e for e in ideals[indikator] if e not in ideals_blacklist]
        elif raw_indikator in [111, 333, 444, 666]:
            for i in range(1, 4):
                P_blacklist = [P_list[raw_indikator][i][e] for e in blacklist]
                ideals_blacklist = [ideals[raw_indikator][i][e] for e in blacklist]
                P_list[raw_indikator][i] = [e for e in P_list[raw_indikator][i] if e not in P_blacklist]
                ideals[raw_indikator][i] = [e for e in ideals[raw_indikator][i] if e not in ideals_blacklist]
        else:
            print(raw_indikator, "fun")
        #if raw_indikator in [1, 3, 4, 6]:      
        #    P_blacklist = [P_list[indikator][e] for e in blacklist]
        #    ideals_blacklist = [ideals[indikator][e] for e in blacklist]
        #    P_list[indikator] = list(set(P_list[indikator]) - set(P_blacklist))
        #    ideals[indikator] = list(set(ideals[indikator]) - set(ideals_blacklist))
        #elif raw_indikator in [111, 333, 444, 666]:
        #    for i in range(1, 4):
        #        P_blacklist = [P_list[raw_indikator][i][e] for e in blacklist]
        #        ideals_blacklist = [ideals[raw_indikator][i][e] for e in blacklist]
        #        P_list[raw_indikator][i] = list(set(P_list[raw_indikator][i]) - set(P_blacklist))
        #        ideals[raw_indikator][i] = list(set(ideals[raw_indikator][i]) - set(ideals_blacklist))
        #else:
        #    print(raw_indikator, "fun")
#main

if __name__ == "__main__":
    try:
        with open(os.path.join('data', 'prepared_final_data.pkl'), 'rb') as f:
               Stab, matricesP, OP, P_list, ideals, all_el, base_c = pickle.load(f)
    except Exception as e:
        print(e)
        print("You need to run `calculate_data.py` first!")
        exit(1)
    print("Started")
    teri = GF(2)['x,y,z,t,u']
    X,Y,Z,T,U = teri.gens()
    bclen = len(base_c)
    brojac = 0
    while all_el != [] and brojac < 10:
        start = time.time()
        brojac = brojac + 1
        P0, Q0, R0 = all_el[0]
        fP0 = P0[0]*X**2 + P0[1]*X*Y + P0[2]*X*Z + P0[3]*X*T + P0[4]*X*U + P0[5]*Y**2 + P0[6]*Y*Z + P0[7]*Y*T + P0[8]*Y*U + P0[9]*Z**2 + P0[10]*Z*T + P0[11]*Z*U + P0[12]*T**2 + P0[13]*T*U + P0[14]*U**2
        fQ0 = Q0[0]*X**2 + Q0[1]*X*Y + Q0[2]*X*Z + Q0[3]*X*T + Q0[4]*X*U + Q0[5]*Y**2 + Q0[6]*Y*Z + Q0[7]*Y*T + Q0[8]*Y*U + Q0[9]*Z**2 + Q0[10]*Z*T + Q0[11]*Z*U + Q0[12]*T**2 + Q0[13]*T*U + Q0[14]*U**2
        fR0 = R0[0]*X**2 + R0[1]*X*Y + R0[2]*X*Z + R0[3]*X*T + R0[4]*X*U + R0[5]*Y**2 + R0[6]*Y*Z + R0[7]*Y*T + R0[8]*Y*U + R0[9]*Z**2 + R0[10]*Z*T + R0[11]*Z*U + R0[12]*T**2 + R0[13]*T*U + R0[14]*U**2
        with open(os.path.join("data", "result.txt"), "a+") as f:
            f.write("[{}, {}, {}], ".format(P0, Q0, R0))
        for k in range(0, bclen):
            m = base_c[k]
            description = "{}%".format(k/bclen * (1 / len(all_el)) * 100)
            print(k)
            h1 = m[0, 0]*fP0 + m[0, 1]*fQ0 + m[0, 2]*fR0
            h2 = m[1, 0]*fP0 + m[1, 1]*fQ0 + m[1, 2]*fR0
            h3 = m[2, 0]*fP0 + m[2, 1]*fQ0 + m[2, 2]*fR0
            P1 = (h1.coefficient(X**2), h1.coefficient(X*Y), h1.coefficient(X*Z), h1.coefficient(X*T), h1.coefficient(X*U), h1.coefficient(Y**2), h1.coefficient(Y*Z), h1.coefficient(Y*T), h1.coefficient(Y*U), h1.coefficient(Z**2), h1.coefficient(Z*T), h1.coefficient(Z*U), h1.coefficient(T**2), h1.coefficient(T*U), h1.coefficient(U**2))
            Q1 = (h2.coefficient(X**2), h2.coefficient(X*Y), h2.coefficient(X*Z), h2.coefficient(X*T), h2.coefficient(X*U), h2.coefficient(Y**2), h2.coefficient(Y*Z), h2.coefficient(Y*T), h2.coefficient(Y*U), h2.coefficient(Z**2), h2.coefficient(Z*T), h2.coefficient(Z*U), h2.coefficient(T**2), h2.coefficient(T*U), h2.coefficient(U**2))
            R1 = (h3.coefficient(X**2), h3.coefficient(X*Y), h3.coefficient(X*Z), h3.coefficient(X*T), h3.coefficient(X*U), h3.coefficient(Y**2), h3.coefficient(Y*Z), h3.coefficient(Y*T), h3.coefficient(Y*U), h3.coefficient(Z**2), h3.coefficient(Z*T), h3.coefficient(Z*U), h3.coefficient(T**2), h3.coefficient(T*U), h3.coefficient(U**2))
            indikator, ind = indicator(P1, Q1, R1)
            if indikator in [1, 3, 4, 6]:
                do_job_i(indikator, ind, P1, Q1, R1, description)
            else:
                do_job_iii(indikator, P1, Q1, R1, description)
        all_el = P_list[1] + P_list[3] + P_list[4] + P_list[6] + P_list[111][1] + P_list[333][1] + P_list[444][1] + P_list[666][1]
        end = time.time()
        print('Operation {} lasts: {} seconds'.format(brojac, end - start))
        print(len(all_el))
        with open(os.path.join('data', 'prepared_final_data.pkl'), 'wb') as f:
            pickle.dump([Stab, matricesP, OP, P_list, ideals, all_el, base_c], f)
    
    print(len(all_el))
    
#needed: lists, Stabs, OPs, matrices 
