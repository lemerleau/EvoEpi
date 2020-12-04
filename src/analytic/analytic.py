


"""
    @author : nono & matteo
    @email : nonosaha@mis.mpg.de/matteo@mis.mpg.de

    TODO

"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import json
from scipy.special import lambertw,factorial, comb
from scipy.stats import multinomial
import constraint
from numpy.linalg import matrix_power
import multiprocess as mp
import time


def m_(s_b, s_d, u,k) :
    return np.exp(-u)*(1+s_b)*((1-s_d)**k)

def f_i(sd, u, i) :
    return (np.exp(-u/sd)*((u/sd)**(i))/factorial(i, exact=True))


def q_(sb, sd, u, k, i) :

    mk = m_(sb, sd,u,k)
    q_kk = (-1./mk)*lambertw(-mk*np.exp(-mk)).real

    if i> k :
        return 1.
    if k == i :
        return q_kk
    else :
        m_i = m_(sb, sd, u, i)
        q_ki = (-1/m_i)*lambertw(-m_i*np.exp(-m_i*(1+sum([((u**(j-i))/factorial((j-i),exact=True))*(1-q_(sb,sd,u,k,int(j)))
            for j in np.arange(i+1,k+1,step=1)])))).real
        return q_ki
def K(s_b, s_d, u):
    if np.exp(-u)*(1+s_b) < 1.:
        return;
    if np.floor((u-np.log(1+s_b))/np.log(1-s_d))==(u-np.log(1+s_b))/np.log(1-s_b) :
        return int(np.floor((u-np.log(1+s_b))/np.log(1-s_d)) -1)
    else :
        return int(np.floor((u-np.log(1+s_b))/np.log(1-s_d)))

def partialExtinctionProb(s_b, s_d, u, k, n) :

    p_ext = 1.

    for j in range(k+1) :
        p_ext = p_ext*(q_(s_b, s_d,u,k,j)**n[j])

    return p_ext

def clickProba(s_b, s_d, u, k,n) :
    if k == 0 :
        return 1. - partialExtinctionProb(s_b, s_d, u,0,n )
    else :
        return partialExtinctionProb(s_b, s_d, u,k-1,n) - partialExtinctionProb(s_b, s_d, u,k,n )

def extinctionProb(s_b, s_d, u, n):
    return partialExtinctionProb(s_b, s_d, u,K(s_b, s_d, u),n)

def survivalProb(s_b, s_d, u, n) :
    return 1- partialExtinctionProb(s_b, s_d, u,K(s_b, s_d, u),n)


def f(s_b, s_d, u, k, i) :


    sum_ = 0.
    for j in range(K(s_b,s_d,u)+1) :
        if j < k :
            continue;
        else:
            sum_ = sum_ + (np.exp(-u/s_d)*((u/s_d)**(j - k)))/factorial(j-k, exact=True)
    if i< k :
        return 0.
    else :

        f_value = (np.exp(-u/s_d)*((u/s_d)**(i - k))/factorial(i-k, exact=True))/sum_
        return f_value

def Q_k (s_b, s_d, u,B, n, k) :
    kk = K(s_b, s_d, u)+1

    if k > 0 :
        f_s = [0]*k
        for i in range(kk-k) :
            f_s.append(f_i(s_d, u,i))
        f_s.append(1-sum(f_s))
    else :
        f_s = []
        for i in range(kk-k) :
            f_s.append(f_i(s_d, u,i))
        f_s.append(1-sum(f_s))

    rvs = multinomial(B,f_s)

    n_v = list(np.copy(n))
    n_v.append(B-sum(n))

    pmf = rvs.pmf(n_v)

    return pmf


def P_mm(s_b, s_d, u, B, m, n) :
    kk = K(s_b, s_d, u) +1
    zeros = list(np.zeros(kk))

    if n==zeros and m==zeros :
        return 1.
    if n!=zeros and m==zeros :
        return 0.
    if n==zeros and m!=zeros :
        sum_ = 0.
        for k in range(kk) :
            sum_ = sum_ +clickProba(s_b, s_d, u, k,m)*Q_k(s_b,s_d,u,B,n,k)
        return sum_+ extinctionProb(s_b, s_d, u, m)
    if n!=zeros and m!= zeros :
        sum_ = 0.
        for k in range(kk) :
            sum_ = sum_ +clickProba(s_b, s_d, u, k,m)*Q_k(s_b,s_d,u,B,n,k)
        return sum_

def stateSpace(N,K) :

    problem = constraint.Problem()

    problem.addVariables(range(K+1), range(N+1))

    problem.addConstraint(lambda *vect: sum(vect)<=N, range(K+1))
    solutions = problem.getSolutions()

    states = dict(zip(range(len(solutions)),sorted([list(solutions[i].values()) for i in range(len(solutions))])))

    return states


def bottleneckTransitionMat(s_b, s_d, u, B ) :
    kk = K(s_b, s_d, u)
    states = stateSpace(B,kk)
    nb_state = len(states)
    matrix = np.zeros((nb_state, nb_state))

    for i in range(nb_state) :
        for j in range(nb_state) :
            matrix[i][j]= P_mm(s_b,s_d,u,B,states[i],states[j])

    return matrix

def bottleneckExtinctionProb(s_b, s_d, u, B, nMax=10):
    M = bottleneckTransitionMat(s_b, s_d, u, B)
    M_ext = matrix_power(M, nMax)
    return M_ext[-1,0]

def pij(i,j,trans_prob) :
    p_ij = 0.

    for k in range(trans_prob.shape[0]) :
        try :
            p_ij = p_ij + trans_prob[i][k] * trans_prob[k][j]
        except IndexError :
            print(i,k,j,trans_prob.shape[0])
    return p_ij


def asymptoticmean_fitness(sb, sd, u, n) :
    return np.exp(-u)*(1.+sb)*sum([((1.-sd)**k)*clickProba(sb,sd,u,k,n) for k in range(0,K(sb,sd,u)+1)])

def conditasymptoticmean_fitness(sb, sd, u, n) :
    return asymptoticmean_fitness(sb, sd, u, n)/survivalProb(sb, sd, u, n)


def getInfectionNum(params) :

    s_b = params[0]['sb']
    s_d = params[0]['sd']
    u = params[0]['mu']
    B = params[1]
    k = params[3]

    kk = K(s_b, s_d, u)
    states = stateSpace(B,kk)
    nb_state = len(states)
    Probs = []
    for i in range(nb_state) :
        Probs.append(Q_k(s_b, s_d, u, B, states[i], k))

    return np.random.choice(range(nb_state), p=np.array(Probs)/sum(Probs) )



def transmit(params) :
    viralParams = params[0]
    bottleneck_size = params[1]
    x = params[2]

    trans_M = bottleneckTransitionMat(viralParams['sb'],viralParams['sd'],viralParams['mu'],bottleneck_size)

    probs = matrix_power(trans_M, 2)[x]

    return np.random.choice(range(trans_M.shape[0]), 1, p=probs)[0]


def incubation_time(sb,sd,u,B,infection_nb,c=400) :
    kk = K(sb, sd,u)

    states = stateSpace(B,kk)
    n = states[infection_nb]

    for i in range(kk+1) :
        if n[i] != 0 :
            k_0 = i
            break;

    p_k = [clickProba(sb,sd,u,k,n) for k in range(k_0,kk+1)]+[extinctionProb(sb, sd,u,n)]


    k = np.random.choice(list(range(k_0,kk+1))+[100], p=np.array(p_k)/sum(p_k))
    mk = m_(sb, sd,u,k)
    qk = (-1./mk)*lambertw(-mk*np.exp(-mk)).real

    if k == 100:
        return 0,k

    if k == k_0 :
        n_t = n[k]
        prob_l = []
        for l in range(1,n_t+1,1) :
            prob_l.append((comb(n_t, l, exact=False)*((1-qk)**l)*(qk**(n_t-l)))/(1-qk**n_t))

        l = np.random.choice(range(1,n_t+1), 1,p=prob_l)[0]

        Wk = np.random.gamma(shape=l, scale=1/(1-qk))

        if Wk > c*np.exp(-u/sd):
            return 0,100

        return np.ceil((np.log(c)-np.log(Wk)-(u/sd))/np.log(mk)),k
    elif k>k_0 :
        sum_ = 0
        for i in range(k_0, k+1) :
            sum_ += ((n[i]*m_(sb, sd, u,i)*(u**(k-i)))/factorial(k-i, exact=True))*(1+q_(sb,sd,u,k,i)*(
                (1-q_(sb,sd,u,k,k))/(q_(sb,sd,u,k-1,i)-q_(sb,sd,u,k,i))))

        n_t = int(np.ceil(sum_))
        prob_l = []
        for l in range(1,n_t+1,1) :
            prob_l.append((comb(n_t, l, exact=False)*((1-qk)**l)*(qk**(n_t-l)))/(1-qk**n_t))

        l = np.random.choice(range(1,n_t+1), 1,p=prob_l)[0]

        Wk = np.random.gamma(shape=l, scale=1/(1-qk))
        if Wk > c*np.exp(-u/sd):
            return 0,100

        return np.ceil((np.log(c)-np.log(Wk)-(u/sd))/np.log(mk))+1,k


def incubation_time2(sb,sd,u,n,c=400) :
    kk = K(sb, sd,u)

    for i in range(kk+1) :
        if n[i] != 0 :
            k_0 = i
            break;

    p_k = [clickProba(sb,sd,u,k,n) for k in range(k_0,kk+1)]+[extinctionProb(sb, sd,u,n)]


    k = np.random.choice(list(range(k_0,kk+1))+[100], p=np.array(p_k)/sum(p_k))
    mk = m_(sb, sd,u,k)
    qk = (-1./mk)*lambertw(-mk*np.exp(-mk)).real

    #print(k)
    if k == 100 :
        return 0, k

    if k == k_0 :
        n_t = n[k]
        prob_l = []
        for l in range(1,n_t+1,1) :
            prob_l.append((comb(n_t, l, exact=False)*((1-qk)**l)*(qk**(n_t-l)))/(1-qk**n_t))

        points = np.random.binomial(1,1-qk,n_t)
        l = np.random.choice(range(1,n_t+1), 1,p=prob_l)[0]

        Wk = np.random.gamma(shape=l, scale=1/(1-qk))

        return np.ceil((np.log(c)-np.log(Wk)-(u/sd))/np.log(mk)),k

    if k>k_0 :
        sum_ = 0
        for i in range(k_0, k+1) :
            #print(sum_,q(sb,sd,u,k,i),q(sb,sd,u,k-1,i),q(sb,sd,u,k,i), i, k)
            sum_ += ((n[i]*m_(sb, sd, u,i)*(u**(k-i)))/factorial(k-i, exact=True))*(1+q_(sb,sd,u,k,i)*(
                (1-q_(sb,sd,u,k,k))/(q_(sb,sd,u,k-1,i)-q_(sb,sd,u,k,i))))

        n_t = sum_
        Wk = np.random.gamma(shape=n_t*(1-qk)/(1-qk**n_t), scale=1/(1-qk))
        return np.ceil((np.log(c)-np.log(Wk)-(u/sd))/np.log(mk))+1, k
