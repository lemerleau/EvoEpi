import analytic as ac
import time
import numpy as np


def main():

    viralParams = {
            'sb': 0.5,
            'sd': 0.9,
            'mu': 0.1
            }
    bottleneck_size = 3
    result = []
    sb = 0.5
    sd = 0.9
    u = 0.3

    data = {}
    for u in np.arange(0.05,0.45, 0.05) :
        mus = {}
        print(u)
        for N_0 in range(1,20,1):
            #print("compute for sd =", sd, 'K = ',K(sb, sd, u), 'N_0 ='+str(N_0))
            m = ac.bottleneckExtinctionProb(sb,sd,u,N_0,5)
            mus[N_0] = m
        data[u] = mus

    print(data)

if __name__ == "__main__" :

   main()
