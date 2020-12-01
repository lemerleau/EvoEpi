"""
    @author : nono & matteo 
    @email : nonosaha@mis.mpg.de/matteo@mis.mpg.de

    TODO

"""

#import necessary libraries 
import numpy 
import pp 
import json
import matplotlib.pyplot as plt 
import pandas
import time

import argparse

import branching as bp

def main() : 

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, argument_default=argparse.SUPPRESS)
    parser.add_argument('-C', type=int, default=1000, help='Carrying capacity.')
    parser.add_argument('-sb', type=float,default=0.5, help="Selective advantage.")
    parser.add_argument('-sd', type=float, default=0.05, help="Deleterious effect of mutation. e.g 0.05 for the robust type and 0.9 for the fragile type.")
    parser.add_argument('-u', type=float, default=0.1, help="Mutation rate. ")
    parser.add_argument('-T', type=int,default=100, help="Maximum number of generations.")
    parser.add_argument('-N_0', type=int,default=5, help="Bottleneck size size.")
    parser.add_argument('-job', type=int,default=10, help="Number of runs.")
    parser.add_argument('-mBN', type=int,default=10, help="Maximum number of bottlenecks.")
    args = parser.parse_args()
    C = args.C
    T = args.T 
    N_0 = args.N_0
    u = args.u
    sb = args.sb
    sd = args.sd
    max_bn = args.mBN 


    number_of_run = args.job

    print ("Starting the BP: ",number_of_run, "job(s)....")
    print ("--------------------------------------")
    print ("Settings : ")
    print ("--------------------------------------")

    print ("N_0 = ", N_0 )
    print ("number of generation = ", T)
    print ("Mutation rate:", u)
    print ("sb = ", sb)
    print ("sd = ", sd)
    print ("--------------------------------------")

    
    pop = numpy.ones(N_0)*(1+sb)

    ppservers = ()


    job_server = pp.Server(10, ppservers=ppservers)
    
    jobs = [(i , job_server.submit(bp.braching_process, (pop, sd, sb, u,N_0,T,C,max_bn), (),( "numpy", "json", "branching"))) for i in range(number_of_run)]
    
    
    data = {
        "Fitness" : [],
        "wtFreq": [],
        "Pop_sizes": [], 
        'number_of_bn': []
    }
    print ("Data will be saved in a json file ../../data/bp/data_robust.json structured as follows:\n", str(data))
    tic = time.time()
    for i,job in jobs : 
        fitnesses,n_0s,Ns,last_pop,nb_bn = job()
        data['Fitness'].append(fitnesses)
        data['wtFreq'].append(n_0s)
        data['Pop_sizes'].append(Ns)
        data['number_of_bn'].append(nb_bn)
    toc = time.time()
    
    if sd < 0.9 : 
        with open("../../data/bp/mu/"+str(numpy.round(u,2))+"/data_robust" +str(N_0)+".json", 'w') as fp : 
            json.dump(data, fp)
    else : 
        with open("../../data/bp/mu/"+str(numpy.round(u,2))+"/data_fragile"+str(N_0)+".json", 'w') as fp : 
            json.dump(data, fp)
    
    print('Branching Process simulation done. in ', toc-tic, 's')

    
if __name__=='__main__': 
    main()
         
