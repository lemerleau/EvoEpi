"""
    @author : nono & matteo 
    @email : nonosaha@mis.mpg.de/matteo@mis.mpg.de

    TODO

"""
import numpy as np
import argparse
import pandas as pd 
import multiprocess as mp 
import ../analytic/analytic  as ac
import matplotlib.pyplot as plt 
import json
from scipy.special import lambertw,factorial


def init(pop_size, int_nb_of_infection, infection_types, params): 
    hosts = []
    for i in range(pop_size-int_nb_of_infection) : 
        hosts.append([i,'S',0,None,0,0])
    init_infectious = []
    j = 0
    for i in range(pop_size-int_nb_of_infection,pop_size): 
        if j <int(int_nb_of_infection/2.) : 
            inf_type = infection_types[0]
        else : 
            inf_type = infection_types[1]
        j = j +1
        init_infectious.append([i,"I",max(params[str(inf_type['sd'])].keys()),inf_type,0,0])
    hosts = hosts + init_infectious
    
    
    return pd.DataFrame(np.array(hosts), columns=["id", "status","infection_nb", "viralParams","incubation_time", "k"], index=None)


def expose(hostPop, params, C) : 
    
    susceptibles = hostPop[hostPop['status']=='S']['id']
    number_to_infect  = int(np.ceil(params['infectionRate']*len(susceptibles)))
    incubation_times = []
    incubation_times_rf = []
    if number_to_infect > 0 :
        
        host_to_infect_index = np.random.choice(susceptibles, size=number_to_infect, replace=False)

        infectious = hostPop[hostPop['status']=='I']['id']
        if len(infectious) > 0 :

            infector_index = np.random.choice(infectious, size=len(host_to_infect_index), replace=True)
        
            vparams = []
            args = []
            
            ks = []
            statuses = []
            infection_nbs = []
            for id_ in list(infector_index) : 
                v_i = hostPop[hostPop['id']==id_]['viralParams'].values[0] 
                arg = (v_i, params['bottleneck_size'],
                        int(hostPop[hostPop['id']==id_]['infection_nb'].values[0] ), int(hostPop[hostPop['id']==id_]['k'].values[0] ))         
                tau,k = ac.incubation_time(v_i['sb'],v_i['sd'],v_i['mu'],params, arg[2], C)
                incubation_times.append(tau)
                incubation_times_rf.append([v_i,tau,arg[2],k])
                #if v_i == {'mu': 0.1, 'sb': 0.5, 'sd': 0.9} : 
                #    print("tao ==== ", v_i, tau, 'n=', n)
                if k ==100 : 
                    statuses.append('S')
                    vparams.append(None)

                else : 
                    statuses.append('E')
                    vparams.append(v_i)
                ks.append(k)
                n = ac.getInfectionNum(arg, params[str(v_i['sd'])])

                infection_nbs.append(n)
        
            
            hostPop.update(pd.DataFrame({'status': statuses, 
                                    'viralParams': vparams, 
                                    'infection_nb': infection_nbs,
                                    'incubation_time': incubation_times,
                                    'k': ks}, index=host_to_infect_index))
     
    return hostPop, incubation_times_rf


def infect(hostPop, params) : 
    
    exposed = hostPop[hostPop['status']=='E']
    id_s = exposed['id'].values 
    updated_incubation_time = np.array(exposed['incubation_time'].values, dtype=float) - 1
    status = []
    for t in updated_incubation_time : 
        if t == 0 : 
            status.append('I')
        else : 
            status.append('E')

        
    hostPop.update(pd.DataFrame({'status': status, 
                                 'incubation_time': updated_incubation_time}, index=id_s))
    
    return hostPop


def recover(hostPop, params) : 
    infected = list(hostPop[hostPop['status']=='I']['id'])
    number_to_recover = int(np.ceil(params['recoverRate']*len(infected)))
    
    if number_to_recover >0 : 
        host_to_recovered = np.random.choice(infected, size=number_to_recover, replace=False)
        
        hostPop.update(pd.DataFrame({'status':['R']*len(host_to_recovered)}, index=host_to_recovered))
    
    return hostPop

def epidemic(host_pop, params,infection_type, T, C) : 
    
    data = []
    infection_number = []
    infected = host_pop[host_pop['status']=='I']
    incubation_times_data = [] 
    
    infection_number.append([len(infected[infected['viralParams']==infection_type[0]]),len(infected[infected['viralParams']==infection_type[1]])])
    
    data.append([len(host_pop[host_pop['status']=='S']),len(host_pop[host_pop['status']=='E']), 
        len(host_pop[host_pop['status']=='I']),len(host_pop[host_pop['status']=='R']), 
        len(host_pop[host_pop['viralParams']==infection_type[0]]),len(host_pop[host_pop['viralParams']==infection_type[1]])])
   
    for i in range(T) : 
        print('phase....',i)
        hts,its = expose(host_pop, params,C)
        infect(host_pop, params)
        recover(host_pop, params)
        incubation_times_data.append(its)
        data.append([len(host_pop[host_pop['status']=='S']),len(host_pop[host_pop['status']=='E']), 
           len(host_pop[host_pop['status']=='I']),len(host_pop[host_pop['status']=='R']), 
           len(host_pop[host_pop['viralParams']==infection_type[0]]),len(host_pop[host_pop['viralParams']==infection_type[1]])])
        infected = host_pop[host_pop['status']=='I']

        infection_number.append([len(infected[infected['viralParams']==infection_type[0]]),len(infected[infected['viralParams']==infection_type[1]])])

    return data, infection_number, host_pop, incubation_times_data



def main() : 
    
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, argument_default=argparse.SUPPRESS)
    parser.add_argument('-C', type=int, default=100, help='Carrying capacity')
    parser.add_argument('-Ir', type=float,default=0.1, help="Infection rate")
    parser.add_argument('-Rr', type=float, default=0.1, help="Recovery rate")
    parser.add_argument('-mu', type=float, default=0.1, help="Mutation rate")
    parser.add_argument('-T', type=int,default=100, help="Time")
    parser.add_argument('-N', type=int,default=1000, help="Population size")
    parser.add_argument('-lf', type=float,default=1, help="Log file name") 
    parser.add_argument('-f', type=int,default=1, help="Log folder name") 
    args = parser.parse_args()
    C = args.C
    T = args.T 
    N = args.N
    mu = args.mu
    lf = args.lf
    params = {
        'recoverRate' : args.Rr, 
        'infectionRate': args.Ir, 
        'bottleneck_size': 3
    }

    infection_types = [
        {
            "sb" : 0.5,
            "sd" : 0.05,
            "mu" : args.mu
        },
        {
            "sb" : 0.5,
            "sd" : 0.9,
            "mu" : args.mu
        }
    ]
    
    init_infected = 100

    hosts = init(N,init_infected,infection_types)
    
    #print(hosts.values)

    epidemic_data, infections, immune_hosts, ic_ts = epidemic(hosts, params,infection_types, T, C)  

    
    data = np.array(epidemic_data, dtype=float)

    dict_data = {
        
        'S': list(data[:,0]),
        'E': list(data[:,1]),
        'I': list(data[:,2]),
        'R': list(data[:,3]),
        'robust' : list(data[:,4]),
        'fragile' : list(data[:,5]),
    }
    print (dict_data)
    with open('seir_data3/'+str(args.f)+'/seir_data'+str(lf)+'.json', 'w') as fp : 
        json.dump(dict_data, fp)
    immune_hosts.to_csv('seir_data3/'+str(args.f)+'/final_host'+str(lf)+'.csv')

    pd.DataFrame(np.array(infections)).to_csv('seir_data3/'+str(args.f)+'/infection_nb'+str(lf)+'.csv')
    for key in ['S', 'E','I', 'R'] : 
        plt.plot(dict_data[key], label=key)
    print (ic_ts)
    plt.xlabel(r'Time($t$)')
    plt.legend()
    plt.savefig('seir_data3/'+str(args.f)+'/seir'+str(lf)+'.pdf') 
    plt.show()


if __name__ == "__main__" : 
    main()