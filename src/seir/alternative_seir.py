"""
    @author : nono & matteo 
    @email : nonosaha@mis.mpg.de/matteo@mis.mpg.de

    TODO

"""
import numpy as np
import argparse
import pandas as pd 
import multiprocess as mp 
import analytic  as ac
import matplotlib.pyplot as plt 
import json
from scipy.special import lambertw,factorial


def init(pop_size, int_nb_of_infection, infection_types, N_0): 
    hosts = []
    for i in range(pop_size-int_nb_of_infection) : 
        hosts.append([i,'S',None,[],0,0])
    init_infectious = []
    j = 0
    for i in range(pop_size-int_nb_of_infection,pop_size): 
        if j <int(int_nb_of_infection/2.) : 
            inf_type = 'F'
        else : 
            inf_type = 'R'
        j = j +1
        init_infectious.append([i,"I",inf_type,[1.+infection_types[inf_type]['sb']]*N_0,0,0])
    hosts = hosts + init_infectious
    
    
    return pd.DataFrame(np.array(hosts), columns=["id", "status", "viralParams","pop","tau", "phase"], index=None)

def reproduce(pop, infection_params, T,  C=100) : 
    
    s_b = float(infection_params['sb'])
    s_d = float(infection_params['sd'])
    mu  = float(infection_params['mu'])
    
    pop_0 = np.copy(pop)
    
    f = lambda i, s_d: ((1-s_d)**i)
    t = 0
    
    while 0<len(pop_0) and t<T: 
        pop_n = np.array(pop_0)*f(np.random.poisson(mu,len(pop_0)),s_d)
        survivors = np.random.binomial(1,pop_n/2.,len(pop_n))
        pop_0 = list(pop_n[np.where(survivors==1)[0]])*2
        t = t + 1 
              
    return pop_0

def reproduce2(pop, infection_params, T,  C=100) : 
    
    s_b = float(infection_params['sb'])
    s_d = float(infection_params['sd'])
    mu  = float(infection_params['mu'])
    
    pop_0 = np.copy(pop)
    
    f = lambda i, s_d: ((1-s_d)**i)
    t = 0
    while 0<len(pop_0) and t<T: 
        ks = np.random.poisson(pop_0, size=len(pop_0))
        pop_n = np.array(pop_0)*f(np.random.poisson(mu,len(pop_0)),s_d)
        pop = []
        for i in range(len(ks)) : 
            pop.append([pop_n[i]]*ks[i])

        pop_0 = list(np.concatenate(pop))
        t = t + 1 
              
    return pop_0


def expose(hostPop, params, phase) : 
    
    susceptibles = hostPop[hostPop['status']=='S']['id']
    number_to_infect  = int(np.ceil(params['infectionRate']*len(susceptibles)))
    
    if number_to_infect > 0 :
        
        host_to_infect_index = np.random.choice(susceptibles, size=number_to_infect, replace=False)

        infectious = hostPop[hostPop['status']=='I']['id']
        if len(infectious) > 0 :


            infector_index = np.random.choice(infectious, size=len(host_to_infect_index), replace=True)
        
            vparams = [] 
            pops_0 = []
            statuses = []
            phases = []

            for id_ in list(infector_index) : 
                p_n = list(hostPop[hostPop['id']==id_]['pop'].values[0])
                
                if np.random.random() < len(p_n)/float(params['N_max']) : 
                    vparams.append(hostPop[hostPop['id']==id_]['viralParams'].values[0])  
                    pop_0 = np.random.choice(p_n, size=params['bottleneck_size'])
                    pops_0.append(pop_0)
                    statuses.append('E')
                    phases.append(phase)
                else : 
                    statuses.append('S')
                    pops_0.append([])
                    vparams.append(None)
                    phases.append(0)
              
        
            
            hostPop.update(pd.DataFrame({'status': statuses, 
                                    'viralParams': vparams, 
                                    'pop': pops_0,
                                    'phase': phases}, index=host_to_infect_index))
    
    return hostPop


def infect(hostPop, params, infect_type) : 
    
    exposed = hostPop[hostPop['status']=='E'] 
    infectious = hostPop[hostPop['status']=='I']  
    id_s = exposed["id"].values
    status = []
    pops_0 = []
    taus = []
    viralParams = [] 

    for host in exposed.values  : 

        if len(host[3])> 0 : 
            pop = reproduce2(host[3],infect_type[str(host[2])], params['tao'], params['C'] )
            vp = host[2]
            if len(pop) >= params['C'] : 
                status.append('I')
                taus.append(host[-1]+1)
                viralParams.append(host[2])
               
            if len(pop)==0 : 
                status.append('S')
                taus.append(0)
                viralParams.append("")
            if 0<len(pop)<params['C']: 
                taus.append(host[-1]+params['tao'])
                status.append('E')
                viralParams.append(host[2])
            
            pops_0.append(pop)
        else : 
           
            pops_0.append(host[3])
            status.append('E')
            taus.append(host[-1])
            viralParams.append(host[2])
    
    hostPop.update(pd.DataFrame({'status': status, 
                                 'pop': pops_0,
                                 'tau': taus,
                                 'viralParams': viralParams}, index=id_s))
    pops_n = []
    id_s1 = []
    for host in infectious.values  : 

        if len(host[3])<params['N_max'] : 
            pop = reproduce2(host[3],infect_type[str(host[2])], params['tao'], params['C'] )
            pops_n.append(pop)
            id_s1.append(host[0])
            
        
    hostPop.update(pd.DataFrame({'pop': pops_n}, index=id_s1 ))

    return hostPop


def recover(hostPop, params) : 
    infected = list(hostPop[hostPop['status']=='I']['id'])
    number_to_recover = int(np.ceil(params['recoverRate']*len(infected)))
    
    if number_to_recover >0 : 
        host_to_recovered = np.random.choice(infected, size=number_to_recover, replace=False)
        
        hostPop.update(pd.DataFrame({'status':['R']*len(host_to_recovered)}, index=host_to_recovered))
    
    return hostPop

def epidemic(host_pop, params,infection_type, T) : 
    
    data = []
    infection_number = []
    infected = host_pop[host_pop['status']=='I']
    
    infection_number.append([len(infected[infected['viralParams']=='R']),len(infected[infected['viralParams']=='F'])])
    
    data.append([len(host_pop[host_pop['status']=='S']),len(host_pop[host_pop['status']=='E']), 
        len(host_pop[host_pop['status']=='I']),len(host_pop[host_pop['status']=='R']), 
        len(host_pop[host_pop['viralParams']=='R']),len(host_pop[host_pop['viralParams']=='F'])])
   
    for i in range(T) : 
        print('phase....',i)
        expose(host_pop, params,i)
        E = len(host_pop[host_pop['status']=='E'])
        infect(host_pop, params, infection_type)
        I = host_pop[host_pop['status']=='I']
        recover(host_pop, params)

        data.append([len(host_pop[host_pop['status']=='S']),E, 
           len(I),len(host_pop[host_pop['status']=='R']), 
           len(host_pop[host_pop['viralParams']=='R']),len(host_pop[host_pop['viralParams']=='F'])])
        
        infected = host_pop[host_pop['status']=='I']
    
        infection_number.append([len(I[I['viralParams']=='R']),len(I[I['viralParams']=='F'])])
 
    return data, infection_number, host_pop



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
    parser.add_argument('-t', type=int, default=10, help="Incubation time")
    args = parser.parse_args()
    C = args.C
    T = args.T 
    N = args.N
    mu = args.mu
    lf = args.lf
    tao = args.t
    params = {
        'recoverRate' : args.Rr, 
        'infectionRate': args.Ir, 
        'bottleneck_size': 3, 
        'C': C, 
        'tao': tao, 
        'N_max': 1000
    }

    infection_types = {

    'R' : 
        {
            'sb' : 0.5,
            'sd' : 0.05,
            'mu' : args.mu
        },

    'F':
        {
            'sb' : 0.5,
            'sd' : 0.9,
            'mu' : args.mu
        }
        }
    
    init_infected = 10
    N_0=5

    hosts = init(N,init_infected,infection_types, N_0)
    
    #print(hosts)

    epidemic_data, infections, immune_hosts = epidemic(hosts, params,infection_types, T)  

    
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
    #with open('../logs/seir/Ir/'+str(args.f)+'/seir_data'+str(lf)+'.json', 'w') as fp : 
    #   json.dump(dict_data, fp)
    #immune_hosts.to_csv('../logs/seir/Ir/'+str(args.f)+'/final_host'+str(lf)+'.csv')

    #pd.DataFrame(np.array(infections)).to_csv('../logs/seir/Ir/'+str(args.f)+'/infection_nb'+str(lf)+'.csv')
    for key in ['S', 'E','I', 'R'] : 
        plt.plot(dict_data[key], label=key)
    
    plt.xlabel(r'Time($t$)')
    plt.legend()
    #plt.savefig('../logs/seir/Ir/'+str(args.f)+'/seir'+str(lf)+'.pdf') 
    plt.show()
    infect_data = np.array(infections)
    plt.plot(infect_data[:,0], label="Robust")
    plt.plot(infect_data[:,1], label="Fragile")
    plt.legend()
    plt.show()
    plt.plot(dict_data['robust'], label="Robust")
    plt.plot(dict_data['fragile'], label="Fragile")
    plt.legend()
    plt.show()

if __name__ == "__main__" : 
    main()