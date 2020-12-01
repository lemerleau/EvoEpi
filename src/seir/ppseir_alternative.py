"""
    @author : nono & matteo 
    @email : nonosaha@mis.mpg.de/matteo@mis.mpg.de

    TODO

"""
import numpy 
import argparse
import pandas 
import json
import pp
import alternative_seir

def run (N,T,params,infection_types,init_infected, log_folder, lf): 

    hosts = alternative_seir.init(N,init_infected,infection_types,params['bottleneck_size'])
    
    #print(hosts)

    epidemic_data, infections, immune_hosts = alternative_seir.epidemic(hosts, params,infection_types, T)  
    
    """
    ic_ts_nempty = []
    for ic in ic_ts : 
        ic_ts_nempty += list(ic)
    pandas.DataFrame(ic_ts_nempty).to_csv(log_folder+"/ic_t"+str(lf)+".csv")
    """
    data = numpy.array(epidemic_data, dtype=float)

    dict_data = {
        
        'S': list(data[:,0]),
        'E': list(data[:,1]),
        'I': list(data[:,2]),
        'R': list(data[:,3]),
        'robust' : list(data[:,4]),
        'fragile' : list(data[:,5]),
    }
    
    with open(log_folder+'/seir_data'+str(lf)+'.json', 'w') as fp : 
        json.dump(dict_data, fp)
    immune_hosts.to_csv(log_folder+'/final_host'+str(lf)+'.csv')

    pandas.DataFrame(numpy.array(infections)).to_csv(log_folder+'/infection_nb'+str(lf)+'.csv')
     
    return dict_data

def main() : 
    
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, argument_default=argparse.SUPPRESS)
    parser.add_argument('-C', type=int, default=100, help='Carrying capacity')
    parser.add_argument('-Ir', type=float,default=0.1, help="Infection rate")
    parser.add_argument('-Rr', type=float, default=0.1, help="Recovery rate")
    parser.add_argument('-mu', type=float, default=0.1, help="Mutation rate")
    parser.add_argument('-T', type=int,default=100, help="Time")
    parser.add_argument('-N', type=int,default=1000, help="Population size")
    parser.add_argument('-lf', type=float,default=1., help="Log file name") 
    parser.add_argument('-f', type=int,default=1, help="Log folder name") 
    parser.add_argument('-t', type=int, default=10, help="Incubation time")
    args = parser.parse_args()
    C = args.C
    T = args.T 
    N = args.N
    mu = args.mu
    lf = args.lf
    tao = args.t
    nb_jobs = 50
    parameters = []
    log_folders = [] 

    for job in range(nb_jobs) : 

        params = {
            'recoverRate' : args.Rr, 
            'infectionRate': args.Ir, 
            'bottleneck_size': 5, 
            'C': C, 
            'tao': tao, 
            'N_max': 1000
        }
        parameters.append(params)
        root ='../logs/simulated_seir/Ir/'     
        if args.Ir==0.1 or args.Ir==0.2:
            root = root+str(args.Ir)+'0/'
            if mu == 0.1 or mu==0.3 or mu==0.2 or mu==0.4 :
                root = root+str(mu)+'0'
            else :
                root = root+str(mu)
        else : 
            root = root+str(args.Ir)+'/'
            if mu == 0.1 or mu==0.3 or mu==0.2 or mu==0.4 :
                root = root+str(mu)+'0'
            else :
                root = root+str(mu)
            
    
        log_folders.append(root)

    
    infection_types = {
            'R':
        {
            "sb" : 0.5,
            "sd" : 0.05,
            "mu" : args.mu
        },

        'F': 

        {
            "sb" : 0.5,
            "sd" : 0.9,
            "mu" : args.mu
        }
        }
    
    init_infected = 100


   
    ppservers = ()

    job_server = pp.Server(20, ppservers=ppservers)
    
    jobs = [(i , job_server.submit(run, (N,T,parameters[i],infection_types,init_infected, log_folders[i],i), (),( "numpy", "random","pandas","os", "json","alternative_seir"))) for i in range(nb_jobs)]
    
    results = []
    for i,job in jobs : 
        results.append(job())
    print(results)

if __name__ == "__main__" : 
    main()