"""
    @author : nono & matteo 
    @email : nonosaha@mis.mpg.de/matteo@mis.mpg.de

    TODO

"""

import numpy 
import branching

def braching_process(pop, s_d, s_b, mu,N_0,max_gen,C, max_bt=5) : 
    
    f = lambda i, s_d: ((1-s_d)**i)
   
    pop_0 = numpy.copy(pop)
    wtFreq = []
    wtFreq.append(len(numpy.where(pop_0==1+s_b)[0])/len(pop_0))
    mean_fitness = []
    mean_fitness.append(numpy.mean(pop_0))
    pop_size = []
    t = 0
    nb_bn = 0
    while t<max_gen and nb_bn <max_bt :
        pop_size.append(float(len(pop_0)))
        ks = numpy.random.poisson(pop_0, size=len(pop_0))
        pop_n = numpy.array(pop_0)*f(numpy.random.poisson(mu,len(pop_0)),s_d)
        new_pop = []
        for i in range(len(ks)) : 
            new_pop.append([pop_n[i]]*ks[i])

        pop_0 = list(numpy.concatenate(new_pop))
        
        try : 
            wtFreq.append(len(numpy.where(numpy.array(pop_0)==1+s_b)[0])/float(len(pop_0)))
            mean_fitness.append(numpy.mean(pop_0))
        except ZeroDivisionError : 
            wtFreq.append(0)
            mean_fitness.append(0)
            pop_size.append(len(pop_0))
            break;
        if len(pop_0) >= C:
            pop_size.append(len(pop_0))
            pop_0 = numpy.random.choice(pop_0,N_0)
            nb_bn +=1
        
        t = t + 1 
    
    return mean_fitness,wtFreq,pop_size, pop_0, nb_bn


