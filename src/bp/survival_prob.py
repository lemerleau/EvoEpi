import numpy  
import pandas
import pp 


def K(s_b, s_d, u): 
    if numpy.exp(-u)*(1+s_b) < 1.: 
        return; 
    if numpy.floor((u-numpy.log(1+s_b))/numpy.log(1-s_d))==(u-numpy.log(1+s_b))/numpy.log(1-s_b) : 
        return int(numpy.floor((u-numpy.log(1+s_b))/numpy.log(1-s_d)) -1)
    else : 
        return int(numpy.floor((u-numpy.log(1+s_b))/numpy.log(1-s_d)))
        

def reproduce_with_bottleneck(pop, s_b,s_d, mu, T,B, nb_bn, C=100) : 
    
    
    pop_0 = numpy.copy(pop)
    
    f = lambda i, s_d: ((1-s_d)**i)
    t = 0
    b = 0
    while b<nb_bn : 
        ks = numpy.random.poisson(pop_0, size=len(pop_0))
        pop_n = numpy.array(pop_0)*f(numpy.random.poisson(mu,len(pop_0)),s_d)
        new_pop = []
        for i in range(len(ks)) : 
            new_pop.append([pop_n[i]]*ks[i])
        if len(new_pop) > 0 : 
            pop_0 = list(numpy.concatenate(new_pop))
        else : 
            return new_pop
        t = t + 1 
        
        if len(pop_0) >= C :
            pop_0 = numpy.random.choice(pop_0, size=B)
            b = b+1
              
    return pop_0




def surv_prob(pop, sb, sd, u, B, b,C, k_classes) : 

    surv = 0

    for i in range(1000) : 
        p = reproduce_with_bottleneck(pop, sb,sd, u, 100,B, nb_bn=b, C=C) 
        if(all(x not in p for x in k_classes)):
            surv += 1
    return 1- surv/1000.


def run(pop,sb,u,B,b,C) : 
    probs = []
    for sd in numpy.arange(0.03,0.2,0.01) : 
        kk = K(sb,sd,u)
        print(kk)
        k_classes = []
        for i in range(kk + 1 ) : 
            k_classes.append((1+sb)*(1-sd)**i)
        print(k_classes)
        probs.append(surv_prob(pop,sb,sd,u,B,b,C,k_classes))
    return probs 


def main() : 

    sb = .2 
    u = 0.05 
    B = 5
    C = 500
    pop = [1+sb]*B

    ppservers = ()


    job_server = pp.Server(10, ppservers=ppservers)
    
    jobs = [(i , job_server.submit(run, (pop, sb, u,B,i,C), (),( "numpy",))) for i in range(1,6)]
    
    for i,job in jobs : 
        prob = job()
        print (i, prob)
