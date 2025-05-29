from sage.all import *
    
    
def floor(x):
    return(int(math.floor(x)))


# number of compositions of 't' in 'a' parts of size in [1..'m']
def comp(t,a,m):  
    if t==0 and a==0:
        return(1)
    elif m < 1:
        return(0)
    else:
        return(sum(((-1)**i)*binomial(a,i)*binomial(t-1-m*i, a-1) for i in range(min(a,floor((t-a)/m))+1)))

    
#q odd: sphere size of radius 't' in Z_q^n
def sphere_odd(q,n,t):  
    m = floor((q-1)/2) 
    return(
        sum(
            sum(
                binomial(n+1, j)*binomial(n+1-j, k)*comp(t,j,m)*comp(t,k,m)
                for k in range(n+2-j)
            ) 
            for j in range(n+2)
        ) 
        + 
        2*sum(
            sum(
                sum(
                    binomial(n+1, a)*binomial(n+1-a, b)
                    *(1/sum(binomial(a,i)*binomial(b,i) for i in range(min(a,b)+1)))
                    * sum(
                        sum(
                            binomial(n+1-a-b, j)*binomial(n+1-a-b-j, k)
                            *comp(t-a*(mu),j,mu-1)*comp(t-b*(q-mu),k,q-mu-1)
                            for k in range(n+2-a-b-j)
                        ) 
                        for j in range(n+2-a-b)
                    ) 
                    for b in range(min(n+1-a,floor(t/(q-mu)))+1)
                )
                for a in range(1,min(n+1,floor(t/mu))+1)
            )
            for mu in range(m+1,q)
        )
    )

 

#q even: sphere size of radius 't' in Z_q^n
def sphere_even(q,n,t):  
    m = floor((q-1)/2) # = q/2 - 1,  q = 2m+2
    return(
        sum(
            sum(
                binomial(n+1, j)*binomial(n+1-j, k)*comp(t,j,m)*comp(t,k,m)
                for k in range(n+2-j)
            ) 
            for j in range(n+2)
        ) 
        + 
        2*sum(
            sum(
                binomial(n+1, a)
                *sum(
                    sum(
                        binomial(n+1-a, j)*binomial(n+1-a-j, k)
                        *comp(t-a*(mu),j,mu-1)*comp(t,k,q-mu-1)
                        for k in range(n+2-a-j)
                    ) 
                    for j in range(n+2-a)
                )
                for a in range(1,min(n+1,floor(t/mu))+1)
            )
            for mu in range(m+1,q)
        )
        + 
        sum(
            (1+int(mu > m+1))*sum(  # (1+0) if mu=m+1, (1+1) if mu>m+1
                sum(
                    binomial(n+1, a)*binomial(n+1-a, b)
                    *(1/sum(binomial(a,i)*binomial(b,i) for i in range(min(a,b)+1)))
                    * sum(
                        sum(
                            binomial(n+1-a-b, j)*binomial(n+1-a-b-j, k)
                            *comp(t-a*(mu),j,mu-1)*comp(t-b*(q-mu),k,q-mu-1)
                            for k in range(n+2-a-b-j)
                        ) 
                        for j in range(n+2-a-b)
                    ) 
                    for b in range(1,min(n+1-a,floor(t/(q-mu)))+1)
                )
                for a in range(1,min(n+1,floor(t/mu))+1)
            )
            for mu in range(m+1,q)
        )
    )


#return the weight distribution, as predicted by formulae sphere_even and sphere_odd 
def weight_distr_form(q,n):
    sphere_list = []
    
    t = 0
    while True:
        if q % 2 == 0:
            sphere_size = sphere_even(q,n,t)
        if q % 2 == 1:
            sphere_size = sphere_odd(q,n,t)
        
        if sphere_size <= 0:
            break
        
        sphere_list += [int(sphere_size)]
        t += 1

    return(sphere_list)


#return ball size of radius 't' in Z_q^n
def ball_size(q,n,t):
    if q % 2 == 0:
        s = int(sum([sphere_even(q,n,i) for i in range(t+1)]))
    if q % 2 == 1:
        s = int(sum([sphere_odd(q,n,i) for i in range(t+1)]))
    return(s)


#return upper sphere-packing bound on size of code of minimum distance d
def sphere_pack_bound(q,n,d):
    t = floor((d-1)/2)
    return((q**n)/ball_size(q,n,t))
