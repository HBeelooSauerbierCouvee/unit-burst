#Unit-burst graphs
#by Hugo Beeloo Sauerbier Couvée
#29-05-2025


###### Several functions to compute bursts, spectra, graphs, weight distributions ######

#create list of positive bursts (0...011...10...0)
def burst_list(length, group):
    cur_list = []
    #add bursts to list for all starting and ending indices i,j
    for i in range(length+1):
        for j in range(i+1,length+1):
            cur_list +=  [group([0]*(i) + [1]*(j-i) + [0]*(length-j))]  
    return(cur_list)    

    
#create sorted list/spectrum (largest to smallest) of unique eigenvalues of graph
def sort_spec(graph):
    return(sorted(list(set(graph.spectrum())))[::-1])

    
#create unit-burst graph as Cayley graph on (Z/qZ)^n with bursts as generators
def UBgraph_gen(q,n):
    #define group (Z/qZ)^n
    Zqn = AbelianGroup([q]*n) 
    #Zqn.semigroup_generators = Zqn.group_generators   # might be needed to generate cayley graph

    #define corresponding unit-burst graph
    #by converting directed edges to undirected, we automatically take care of negative bursts (0,...,0,-1,-1,...,-1,0,...,0)
    return((Zqn.cayley_graph(generators=burst_list(n, Zqn))).to_undirected())


#check whether list of weights satisfies necessary conditions for a valid weight distribution 
def check_weights(q,n,weights):
    satisfied = True
    if not all(int(weight)==weight for weight in weights):
        print('weights not integral')
        satisfied = False
    if any(weight < 0 for weight in weights):
        print('weights negative')
        satisfied = False
    if weights[0] != 1:
        print('zero weight not correct')
        satisfied = False
    if sum(weights) != q**n:
        print('sum of weights not correct')
        satisfied = False
    return(satisfied)


#Give the weight distribution, i.e. the sphere sizes centered around 0, with different methods
#
#   optional arguments: 
#      'method=method_name' with method_name from {'formula', 'graph', 'dynamic'}. Default method is 'formula' 
#      'graph=G'   with a specified graph G, being used with method='graph'. Default graph is UBgraph_gen(q,n) 
#   methods: 
#      'formula': fastest and uses theoretic formulas
#      'graph':   experimentally verified but slowest for large parameters 
#      'dynamic': experimentally verified and speed between 'formula' and 'graph'
def weight_distr(q,n,**kwargs):
    method = kwargs.get('method', 'formula')
    print(f"Method '{method}' selected")
    
    #formula method using formulae from Unit_Burst_sphere_sizes package
    if method == 'formula':
        from Unit_Burst_sphere_sizes import weight_distr_form
        weights = weight_distr_form(q,n)
        check_weights(q,n,weights)
        return(weights)
    
    #graph method that uses built-in distances_distribution() function
    elif method == 'graph':
        graph = kwargs.get('graph', UBgraph_gen(q,n))
        weights = [1]+[(graph.order()-1)*val for val in graph.distances_distribution().values()]
        check_weights(q,n,weights)
        return(weights)
    
    #dynamic programming method
    elif method == 'dynamic':
        from bisect import bisect_left
        
        #initialization
        Zqn = AbelianGroup([q]*n) 
        weights = [[Zqn(1)]]
        i = 1        
        bursts = burst_list(n,Zqn)
        if q >= 3:
            bursts = bursts + [burst^(-1) for burst in bursts]
        
        while True:
            #'weights' will be a list of lists, with at each index i a list of weight i vectors
            #when no longer needed for computation, the list of weight i-3 vectors can be replaced by its length
            weights = weights + [[]]
            if i >= 3:
                weights[i-3] = len(weights[i-3])
                
            #define a sliding window 'wt_window', a sorted list consisting of the vectors of weight i-1 and i-2
            wt_window = weights[i-1]
            if i >= 2:
                wt_window = wt_window + weights[i-2]
            wt_window = sorted(wt_window)
            
            #add single bursts to weight i-1 vectors. When a sum is not in the window, it must be of weight i and is added to weights[i].
            for v in weights[i-1]:
                for b in bursts:
                    elem = v*b
                    #as the window is sorted, we can efficiently check if 'elem' is in the window; if not, we insert it respecting the order. 
                    index = bisect_left(wt_window,elem)
                    if index == len(wt_window) or (index != len(wt_window) and wt_window[index] != elem):
                        weights[i] = weights[i]+[elem]
                        wt_window.insert(index,elem)
            
            #when we encounter no more vectors of larger weight, we stop. Alternative halting condition can be when the sum of weights exceeds q^n
            if len(weights[i]) <= 0:
                break
            i += 1
            
        #remove trailing [] at last index in 'weights', and replace remaining lists by their lengths
        weights = weights[:-1]
        weights[-2] = len(weights[-2])
        weights[-1] = len(weights[-1])
        
        check_weights(q,n,weights)
        return(weights)
    

    else: 
        print('unkown method')

        
#return k-independence number of graph, i.e. size of largest set of vertices such that any two vertices in the set are at distance larger than k
def k_independence_number(graph, k):
    k_graph = graph.distance_graph(list(range(1,k+1)))
    return(k_graph.independent_set(algorithm='Cliquer', value_only=True))





###### main commands to compute unit-burst graphs and their properties ######

from Unit_Burst_sphere_sizes import *
import time

#set parameters: q >= 2 modulus for cyclic ring Z/qZ; n >= 2, dimension for (Z/qZ)^n; d >= 1 for computing maximum d-code
q = 3
n = 4    
d = 3

#define group (Z/qZ)^n
Zqn = AbelianGroup([q]*n) 


#define corresponding unit-burst graph
start = time.time()
UBgraph = UBgraph_gen(q,n) 
end = time.time()
print(f'Time to generate graph: {end - start} sec \n')

#print basic properties of graph
print(f'Order: {UBgraph.order()} = q^n = {q}^{n}')
print(f'Regular: '+('yes' if UBgraph.is_regular() else 'no :(')
      +f' with Degree: {UBgraph.average_degree()}'
      +(f' = 2*binom(n+1,2) for q={q}' if q>2 else ' = binom(n+1,2) for q=2') 
     +'\n')

#print(UBgraph.graph6_string()) #for exporting to graph6 format

#print weight distribution, i.e. sphere sizes around 0
start = time.time()
UBweights = weight_distr(q,n,method='formula')
print('Predicted weight distribut:', UBweights)
end = time.time()
print(f'Time: {end - start} sec \n')

start = time.time()
print('Actual weight distribution:', weight_distr(q,n,method='dynamic'))
end = time.time()
print(f'Time: {end - start} sec \n')

#calculate diameter of graph
print(f'Diameter: {UBgraph.diameter()} (graph method) = {len(UBweights)-1} (formula method)\n')

#calculate maximum size d-code, compared to sphere-packing bound
print(f'Maximum size {d}-code or ({d}-1)-independent set:')
time.sleep(1.0) #to ensure output above is printed before computing
print(f'{k_independence_number(UBgraph, d-1)} <= {N(sphere_pack_bound(q,n,d))} (SP bound) \n')


#print spectrum if not too large (usually order < 1024)
print('Spectrum:')
if q**n <= 1024 and (q,n)!=(4,5):   #maximal (q,n) feasible: (32,2), (10,3), (5,4), (3,5), (3,6), (2,7),...,(2,10)  
    if q**n >= 500:
        time.sleep(1.0) #to ensure output above is printed before computing spectrum for large graphs
    print(sort_spec(UBgraph))
else:
    print('order might be too large to compute full spectrum directly')


#plot graph if order not too large
if q**n <= 70:
    show(UBgraph.plot(layout='circular'))

    
