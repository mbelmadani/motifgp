from math import log

class MOEATerminationDetection(object):
    """
    __author__ = "Manuel Belmadani"
    __copyright__ = "Copyright 2016, MotifGP"
    __credits__ = ["Manuel Belmadani"]

    __license__ = "LGPL v3.0"
    __version__ = "1.0.0"
    __maintainer__ = "Manuel Belmadani"
    __email__ = "mbelm006@uottawa.ca"
    __status__ = "Development"
    
    """

    """
    This class is a generic implementation of the Entropy-Based Termination Criterion for Multiobjective Evolutionary Algorithms by Saxena et al.

    Publication: 
    Saxena, D. K., Sinha, A., Duro, J. A. & Zhang, Q. Entropy-Based Termination Criterion for Multiobjective Evolutionary Algorithms. IEEE Trans. Evol. Comput. 20, 485-498 (2016).

    The class is initialized with the runtime parameters and requires a call to update() every generation of the algorithm.
    Once termination is satisfied, self.terminate is set to True, returned by done().
    """
    
    def __init__(self, n_s=20, n_p=3, n_b=10, DEBUG=False):
        """
        Initialize termination detection algorithms
        Input:
        n_s: Number of generations to compare between
        n_p: Decimal places to keep after rounding mean M and std S during dissimilarity
        n_b: Number of bins for the multidimensional histogram

        """
        self.DEBUG = DEBUG
        
        self.n_s = n_s # Number of generations to compare between
        self.n_p = n_p # Decimal places to keep after rounding mean M and std S during dissimilarity
        self.n_b = n_b # Number of bins for the multidimensional histogram

        self.c_1 = False # Termination criteria: mean dissimilarity has stabilized
        self.c_2 = False # Termination criteria: standard deviation has stabilized
        self.terminate = False # Termination state
        
        self.D = [] # Dissimilarity vector
        self.M = [] # Mean Dissimilarity vector
        self.S = [] # Mean Standard deviation vector
        
    def done(self):
        """ Returns true when termination has been reached """
        return self.terminate

    def status(self, t):
        """ Print vectors D, M and S for a generation t"""
        #if self.DEBUG:
        print "t:", t, "| D:",self.D[t], "M:",self.M[t], "S:",self.S[t]

    def update(self, P, Q, t):
        """
        MOEA Termination Dectection Algorithm (Algorithm 2)
        Input:
        P - Population at time t
        Q - Population at time t+1
        t - Generation counters                   
       
        For generation t, calculates dissimilarity between non-dominated pareto fronts P and Q 
        Update vectors for dissimilarity, mean dissimilarity and std dissimilarity (self.D, self.M and self.S)
        Check if last n_s generations have stabilized.
        """
        
        C,C_q,P_c,Q_c,Q_cq = self.MultiHistogram(P,Q)
        
        self.D.append( 0.0 ) # Dissimilarity
        for i in xrange(len(C)):
            p = float(P_c[i])/len(P)
            q = float(Q_c[i])/len(Q)
            if p == 0: continue
            
            if q > 0:
                self.D[t] = self.D[t] - ( ( (p/2)*log(p/q)  ) + ( (q/2)*log(p/q) ) )
            elif q == 0:
                self.D[t] = self.D[t] - p*log(p)
        for i in xrange(len(C_q)):
            if q == 0: continue
            q = float(Q_cq[i])/len(Q)
            self.D[t] = self.D[t] - (q * log(q))

        if t <= 1:
            M_t = self.D[t]
            S_t = M_t
        else:
            # Determine M_t and S_t
            M_t = round((1/float(t)) * sum([D_i for D_i in self.D[-self.n_s:] ]), self.n_p)
            S_t = round((1/float(t)) * sum([ (D_i - M_t)**2 for D_i in self.D[-self.n_s:] ]), self.n_p)
        self.M.append(M_t)
        self.S.append(S_t)

        if t > self.n_s:
            if self.DEBUG:
                print " Testing termination:"
                print "Matching M -> ", self.M[t]
                print "against -> ", self.M[-self.n_s:-1]
            if all(x == self.M[t] for x in self.M[-self.n_s:-1]): self.c_1 = True
            
            if self.DEBUG:
                print " Testing termination:"
                print "Matching S -> ", self.S[t]
                print "against -> ", self.S[-self.n_s:-1]
            if all(x == self.S[t] for x in self.S[-self.n_s:-1]): self.c_2 = True
        
        if all( [self.c_1, self.c_2] ):
            self.terminate = True
            if self.DEBUG:
                print "Termination reached."
        else:
            self.c_1 = False
            self.c_2 = False
            # And implicitely, t += 1 and P <- Q
        return
    
    def MultiHistogram(self, P, Q):
        """
        Multidimensional Histogram Algorithm for Two MOEA Population (Algorithm 1)

        Input: 
        P - non-dominated population at time t
        Q - non-dominated population at time t+1
        self.n_b - Number of bins used to partition the search space

        Output:
        C - List of cells with solutions from the intersection between P and Q and those that are in P but not Q.
        C_q - List of cells representing contained solutions of Q but not P

        P_c - List of numbers of solutions in each cell of C belongning to population P
        Q_c - List of number of solutions for each cell of C belongning to population Q
        Q_cq - List of number of solutions for each cell of C_q belongning to population Q

        Description:
        Assembles the data structure to partition the M-dimensional objective space of solutions between two populations.
        See Section IV. Mutidimensional Histogram Algorithm for MOEA Populations
        """
        C = []
        C_q = []
        P_c = []
        Q_c = []
        Q_cq = []
         
        for s in P:
            c = self.GetCell_id(s, P, Q)
            if c in C:
                k = C.index(c)
                P_c[k] += 1
            else:
                C.append(c)
                P_c.append(1)
                Q_c.append(0)

        for s in Q:
            c = self.GetCell_id(s, P, Q)
            if c in C:
                k = C.index(c)
                Q_c[k] += 1
            else:
                if c in C_q:
                    k = C_q.index(c)
                    Q_cq[k] += 1
                else:
                    C_q.append(c)
                    Q_cq.append(1)

        return C,C_q,P_c,Q_c,Q_cq

    def GetCell_id(self, s, P, Q):
        """
        Unique cell identification ( Equation 7 )

        Input:
        s - a fitness tuple of dimension M
        P - Population at generation t
        Q - Population at generation t+1

        Output:
        c - a unique cell identifier representing a cell of the histogram partitioning. 

        Description:
        Identifies a cell position for  s in the M-dimension space from (P union Q)
        See Section IV. Mutidimensional Histogram Algorithm for MOEA Populations: B Assignment of Unique Identification to Cell
        """
        union = P + Q
        _s = []
        B = [x/self.n_b for x in xrange(self.n_b+1)]
        def get_k(_sj, B):
            k = 0
            for b in B:
                if _sj > b:
                    k+=1
                else: break
            return k
            
        #for s in union:
        for _ in [1]:
            for j in xrange(len(s)):
                Omin_j = min([x[j] for x in union])
                Omax_j = max([x[j] for x in union])
                if Omax_j == Omin_j:
                    _sj = 1.0
                else:
                    _sj = (s[j] - Omin_j)/(Omax_j - Omin_j)
                _s.append(_sj)

            c = 0
            for j in xrange(len(s)):                
                k_j = get_k(_s[j], B)
                c += k_j * (self.n_b**j) # Paper has n_b ** (j-1) but numbers j = 1 .. M
        return c

#================================================================#
#================================================================#
if __name__ == "__main__":
    """
    Example
    """
    import random
    import sys

    if len(sys.argv) == 1:        
        n_s=10
        n_p=3
        n_b=10
    elif sys.argv[1] == "-h":
        print "Usage:"
        print "python",sys.argv[0], "n_s", "n_p", "n_b"
        print "n_s : Number of generations to compare between."
        print "n_p : Decimal places to keep after rounding mean M and std S during dissimilarity"
        print "n_b : Number of bins for the multidimensional histogram"
        exit(0)
    else:
        n_s,n_p,n_b = [int(x) for x in sys.argv[1:]]
        
    terminator = MOEATerminationDetection( n_s=n_s, n_p=n_p, n_b=n_b )
    #terminator.DEBUG = True
    
    a_generation = [    [0.9, 0.9],
                        [0.2, 0.5],
                        [0.2, 0.1]
    ]
    
    b_generation = [    [0.1, 1.0],
                        [0.5, 0.5],
                        [1.0, 0.1]
    ]

    g = []
    for i in range(50):
        # Simulated 3 M=2 fitness tuples.
        g.append( [[random.random(),random.random()],
                   [random.random(),random.random()],
                   [random.random(),random.random()]]  ) 

    # Add many of the same population
    generations = g + [ b_generation ]*100

    t = 0    
    print "Input generations:"
    print "\n".join([str(x) for x in generations])
    while not terminator.done() and t != len(generations):
        terminator.update(generations[t-1], generations[t], t)
        terminator.status(t)

        """
        print "D", terminator.D
        print "M", terminator.M
        print "S", terminator.S    
        """
        
        if terminator.done():
            print "Done. Exiting"
            exit(0)
            break
        t+=1
print "Automatic termination was not reached."

