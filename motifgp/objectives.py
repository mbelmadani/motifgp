"""
Objective functions that require a number of positive successes (p), negative success (n), and the total in each class (P and N)

Each function should support p, P, n, N as parameter signature.
"""
from hypergeometric import getLogFETPvalue
from math import log
from math import sqrt
from utils import Utils

SCIPY = False
try:
    from scipy.stats import fisher_exact
    SCIPY = True
except:
    # Cannot use scipy-based objectives
    SCIPY = False
    pass
    

class Objective(object):
    def __init__(self):
        self.f = None
        self.w = None
        self.tag = "NoName"

    def get(self):
        return [self.tag, self.f, self.w]

class Discrimination(Objective):
    def __init__(self):
        self.tag = "Discrimination"
        self.w = 1.0

        def fit_discrimination(p, _, n, __, ___):
            if p + n == 0:
                return 0
            else:
                return p/(1.0 * (p+n))
        self.f = fit_discrimination

class Support(Objective):
    def __init__(self):
        self.tag = "Support"
        self.w = 1.0

        def fit_support(p, P, _, __, ___):
            return p/(1.0*P)

        self.f = fit_support

class AbstractFisher(Objective):
    pvalue = None
    OR = None        

    def __init__(self):
        super(AbstractFisher, self).__init__()
        global SCIPY
        if not SCIPY:
            print "WARNING: Scipy-based objectives will not work if Scipy is not installed."
            print "Continue? (y/n)"
            r = raw_input()
            if r.lower() != 'y':
                print "Abort."
                exit(0)
                
        
    def refresh_stats(self, p, P, n, N):
        if p == 0:
            AbstractFisher.OR, AbstractFisher.pvalue = 0.0, 1.0
            return
        else:
            fp = P - p
            fn = N - n
        try:
            AbstractFisher.OR, AbstractFisher.pvalue = fisher_exact([[p,fp],[n,fn]], "greater")
        except Exception as e:
            print "Error!"
            print "p, fp, n, fn"
            print p, fp, n ,fn
            print P, p
            raise e
        
        #Handle ORs
        if AbstractFisher.OR == float('inf'):
            AbstractFisher.OR = float(P + p)
            #AbstractFisher.OR = float(P)
            

class ScipyFisher(AbstractFisher):
    def __init__(self):
        super(ScipyFisher, self).__init__()
        self.tag = "ScFisher"
        self.w = -1.0
                
        def fit_fisher(p, P, n, N, _):
            self.refresh_stats(p,P,n,N)
            return AbstractFisher.pvalue
            
        self.f = fit_fisher

class ScipyOddsRatio(AbstractFisher):
    """
    Currently, this assumes Fisher is also used and precomputed before 
    """
    def __init__(self):
        super(ScipyOddsRatio, self).__init__()
        self.tag = "ScOR"
        self.w = 1.0
                
        def fit_OR(p, P, n, N, _):
            # Assuming we don't need to compute Fisher before
            return self.OR
            
        self.f = fit_OR


class Fisher(Objective):
    def __init__(self):
        self.tag = "Fisher"
        self.w = -1.0
        
        def fit_fisher(p, P, n, N, _):
            #pvalue_threshold = 0.5
            pvalue_threshold = False # Modification to allow p-value above 0.5
            if p < 1:
                return 1
            try:                
                return getLogFETPvalue(p, P, n, N, pvalue_threshold)
            except Exception as e:
                try:
                    # Weird bug with : python hypergeometric.py 5705 5706 2238 2238
                    return getLogFETPvalue(p+1, P, n, N, pvalue_threshold)
                except Exception as e:
                    print p,P,n,N
                    print e
                    raise e
            
        
        self.f = fit_fisher


class FalseDiscoveryRate(Objective):
    def __init__(self):
        self.tag = "FDR"
        self.w = -1.0
        
        def fit_fdr(p, P, n, N, _):
            T = p+n
            if T == 0:
                return 1.0
            else:
                return n/float(p+n)
        self.f = fit_fdr

    
class OddsRatio(Objective):
    def __init__(self):
        self.tag = 'OR'
        self.w = 1.0

        def fit_OR(p, P, n, N, _):
            _p,_n = float( P-p ), float(N-n)
            if _n == 0 or p == 0:
                return 0.0
            elif _p == 0 or n == 0:
                return float(p + P)
            else:
                return (p/_p)/(n/_n)

        self.f = fit_OR

class MCC(Objective):
    """
    Matthew's correlation coefficient.
    """
    def __init__(self):
        self.tag = 'MCC'
        self.w = 1.0

        def fit_MCC(p, P, n, N, _):
            tp = p
            tn = N - n
            fp = n
            fn = P - p
            
            if not all([tp, tn]):
                return -1.0
            else:                            
                mcc = (tp*tn - fp*fn) / sqrt( (tp + fp)*(tp + fn)*(tn + fp)*(tn + fn) )
            return mcc

        self.f = fit_MCC    


class EntropyBalance(Objective):
    """
    Optimizes for the highest amount of information content on each side of non-network expresion substrings.
    """
    # TODO: Fix for multiple spacers and non-DNA-nucleotide alphabet.
    # account for: { } , . * ( ) + N A C T G RNAAlphabet ProteinAlphabet

    def __init__(self):
        self.tag = 'EntropyBalance'
        self.w = 1.0
        u = Utils()

        def entropy(list_of_tokens):
            # Information content normalized by length of the motif
            if len(list_of_tokens) == 0:
                return 0.0
            # TODO: Fix for variable alphabet size
            return sum( [ 4.0 / len(token) for token in list_of_tokens ]  ) 
        
        def fit_EntropyBalance(_,__,___,____,pattern):
            # TODO: Fix for multiple breakable chunks
            left = pattern.split(".")[0]
            right = pattern.split("}")[-1]
            list_pattern_left = u.motif_str_to_list(left)
            list_pattern_right = u.motif_str_to_list(right)
            
            # Return the minimum of both side's normalized entropy
            v = min(entropy(list_pattern_left),
                    entropy(list_pattern_right))
            return v

        self.f =  fit_EntropyBalance
    
