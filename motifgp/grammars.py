"""
Collection of grammars used in MotifGP's Engine
"""
import operator
import deap
import random


ALPHABET=["A", "C", "G", "T"]

class Grammar(object):
    def __init__(self):
        self.pset = None

        # Alphabet
        self.ALPHABET=ALPHABET

    def get_pset(self):
        try:
            deap.creator.Individual
        except:
            deap.creator.create("Individual", deap.gp.PrimitiveTree, fitness=deap.creator.FitnessMax, pset=self.pset)

        return self.pset

class AlphaGrammar(Grammar):
    def __init__(self):
        super(AlphaGrammar, self).__init__()
        print "Creating pset for Regex"
        self.pset = deap.gp.PrimitiveSetTyped("MAIN", [], str, "IN")
        self.pset.addPrimitive(operator.add, [str,str], str)        

        # Alphabet of CHARs
        for c in self.ALPHABET:
            self.pset.addTerminal(c, str)
            # Alphabet of Booleans

class IUPACGrammar(AlphaGrammar):
    def __init__(self):
        super(IUPACGrammar, self).__init__()
        self.pset.addPrimitive(primitive_charclass, [bool, bool, bool, bool], str)
        
        self.pset.addPrimitive(operator.and_, [bool, bool], bool)
        self.pset.addPrimitive(operator.or_, [bool, bool], bool)
        self.pset.addPrimitive(operator.not_, [bool], bool)

        self.pset.addTerminal(True, bool)
        self.pset.addTerminal(False, bool)

class NetworkExpressionGrammar(AlphaGrammar):
    def __init__(self):
        # Adds nucleotides and concatenantions
        super(NetworkExpressionGrammar, self).__init__()
        # CharClasses
        self.pset.addPrimitive(primitive_str_charclass, [str,str,str,str], str)
        
 
class AlphaRangeGrammar(AlphaGrammar):
    def __init__(self):
        super(AlphaRangeGrammar, self).__init__()
        self.pset.addPrimitive(primitive_range, [str,int,int,str], str)
        RANGE_MIN, RANGE_MAX = 1,10
        self.pset.addEphemeralConstant("rangeInt", lambda: random.randint(RANGE_MIN, RANGE_MAX), int)
        self.pset.addPrimitive(primitive_addrange, [int,int], int)  
 
class IUPACRangeGrammar(AlphaRangeGrammar, IUPACGrammar):
    def __init__(self):
        super(IUPACRangeGrammar, self).__init__()

class PSSMGrammar(Grammar):
    def __init__(self):
        super(PSSMGrammar, self).__init__()
        print " Creating pset for PSSMs"
        self.pset = deap.gp.PrimitiveSetTyped("MAIN", [], list, "IN")
        self.pset.addPrimitive(operator.add, [list,list], list)   
        #self.pset.addPrimitive(operator.add, [int,int], int)  
        
        self.pset.addPrimitive(primitive_position, [int, int, int, int], list)
        self.pset.addEphemeralConstant("randWeight", lambda: random.randint(0, 10), int)
        #pset.addPrimitive(operator.add, [list,list], list)
        #pset.addPrimitive(operator.add, [primitives.PrimitivePSSM, list], primitives.PrimitivePSSM)
        #pset.addPrimitive(operator.add, [list, primitives.PrimitivePSSM], primitives.PrimitivePSSM)
        #pset.addTerminal([[0,0,0,0]], list)
        self.pset.addTerminal([], list)
        #pset.addTerminal(True, bool)
        #pset.addTerminal(False, bool)
        #pset.addTerminal(0, int)
        #pset.addTerminal(1, int)

## Defining Primitives ##
def primitive_addrange(a, b):
    MAX = 10
    MIN = 8
    c = a + b
    if c < MIN:
        c = MIN
    elif c > MAX:
        c = MAX
    return c

def primitive_not(boolean):
    return not boolean

def primitive_str_charclass(one,two,three,four):
    args = [one,two,three,four]
    charclass = ""

    for x in ALPHABET:
        if x in args:
            charclass += x
    if len(charclass)==0:
        return ""
    return "["+charclass+"]"
            
def primitive_charclass(A,C,T,G):
    charclass = ""
    if A: charclass += "A"
    if C: charclass += "C"
    if G: charclass += "G"
    if T: charclass += "T"

    if len(charclass) <= 1:
        return charclass
    else:
        return "["+charclass+"]"

    
def primitive_position(A,C,G,T):
    """
    Creates a position object for a PSSM
    Assuming A,C,T,G for PSSM order
    """
    position = [[A,C,T,G]]
    return position

"""
    def primitive_get_seeds():

    #input: None
    #description: Access the list of of gathered seeds, picks one randomly
    #returns: str

    global epheremal_seeds
    idx = random.randint(0, len( epheremal_seeds)-1)
    return epheremal_seeds[idx]

"""


def primitive_range(pre, a, b, post):
    """
    input: ints, range start and range length, respectively
    """
    start,stop = str(int(a)), str(int(a+b))
    str_range = pre+ \
                ".{"+ \
                start+ \
                ","+ \
                stop+ \
                "}"+ \
                post 
    #Regex token should be something like ".{1, 5}"
    return str_range
