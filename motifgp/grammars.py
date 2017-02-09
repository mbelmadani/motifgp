"""
Collection of grammars used in MotifGP's Engine
"""
import operator
import deap
import random

from primitives import *

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
        self.pset = deap.gp.PrimitiveSetTyped("MAIN", [], str, "IN")
        self.pset.addPrimitive(operator.add, [str,str], str)        

        # Characters of selected nucleotide alphabet
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
            
class SingleSpacerIUPACGrammar(Grammar):
    def __init__(self):
        super(SingleSpacerIUPACGrammar, self).__init__()
        self.pset = deap.gp.PrimitiveSetTyped("MAIN", [], SingleSpacerIUPACExpression, "IN")
        self.pset.addPrimitive(operator.add, [str,str], str)        
        
        # Alphabet of CHARs
        for c in self.ALPHABET:
            self.pset.addTerminal(c, str)
            # Alphabet of Booleans

        # IUPAC Tokens
        self.pset.addPrimitive(primitive_charclass, [bool, bool, bool, bool], str)        
        self.pset.addPrimitive(operator.and_, [bool, bool], bool)
        self.pset.addPrimitive(operator.or_, [bool, bool], bool)
        self.pset.addPrimitive(operator.not_, [bool], bool)
        self.pset.addTerminal(True, bool)
        self.pset.addTerminal(False, bool)

        # Spacer Tokens
        ##Hardcoded into spacer        
        RANGE_MIN, RANGE_MAX, INDEX_MIN, INDEX_MAX = 8, 10, 1, 30 # TODO: Problem specific; Parameterize this for future use.
        self.pset.addPrimitive(typed_singlespacer_iupac_expression, [NetworkExpression, Spacer], SingleSpacerIUPACExpression)
        self.pset.addEphemeralConstant("rangeInt", lambda: random.randint(RANGE_MIN, RANGE_MAX), rangeInt)
        self.pset.addEphemeralConstant("indexInt", lambda: random.randint(INDEX_MIN, INDEX_MAX), indexInt)
        self.pset.addPrimitive(typed_spacer, [indexInt, rangeInt, rangeInt], Spacer)
        self.pset.addPrimitive(typed_network_expression, [str], NetworkExpression)

        self.pset.addTerminal(NetworkExpression(""), NetworkExpression)
        self.pset.addTerminal(Spacer(0,0,0), Spacer)

        ## TODO:FIXME: Can cause bloat
        def nothing(number):
            return number
        self.pset.addPrimitive(nothing, [rangeInt], rangeInt)
        self.pset.addPrimitive(nothing, [indexInt], indexInt)
        ##                              
    
class ConditionalIUPACGrammar(Grammar):
    def __init__(self):
        super(ConditionalIUPACGrammar, self).__init__()
        self.pset = deap.gp.PrimitiveSetTyped("MAIN", [], ConditionalIUPACExpression, "IN")
        self.pset.addPrimitive(operator.add, [str,str], str)        

        # Alphabet of CHARs
        for c in self.ALPHABET:
            self.pset.addTerminal(c, str)
            # Alphabet of Booleans
            
        #self.pset.addPrimitive(typed_conditional_iupac_expression, [ConditionalIUPACExpression], ConditionalIUPACExpression)
        self.pset.addPrimitive(typed_conditional_iupac_expression, [str], ConditionalIUPACExpression)
        self.pset.addPrimitive(primitive_conditional_immediate, [str, str], str) # Conditional, immediate.
        self.pset.addPrimitive(primitive_conditional, [str, str, str], str) # Conditional, with gap.
        self.pset.addTerminal("[ACGT]", ConditionalIUPACExpression)
        
        
class ConditionalGrammar(Grammar):
    """
    A grammar for regular expressions with Lookbehind assertions
    """
    def __init__(self):
        super(ConditionalGrammar, self).__init__()        
        
        self.pset = deap.gp.PrimitiveSetTyped("MAIN", [], ConditionalExpression, "IN")

        self.pset.addPrimitive(operator.and_, [bool, bool], bool)
        self.pset.addPrimitive(operator.or_, [bool, bool], bool)
        self.pset.addPrimitive(operator.not_, [bool], bool)

        self.pset.addPrimitive(operator.add, [Nucleotide,Nucleotide], Nucleotide)
                
        self.pset.addTerminal(True, bool)
        self.pset.addTerminal(False, bool)
        
        # Alphabet terminals
        for c in self.ALPHABET:            
            C = typed_conditional_anchor(c)
            self.pset.addTerminal(C, Anchor)
            self.pset.addTerminal(Nucleotide(c), Nucleotide)

        CC = CharacterClass("[ACGT]")
        CONDITIONAL_POSITION = typed_conditional_position(CC, CC)
        CONDITIONAL_PORTION =  CONDITIONAL_POSITION
        
        CPos = ConditionalPortion(CONDITIONAL_POSITION)
        CPor = ConditionalPosition(CPos)

        self.pset.addTerminal(CC, CharacterClass)
        self.pset.addTerminal(CPos, ConditionalPosition)
        self.pset.addTerminal(CPor, ConditionalPortion)

        self.pset.addPrimitive(typed_conditional_expression,
                               [Anchor, ConditionalPortion],
                               ConditionalExpression)
        
        self.pset.addPrimitive(typed_conditional_portion,
                               [ConditionalPosition, ConditionalPosition],
                               ConditionalPortion)

        self.pset.addPrimitive(typed_conditional_position,
                               [CharacterClass, CharacterClass],
                               ConditionalPosition)
        
        self.pset.addPrimitive(typed_conditional_charclass,
                               [bool, bool, bool, bool],
                               CharacterClass)

        self.pset.addPrimitive(typed_conditional_anchor,
                               [Nucleotide],
                               Anchor)

        #self.pset.addPrimitive(typed_conditional_nucleotide,
        #                       [str],
        #                       Nucleotide)
        
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


