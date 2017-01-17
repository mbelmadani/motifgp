class PrimitivePSSM():
    """
    A class representing multiple positional-probability array values
    Basically this class is only here to prevent the GP grammar from creating
    recursive matrices.
    """

    def __init__(self,
                 positions):
        self.positions = positions

    def __add__(self, x):
        if len(self.positions) < 0:
            self.positions = [x]
        else:
            self.positions.append(x)

class Anchor(str): pass

class ConditionalExpression(str): pass
class ConditionalIUPACExpression(str): pass
class ConditionalPortion(str): pass
class ConditionalPosition(str): pass
    
class CharacterClass(str): pass
class Nucleotide(str): pass

class NetworkExpression(str):pass
class SingleSpacerIUPACExpression(str):pass
class rangeInt(int):pass
class indexInt(int):pass

class Spacer(object):

    def __init__(self, index, minimum, maximum, optional=False):
        self.index = index
        self.minimum = minimum
        self.maximum = maximum
        self.optional = optional

    def __repr__(self):
        #return "<Spacer(): index:%s, minimum:%s, maximum:%s>" % (self.index, self.minimum, self.maximum)
        return "typed_spacer(%s,%s,%s)" % (self.index, self.minimum, self.maximum)
    
    def __str__(self):
        return self.toString()

    def toString(self):
        return ".{"+str(self.minimum)+","+str(self.maximum)+"}"

    def insertInto(self, string):
        if len(string) <= 1 : return string
        print "STRING",string
        print "SELF.TOSTRING", self.toString()
        optional = "?" if self.optional else ""
        space = str(self.toString())+optional
        safeIndex = (self.index % (len(string)) ) + 1
        expression = string[:safeIndex]+space+string[safeIndex:]
        print "EXPRESSION",expression, self.index , "out of", safeIndex, self.minimum, self.maximum
        return expression
    
    
## Conditionals
def typed_conditional_iupac_expression(E):
    """ returns a string as a type Nucleotide """
    #print "CREATING REGEX ", E
    c = E.count("@")
    corrected_str = E
    for i in xrange(1, c+1):
        corrected_str = corrected_str.replace("@", str(i), 1)
    #print "REGEX CREATED ", corrected_str
    #raw_input()
    return ConditionalIUPACExpression(corrected_str)

def typed_conditional_charclass(A,C,G,T):
    """ primitive_charclass functions that returns a type CharacterClass """
    return CharacterClass(primitive_charclass(A,C,G,T))
    
def typed_conditional_nucleotide(N):
    """ returns a string as a type Nucleotide """    
    return Nucleotide(N)

def typed_conditional_expression(anchor, rest):
    expression = anchor + rest
    c = expression.count("@")

    for i in xrange(1, c+1):
        expression = expression.replace("@", str(i), 1)
    return expression

def typed_conditional_portion(C1, C2):
    return ConditionalPortion(C1 + C2)

def typed_conditional_position(CC1, CC2):
    return ConditionalPosition("((?(@)"+CC1+"|"+CC2+"))")

def typed_conditional_anchor(N):
    return Anchor("("+N+")?")

## Spacers
def typed_singlespacer_iupac_expression(string, spacer):
    print "spacer:",spacer, " into ", string
    return SingleSpacerIUPACExpression(spacer.insertInto(string))

def typed_network_expression(string):
    return NetworkExpression(string)

def typed_spacer(index, minimum, maximum):
    first = min([minimum, maximum])
    second = max([minimum, maximum])
    return Spacer(index, first, second)
    

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

    if len(charclass)==0:
       args = ['A','C','G','T']		

    for x in ALPHABET:
        if x in args:
            charclass += x
    return "["+charclass+"]"
            
def primitive_charclass(A,C,T,G):
    charclass = ""
    if not any([A,C,G,T]):
        charclass = "ACGT"
    if A: charclass += "A"
    if C: charclass += "C"
    if G: charclass += "G"
    if T: charclass += "T"        

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

def primitive_conditional_immediate(condition, consequence):
    pre = "("+condition+")?"
    consequence = "(?(@)"+consequence+")"
    return pre + consequence

def primitive_conditional(condition, string, consequence):
    pre = "("+condition+")?"
    consequence = "(?(@)"+consequence+")"
    return pre + string + consequence
