# Generator to create multiple patterns from a regex
import altschulEriksonDinuclShuffle #dinuclShuffle
import re
import sys

class Utils():
    
    def __init__(self):
        self.complement = {'A': 'T', 
                           'C': 'G', 
                           'G': 'C', 
                           'T': 'A',
                           '[': ']',
                           ']': '['
                       } 
    
    def readfasta(self, filepath, HARDMASK=False): 
        with open(filepath, 'r') as f: 
            buff = ""
            for line in f:
                if HARDMASK:
                    text = line.strip()
                else:
                    # FIXME: Legacy behavior
                    text = line.upper().strip()

                #print "text",text, len(text)
                if len(text) == 0:
                    continue
                if text[0] == '>':
                    if len(buff) > 0:                 
                        yield buff
                    buff = ""
                    yield text
                elif text[0].upper() in ["A","C","G","T","N"]: 
                    #only add to buffer lines that are sequences
                    buff += text

            if len(buff) > 0:                 
                yield buff

    def produce_sequences(self, regex):
        regex = regex.replace("*", "[ACTG]") # Easier that way...
        regex = regex.strip()
        regex_as_list = []
        current_charclass = ""
        charclass_flag = 0
        for c in regex:
            if c == '[':
                charclass_flag = 1
            elif c == ']':
                charclass_flag = 0
                regex_as_list.append(current_charclass) 
                current_charclass = ""
            elif len(str(c)) == 1:
                if charclass_flag:
                    current_charclass += c
                else:
                    regex_as_list.append(c)
            else:
                print "Rejecting token",c

        list_of_sequences = [""]
        for token in regex_as_list:
            if len(token) > 1:
                old_list = list_of_sequences[:] # Copy old list by value
                old_list_count = len(old_list)
                for i in range(len(token)-1):
                    list_of_sequences.extend(old_list)
                for i in range(old_list_count):  #For each sublist
                    for j in range(len(token)):  #For each item in the charclass..
                        index = j*old_list_count + i
                        list_of_sequences[index] += token[j]
            else:
                for i in range(len(list_of_sequences)):
                    list_of_sequences[i] += token
        return list_of_sequences

    # More or less deprecated...
    def JDonnerRegex(pattern):
        """
        Used for compatibility with J. Donner's Suffix Tree
        """
        return utils.produce_sequences(pattern)

    def list_to_regex_obj(self, a_list):
        """
        Creates a regex objects out of a list of regex tokens
        """
        expr = self.list_to_regex_string(a_list)
        return re.compile(expr)

    def list_to_regex_string(self, a_list):
        """
        Creates a regex string out of a list of regex tokens
        """
        expr = "".join(["["+x+"]" if len(x) > 1 else x for x in a_list])        
        return expr


    def load_datasets(self, 
                      path_pos, 
                      path_neg, 
                      background_alg,
                      SHOULD_SHORT=0,
                      HEADERS=False,
                      HARDMASK=False
                      ):
        positive = []
        negative = []

        header_pos = []
        header_neg = []
        print "HARDMASK:", HARDMASK
        
        #Load positive for path_pos FASTA
        for line in self.readfasta(path_pos, HARDMASK):            
            if HEADERS and line[0] == '>':
                header_pos.append(line)
            else:
                positive.append(line)
            if SHOULD_SHORT > 0 and int(SHOULD_SHORT) == len(positive): 
                break            

        #Load background as shuffled copy of pisitive
        if path_neg != None:
            #Load background from path_neg FASTA
            for line in self.readfasta(path_neg): # Should bg determi hardmasking with the same variable as positive?
                if HEADERS and line[0] == '>':
                        header_neg.append(line)
                else:
                    negative.append(line)
                if SHOULD_SHORT > 0 and int(SHOULD_SHORT) == len(negative): 
                    break               
        elif background_alg == "dinuclShuffle":
            for line in positive:            
                sh_line =  altschulEriksonDinuclShuffle.dinuclShuffle(line.upper())
                negative.append(sh_line)
        else:
            #No background
            negative.append("NNNNNNNNNNNNNNNNNNNNNNNNNNN")

        print "Lengths:", len(positive), len(negative)
        if HEADERS:
            return positive, negative, header_pos, header_neg
        else:
            return positive, negative

    def motif_str_to_list(self, string):
        """
        Convert str to a list of of chars or charclasses
        """
        motif = []
        flag = 0 # Charclass parsing
        buff = ""
        string = string.strip()
        for c in string:
            if flag == 1:
                if c == ']':
                    flag = 0
                    motif.append(buff)
                    buff = ""
                else:
                    buff += c                
            else:
                if c == '[':
                    flag = 1
                else:
                    motif.append(c)
        motif = filter(bool, motif)
        return motif

    """        
            
        string = string.replace('[', ']').split(']')
        string = [ x.strip() for x in string]
        a_list = filter(bool, string)
        return a_list
    """
    def reverse_complement(self, pattern):    
        revcomp = ""
        revcomp = "".join([self.complement.get(c) for c in reversed(pattern)])
        return revcomp

    def add_reverse_complement(self, pattern):
        # Append reverse complment to the pattern
        # Sorting to return consistant key on cache, even if complement is given next.
        merge = "|".join(sorted([pattern, self.reverse_complement(pattern)]))
        return merge

    def motif_eraser(self, dataset, patterns ):
        """
        EXPERIMENTAL: Takes a list of sequences a removes occurences of patterns
        
        Consideration: The order in which patterns are deleted matters.
        """
        for _p in patterns:
            p = re.compile(_p)
            erased = []
        
            for idx in range(len(dataset)):
                m = p.search(dataset[idx])
                if m:
                    replacement = 'N' * len(m.group(0))

                    if p.sub(replacement, dataset[idx]):
                        erased.append(p.sub(replacement, dataset[idx]))
                else: 
                    erased.append( dataset[idx] )

            dataset = erased
        return dataset


    def get_motifs_from_nef(self, nef_path, tag="NE"):
        """
        Finds motifs in a .tsv format (such as .nef or .neft) stored under the field tag.
        Returns a list of motifs
        """

        NETWORK_EXPRESSION_TAG = tag
        motifs = []

        with open(nef_path) as f:
            idx = None
            for line in f:
                fields = line.split("\t")
                if idx:
                    network_expression = fields[idx]
                    motifs.append(network_expression)
                elif NETWORK_EXPRESSION_TAG in fields:
                    idx = fields.index(NETWORK_EXPRESSION_TAG)                   
                else:
                    print "Format error!"
                    print "line:", line
                    sys.exit(-1)
        return motifs
