import re
import numpy 

def sxor(s1,s2):    
    # convert strings to a list of character pair tuples
    # go through each tuple, converting them to ASCII code (ord)
    # perform exclusive or on the ASCII code
    # then convert the result back to ASCII (chr)
    # merge the resulting array of characters as a string
    return [ord(a) ^ ord(b) for a,b in zip(s1,s2)]

def hamming_pre_string(regex, sequence):
    """
    To compute the hamming distance, we need to match de regex on sequence and then replace the match with "1"
    """
    match = re.search(regex, sequence)
    if match:
        match = match.group(0)
    else:
        #match = ""
        #"0" * len(sequence)
        return None 

        
    placeholder = "1" * len(match)
    pre_string = list(sequence.replace(match, placeholder))

    for i in range(len(pre_string)):
        if pre_string[i] != '1':
            pre_string[i] = '0'

    return "".join(pre_string)


def compute_hamming(list_of_regex, template, sequence):
    """
    For each regex, create a weighted average from the list of regexs given.


    Matches each regex vs the template on the sequence, calculate the
    hamming distance on the template and adjusts the weight of the 
    result on the length of list_of_regex.

    Sums weighted hammings strings

    Return the overall performance of list_of_regexs vs. template on
    sequence
    """
    hamming_template = hamming_pre_string(template, sequence)

    regexs = None
    if type(list_of_regex) == str:
        regexs = list(list_of_regex)
    else:
        regexs = list_of_regex

    output = None
    for regex in regexs:
        hamming_bs = hamming_pre_string(regex, sequence)
        #print bs1+"$", "\n", bs2+"$"
        #print "".join([str(x) for x in sxor(bs1, bs2)])
        if hamming_bs == None:
            xor_string = [float(x) for x in hamming_template]
            #"".join([str(x) for x in sxor(hamming_template, str("0"*len(hamming_template)))]) # Invert template because no match was found. So match == everything but the template motif.
        else:
            #print hamming_bs, hamming_template
            xor_string = sxor(hamming_bs, hamming_template)
        xor_string = [x/float(len(regexs)) for x in xor_string]
        """
        print ">"
        print xor_string
        print "< "
        """
        if output:
            output = [x + y for x,y in zip(output, xor_string)]
        else:
            output = xor_string
        
    return output
    
def score_hamming(floatstring):
    """
    Converts the weigthed hamming distance string to a numerical value
    """
    return sum( floatstring ) / float(len(floatstring))


class HammingBenchmark():
    """
    Class to contain a benchmark of hamming distance against a synthetic dataset
    """
    def __init__(self):
        self.scores = {}
        """
        self.max = -1
        self.min = -1
        self.mean = -1
        self.std = -1
        self.rmax = -1
        self.rmin = -1
        self.rmean = -1
        self.rstd = -1
        """

    def __repr__(self):
        return "HammingBenchmark()"
    def __str__(self):
        output=""
        for each in self.scores:
            output += each+":\n"
            benchmark_str = [
                #self.scoxres, "\n", 
                "max:",self.scores[each]["max"],
                "min:",self.scores[each]["min"], 
                "mean:",self.scores[each]["mean"],
                "std:",self.scores[each]["std"],

                ]
            output += ",".join([str(x) for x in benchmark_str]) + "\n"
        
        #print benchmark_str
        #output = ",".join(str(x) for x in benchmark_str)
        
        return output

    def compile(self, candidates, sequence_tuples):
        """
        candidates; a batch of regular expression that are to be evaluated 
        sequence_tuples: a list of pairs of templates-sequences
        
        """
        bins = {}
        for each in candidates: # Slice candidates one by one. This can be changes to have a real bin behavior
            sequence_scores = []
            candidate = [each] #TODO:CHANGEME; Quick hack to evaluate each candidate on its own versus the sequence set            
            bins[each] = {}
            bins[each]["score"] = []
            bins[each]["max"] = -1
            bins[each]["min"] = -1
            bins[each]["std"] = -1
            bins[each]["mean"] = -1
            
            for template, sequence in sequence_tuples:
                hamming_str_score = compute_hamming(candidate, template, sequence) 
                candidates_score = tuple((sum(hamming_str_score), score_hamming(hamming_str_score) , hamming_str_score ))
                
                bins[each]["score"].append(candidates_score)

        self.scores = bins
        self.update()

    def update(self):
        for each in self.scores.keys():
            numeric_scores = [x[0] for x in self.scores[each]["score"]]
            #if not numeric_scores:
            #    numeric_scores.append(0)
            self.scores[each]["max"] = max(numeric_scores)
            self.scores[each]["min"] = min(numeric_scores)
            self.scores[each]["std"] = numpy.std(numeric_scores)
            self.scores[each]["mean"] = numpy.mean(numeric_scores)
               
    def flush_data_points(self, ks, xs, outpath, seed, CLEAR=True):
        """
        Prints a data point y such that k[x] = y
        k is an individual. x is the mapping value
        seed will be used to color the datapoint
        outpath is where to append the datapoint. CLEAR overwrites instead of appending.
        """
        
        if CLEAR:
            f = open(outpath, 'w')
        else:
            f = open(outpath, 'a')
        
        for idx in range(len(ks)):
            each = ks[idx]
            x = xs[idx]
            scores = self.scores[each]
            y = scores["mean"],scores["std"],
            output = [str(x) ,  str(y) ,  str(seed)]
            output = "\t".join(output)
            output += "\n"
            print output
            f.write(output)
