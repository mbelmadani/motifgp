"""
An interface for STGP Motif Discovery algorithms
"""

import functools
import objectives 
import os
import deap

from utils import Utils
from instancer import Instancer

from subprocess import Popen, PIPE
import commands

class memoize(object):
    # From codysoyland / methodmemoize.py - https://gist.github.com/267733/8f5d2e3576b6a6f221f6fb7e2e10d395ad7303f9
    def __init__(self, func):
        self.func = func
        self.memoized = {}
        self.method_cache = {}

    def __call__(self, *args):
        return self.cache_get(self.memoized, args,
                              lambda: self.func(*args)
        )
    def __get__(self, obj, objtype):
        return self.cache_get(self.method_cache, obj,
                              lambda: self.__class__(functools.partial(self.func, obj
                                                                   )
                                                 )
        )

    def cache_get(self, cache, key, func):
        try:
            v = cache[key]
        except KeyError:
            cache[key] = func()
            v = cache[key]
        return v

"""
class multimemoize(memoize):
    def cache_get(self, cache, key, func):
        try:
            v = cache[key]
            #print 0
            return v
        except KeyError:
            #print 1
            cache[key] = func()
            v = cache[key]
            return v
   ""
        get = []
        new = []
        for k in key:
            try:
                v = cache[k]
                get.append(v)
            except KeyError:
                new.append(v)
            print key
            raw_input()
        if new:
            new_results = func()
            ## Todo ..... continuie here
   ""
"""

class STGPFitness(object):

    def __init__(self, fitness):

        self.memoize_count = 0
        self.precompute = {}
        self.FITNESS_TAGS = []
        self.FITNESS_FUNCTIONS = []
        self.FITNESS_WEIGHTS = []
        self.FITNESS_0 = []
        self.FITNESS_PENALTIES = []
        self.add_objectives_from_str(fitness)         
        self.REVCOMP=False

        self.utils = Utils()
        self.options = {} # Expected to be overriden
        
        self.grep_env = dict(os.environ)
        self.grep_env['LC_ALL'] = 'C'
        
        print self.FITNESS_FUNCTIONS
        print self.FITNESS_WEIGHTS

        try:
            deap.creator.FitnessMax
        except:
            deap.creator.create("FitnessMax", deap.base.Fitness, weights=self.FITNESS_WEIGHTS)    
            pass
        

    def add_objectives_from_str(self, objectives_str):
        """
        Creates objectives from the configuration map in "config/objectives.map"
        Objectives must be defined in objective.py
        """
        instance_objectives = []
        I = Instancer("config/objectives.map")
        
        for objective in [x.upper() for x in objectives_str]:
            instance_objectives.append( I.construct(objective)().get() )
            """
            if objective == 'F':
                instance_objectives.append(objectives.Fisher().get())
            elif objective == 'I': #Scipy Fisher
                instance_objectives.append(objectives.ScipyFisher().get())
            elif objective == 'D':
                instance_objectives.append(objectives.Discrimination().get())
            elif objective == 'S':
                instance_objectives.append(objectives.Support().get())        
            elif objective == 'O':
                instance_objectives.append(objectives.OddsRatio().get())
            elif objective == 'R':
                instance_objectives.append(objectives.ScipyOddsRatio().get())
            elif objective == 'Q':
                instance_objectives.append(objectives.FalseDiscoveryRate().get())
            """
        self.add_objectives(instance_objectives)

    def add_objectives(self, objective_list):
        print objective_list
        for t,f,w in objective_list:
            self.FITNESS_TAGS.append(t)
            self.FITNESS_FUNCTIONS.append(f)
            self.FITNESS_WEIGHTS.append(w)

            # TODO: Add FITNESS_0 to Objective class
            if w < 0.0:
                self.FITNESS_0.append(1.0)
            else:
                self.FITNESS_0.append(0.0)

    def memoize_grep_match(self, program):
        """
        Wrapper for grep system call using memoizing
        """
        if len(program) < 1:
            return self.FITNESS_0
        
        pattern = self.toolbox.compile(expr=program)        
        if len(pattern) < 1:
            return self.FITNESS_0

        if 'X' in pattern: # The 'pathogenic' token
            return self.FITNESS_0 
        
        if self.REVCOMP:
            pattern = self.utils.add_reverse_complement(pattern)
            
        fitness = self.memoize_or_grep(pattern)
        
        return fitness


    def do_grep(self, params):
        args = ["egrep", "-c"] + params
        output, error = Popen(args, stdin=PIPE, stdout=PIPE, stderr=PIPE ).communicate()
        return int(output.strip())
        
    def test_grep(self):
        """
        Test to see if grep works
        """
        pos_matchs = self.do_grep(["MotifGP", "motifgp.py"])
        
        if pos_matchs > 0:
            return True
        else:
            raise Exception("grep did not work as expected. Use --matcher=python")
    
    @memoize
    def memoize_or_grep(self, pattern):
        """
        Do a system call to get matches by grep 
        """
        self.memoize_count+= 1
        ppath = self.options.training_path
        npath = self.options.background_path
        
        pos_matchs = self.do_grep([pattern, ppath])
        neg_matchs = self.do_grep([pattern, npath])
        
        fitnesses = []

        vals = [ pos_matchs, 
                 len(self.re_positive_dataset),
                 neg_matchs, 
                 len(self.re_negative_dataset),
                 pattern
        ]
        for f in self.FITNESS_FUNCTIONS:
            fitnesses.append( f(*vals) ) 

        return tuple(fitnesses)


    def memoize_singlespace_python_matcher(self, program):
        """
        Wrapper for the singlespacer python regex matcher using memoizing        
        """
        if len(program) < 1:
            return self.FITNESS_0
        
        compiled = self.toolbox.compile(expr=program)
        pattern = compiled.spacer.insertInto(compiled.NE)
        #print "PATTERN", pattern
        if len(pattern) < 1:
            return self.FITNESS_0
        
        if 'X' in pattern: # The 'pathogenic' token
            return self.FITNESS_0         
        
        if self.REVCOMP:
            try:
                pattern = self.utils.add_reverse_complement(pattern)
            except Exception as e:
                print "Exception in stgpfitness.py"
                print e, e.message
                print "pattern:", pattern
                raise e

        fitness = self.memoize_or_python_match(pattern)
        
        return fitness
    
    def memoize_python_matcher(self, program):
        """
        Wrapper for the python regex matcher using memoizing        
        """
        if len(program) < 1:
            return self.FITNESS_0
        
        pattern = self.toolbox.compile(expr=program)        
        if len(pattern) < 1:
            return self.FITNESS_0

        if 'X' in pattern: # The 'pathogenic' token
            return self.FITNESS_0 
        
        if self.REVCOMP:
            pattern = self.utils.add_reverse_complement(pattern)

        fitness = self.memoize_or_python_match(pattern)
        
        return fitness
    
    @memoize
    def memoize_or_python_match(self, pattern, DEBUG=False):
        """
        Returns either the cached results or the python matcher's results.
        """
        self.memoize_count+= 1
        regex = self.utils.list_to_regex_obj(pattern)
        pos_matchs = 0
        neg_matchs = 0
        
        for s in self.re_positive_dataset:
            m = regex.search(s) 
            if m:
                pos_matchs += 1
                
        for s in self.re_negative_dataset:
            m = regex.search(s)
            if m:
                neg_matchs += 1
        
        fitnesses = []
        vals = [ pos_matchs, 
                 len(self.re_positive_dataset),
                 neg_matchs, 
                 len(self.re_negative_dataset),
                 pattern
        ]

        if DEBUG:
            print vals

        for f in self.FITNESS_FUNCTIONS:
            fitnesses.append( f(*vals) ) 

        return tuple(fitnesses)    


    @memoize
    def memoize_or_match(self, pattern):
        self.memoize_count+= 1

        pos_matchs = 0
        support = []

        for idx in xrange(len(self.re_positive_dataset)):
            s = self.re_positive_dataset[idx]
            regex = self.utils.list_to_regex_obj(pattern)
            m = regex.search(s) 
            if m:
                pos_matchs += 1

        neg_matchs = 0
        for s in self.re_negative_dataset:
            regex = self.utils.list_to_regex_obj(pattern)
            m = regex.search(s)
            if m:
                neg_matchs += 1
        
        fitnesses = []
        vals = [ pos_matchs, 
                 len(self.re_positive_dataset),
                 neg_matchs, 
                 len(self.re_negative_dataset) 
        ]

        for f in self.FITNESS_FUNCTIONS:
            fitnesses.append( f(*vals) ) 

        return tuple(fitnesses)

    def reset_hits():
        """ Resets the cache hit count """
        self.cache_hits = 0

    def eval_pssm_match(self, program):
        """
        RETIRED
        """
        if len(program) < 1:
            return (0.0, 0.0,)
        
        pre_pssm = self.toolbox.compile(expr=program)
        pssm = self.transpose_pssm(pre_pssm)

        if len(pssm) < 1:
            return (0.0, 0.0,)
        
        pos_matchs = 0
        matrices = [pssm]
        support, score = self.ME.match(matrices)

        return ( 
            support,
            score,
            )


    @memoize
    def memoize_pfm_helper(self, pattern):
        """
        RETIRED
        """
        self.nevals+= 1

        pos_matchs = 0
        support = []
        lengths = []
        seq = ""
        pfm = self.ME.regex_to_pfm(pattern)
        pfm = [pfm]
        #print pfm
        p = self.ME.match(pfm)
        n = self.ME.match(pfm, negative=True)

        return ( 
            self.accuracy(p,n), 
            self.coverage(p,self.ME.count()), 
            )

    def memoize_pfm_match(self, program):
        """
        RETIRED
        """
        if len(program) < 1:
            return (0.0, 0.0,)
        
        pattern = self.toolbox.compile(expr=program)

        if len(pattern) < 1:
            return (0.0, 0.0,)

        return self.memoize_pfm_helper(pattern)


    
    def _eval_pssm_match(self, program):
        """
        RETIRED
        """
        if len(program) < 1:
            return (0.0, 0.0,)
        
        pre_pssm = self.toolbox.compile(expr=program)
        pssm = self.transpose_pssm(pre_pssm)
        
        """
        print "program", program
        print "pre_pssm", pre_pssm
        print "pssm", pssm
        raw_input()
        """

        if len(pssm) < 1:
            return (0.0, 0.0,)
        
        pos_matchs = 0
        support = []
        #lengths = []
        for idx in range(len(self.re_positive_dataset)):
            s = self.re_positive_dataset[idx]
            m = MOODS.search(s, [pssm], 0.001) #30, convert_log_odds=True, threshold_from_p=True)
            if m[0]:
                pos_matchs += 1
                support.append(idx)
                #span =  m.span()
                #lengths.append( span[1] - span[0] )

        neg_matchs = 0
        for s in self.re_negative_dataset:
            m = MOODS.search(s, [pssm], 0.001) #, convert_log_odds=True, threshold_from_p=True)
            if m[0]:
                neg_matchs += 1

        acc, cov =  self.accuracy(pos_matchs, neg_matchs), self.coverage(pos_matchs, len(self.re_positive_dataset)), 
            
        
        return ( 
            acc,
            cov,            
            )


    def eval_pfm_match(self, program):
        """
        RETIRED
        """
        if len(program) < 1:
            return (0.0, 0.0,)
        
        pattern = self.toolbox.compile(expr=program)
        if len(pattern) < 1:
            return (0.0, 0.0,)

        #print pattern
        #raw_input()

        pos_matchs = 0
        support = []
        lengths = []
        seq = ""
        pfm = self.ME.regex_to_pfm(pattern)
        pfm = [pfm]
        #print pfm

        p = self.ME.match(pfm)
        n = self.ME.match(pfm, negative=True)

        return ( 
            self.accuracy(p,n), 
            self.coverage(p,self.ME.count()), 
            )



    def precompute_pfm_matches(self, population):
        """
        RETIRED
        """
        regexs = []
        idx = 0

        for program in population:
            if len(program) < 1:
                #regexs.append("[ACGT]")
                continue
            else:        
                pattern = self.toolbox.compile(expr=program)
                if pattern not in self.precompute:
                    regexs.append(pattern)
            idx+=1
        
        pfm = [self.ME.regex_to_pfm(pattern) for pattern in regexs]

        # TODO
        for i in range(len(pfm)):
            pass
            """
            if regexs[i] not in self.precompute:
                p = self.ME.match_multiple_pfm[i])
                n = self.ME.match_multiple_pfm[i], negative=True)
                self.precompute[regex[i]] =  (         
                    self.accuracy(p,n), 
                    self.coverage(p,self.ME.count()), 
                )
            """

    def eval_precomputed_pfm_matches(self, pattern, population):
        """
        RETIRED
        """
        try:
            return self.precompute[pattern]
        except:
            self.precompute_pfm_matches(population)
            return self.precompute[pattern]
            


    def eval_regex_match(self, program):
        """
        RETIRED
        """
        if len(program) < 1:
            return (0.0, 0.0,)
        
        pattern = self.toolbox.compile(expr=program)
        
        if len(pattern) < 1:
            return (0.0, 0.0,)
            
        self.nevals+= 1
        
        pos_matchs = 0
        support = []
        lengths = []
        for idx in range(len(self.re_positive_dataset)):
            s = self.re_positive_dataset[idx]
            regex = self.utils.list_to_regex(pattern)
            m = regex.search(s) 
            if m:
                pos_matchs += 1
                support.append(idx)
                span =  m.span()
                lengths.append( span[1] - span[0] )

        neg_matchs = 0
        for s in self.re_negative_dataset:
            regex = self.utils.list_to_regex(pattern)
            m = regex.search(s)
            if m:
                neg_matchs += 1

        return ( 
            self.accuracy(pos_matchs, neg_matchs), 
            self.coverage(pos_matchs, len(self.re_positive_dataset)), 
            )

