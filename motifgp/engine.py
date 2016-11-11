import copy
import math
import numpy
import operator
import random
import deap
import uuid
import os

from utils import Utils

import evoalgo
import nsgafortin
import reseed
import primitives
from moead import MOEAD
#import MOODS

#import matrixevaluator
from stgpfitness import STGPFitness
from grammars import AlphaGrammar
from grammars import IUPACGrammar
from grammars import AlphaRangeGrammar
from grammars import IUPACRangeGrammar
from grammars import NetworkExpressionGrammar
from grammars import PSSMGrammar
from grammars import ConditionalGrammar

PATHOS = False
try:
    from pathos.multiprocessing import ProcessingPool as Pool
    from pathos.multiprocessing import cpu_count
    PATHOS = True
except:
    # Cannot use pathos for multiprocessing
    PATHOS = False
    pass

class Engine(STGPFitness):

    NUM_OFFSPRING=75

    def __init__(self, 
                 experiment=None, 
                 options=None,
                 fitness=None, 
                 mstats=None, 
                 logbook=None,
                 pset=None, 
                 creator=None,
                 toolbox=None, 
                 population=None, 
                 invalid_ind=None, 
                 hof=None, 
                 START_GEN=0,
                 NGEN=0,
                 OUTPUT_PATH=None,
                 POP_SIZE=None, 
                 SEED=42,
                 CXPB=0.7,
                 MUTPB=0.3
                 ):

        super(Engine, self).__init__(fitness)

        self.experiment=experiment 
        self.options=options 
        self.mstats=mstats 
        self.logbook=logbook
        self.creator=creator
        self.pset=pset
        self.toolbox=toolbox
        self.population=population
        self.invalid_ind=invalid_ind
        self.hof=hof
        self.POP_SIZE=POP_SIZE 
        self.SEED=SEED

        self.START_GEN = START_GEN
        self.NGEN=0

        self.OUTPUT_PATH = OUTPUT_PATH

        self.positive_path=None
        self.negative_path=None
        self.re_positive_dataset=None
        self.re_negative_dataset=None
        
        self.initial_seeds = []

        self.utils = Utils()


        self.CXPB=CXPB
        self.MUTPB=MUTPB
        
        self.EXPERIMENT="regex"
        #print deap.creator
        #try:
        #    """
        #    Recursive calls to engine may already have created the deap.creator types
        #    """       
        #    deap.creator.FitnessMax
        #except:
        #deap.creator.create("Individual", deap.gp.PrimitiveTree, fitness=deap.creator.FitnessMax, pset=self.pset)    
        #deap.creator.create("FitnessMax", deap.base.Fitness, weights=self.FITNESS_WEIGHTS)    
        
    def boot(self, options):
        self.dna_pset(options)
        self.load_dataset(options)
        if not options.no_seed:
            pass
            #self.seed_population()

    def run(self, NGEN=None):            

        if NGEN:
            self.NGEN=NGEN
        offspring, log = None, None
        if self.options.moo == "NSGAR":
            # Begin the generational process
            
            offspring, log = evoalgo.eaFortin(self.population, 
                                              self.toolbox, 
                                              mu= self.POP_SIZE, #Engine.NUM_OFFSPRING #int(math.ceil(self.POP_SIZE/2.0)),
                                              cxpb=self.CXPB,
                                              mutpb=self.MUTPB,
                                              start_gen=self.START_GEN,
                                              ngen=self.NGEN,
                                              stats=self.mstats,
                                              halloffame=self.hof,
                                              invalid_ind=self.invalid_ind,
                                              logbook=self.logbook,
                                              CPOUT=self.OUTPUT_PATH,
                                              checkpoint=self.options.checkpoint_path,
                                              verbose=True,
                                              timelimit=self.options.timelimit,
                                              termination=self.options.termination
                                              ) 
            print "Memoized",self.memoize_count,"calls."
        elif self.options.moo in ["SPEA2", "NSGA2"]:
            offspring, log = evoalgo.checkpoint_eaMuPlusLambda(self.population, self.toolbox, 
                                                               mu=self.POP_SIZE,
                                                               #lambda_=int(math.ceil(self.POP_SIZE/2.0)),
                                                               #lambda_=int(math.ceil(self.POP_SIZE/2.0)),
                                                               lambda_=self.POP_SIZE,
                                                               cxpb=self.CXPB,
                                                               mutpb=self.MUTPB,
                                                               start_gen=self.START_GEN,
                                                               ngen=self.NGEN, 
                                                               stats=self.mstats, 
                                                               halloffame=self.hof,
                                                               logbook=self.logbook,
                                                               CPOUT=self.OUTPUT_PATH,
                                                               checkpoint=self.options.checkpoint_path,
                                                               verbose=True,
                                                               timelimit=self.options.timelimit,
                                                               termination=self.options.termination)
            print "Memoized",self.memoize_count,"calls."
        elif self.options.moo == "MOEAD":
            MU = self.POP_SIZE
            LAMBDA = 2 #TODO: Make this an input parameter
            
            ea = MOEAD(self.population,
                       self.toolbox,
                       MU,
                       self.CXPB,
                       self.MUTPB,
                       self.logbook,
                       ngen=self.NGEN,
                       timelimit=self.options.timelimit,
                       stats=self.mstats,
                       halloffame=self.hof,
                       termination=self.options.termination,
                       nr=LAMBDA)
            offspring, log = ea.execute()                              
            
            
        elif self.options.moo == "ESCMA":
            offspring, log = evoalgo.cmaES(self.population,
                                           self.toolbox,
                                           mu=self.POP_SIZE,
                                           lambda_=self.POP_SIZE,
                                           ngen=self.NGEN,
                                           stats=self.mstats,
                                           halloffame=self.hof,
                                           logbook=self.logbook)
        return offspring, log



    def dna_pset(self, options):
        """ 
        Calls a grammar depending on the parameter
        """
        try: # FIXME: Limitation; Grammar can only be created once...
            grammar = None
            if options.grammar in ["alpha"]:
                grammar = AlphaGrammar()
            elif options.grammar in ["iupac"]:
                grammar = IUPACGrammar()
            elif options.grammar in ["alpharange"]:
                grammar = AlphaRangeGrammar()
            elif options.grammar in ["iupacrange"]:
                grammar = IUPACRangeGrammar()
            elif options.grammar in ["ne"]:
                grammar = NetworkExpressionGrammar()
            elif options.grammar in ["pssm"]:
                grammar = PSSMGrammar()
            elif options.grammar in ["conditional"]:
                grammar = ConditionalGrammar()
                
            self.pset = grammar.get_pset()

        except Exception as e:
            print e.message            
            print "Failed to create grammar. Using", self.pset
            raw_input()
                
        return self.pset

    def pareto_similar(self, a,b):
        """ Crowding function for the pareto front """
        return (self.toolbox.compile(expr=a) in self.toolbox.compile(expr=b))

    ### Objectives evalution ###
    def accuracy(self, pos_matchs, neg_matchs):
        if pos_matchs == 0: 
            return 0
        else:
            return (pos_matchs)/(1.0*pos_matchs+neg_matchs)

    def coverage(self, pos_matchs, positive_count):
        return pos_matchs/(1.0*positive_count)
    def support_score(pos_matchs):
        return coverage(pos_matchs)
    def length_score(self, lengths):
        if len(lengths) < 1:
            return 0
        return numpy.mean(lengths)

    def fishers_exact_test(self, p, P, n, N):
        pvalue_threshold = 0.5
        return getLogFETPvalue(p, P, n, N, 0.5)

    def transpose_pssm(self , positions):
        """
        Takes a list of N lists that are of lenght 4, representing 
        a position's alphabet probability. This creates a list of
        4 list of length N. Elements of output 'pssm' are vectors
        of probabilities for given nucleotides.
        """
        pssm = []
        for x in zip(*positions):
            nc_vector = []
            for y in x:
                nc_vector.append(y)
            pssm.append(nc_vector)
        return pssm

    def load_dataset(self, options):
        if not self.options:
            self.options = options
        
        """
        self.ME = matrixevaluator.MatrixEvaluator(options.training_path, options.background_path)
        print "Foreground ME loaded with", self.ME.count(), "sequences"
        """
        if options.training_path == None:
            self.re_positive_dataset, self.re_negative_dataset, self.header_positive_dataset,self.header_negative_dataset = [],[],[],[]            
        else:            
            self.re_positive_dataset, self.re_negative_dataset, self.header_positive_dataset,self.header_negative_dataset = self.utils.load_datasets(
                self.options.training_path,
                self.options.background_path,
                self.options.bg_algo,
                self.options.short,
                HEADERS=True,
                HARDMASK=options.hardmask
            )
            #self.re_positive_dataset, self.re_negative_dataset, self.header_positive_dataset,self.header_negative_dataset = list(self.re_positive_dataset), list(self.re_negative_dataset), list(self.header_positive_dataset),list(self.header_negative_dataset)
            
            if not self.options.background_path or self.options.matcher == "grep":
                #TODO: It would be better not to replace the options paths
                runtime_dir = "runtime_tmp"
                if not os.path.isdir(runtime_dir):
                    os.mkdir(runtime_dir)
                unique_id = str(uuid.uuid4())
                
                unique_filename = os.path.join(runtime_dir, os.path.basename(self.options.training_path)+"_background_"+unique_id)
                with open(unique_filename, "w") as f:
                    print "Writing background to" , unique_filename
                    for header, line in zip(self.header_positive_dataset, self.re_negative_dataset):
                        #f.write(header+"_SHUFFLED\n")
                        f.write(line+"\n")
                    
                self.options.background_path=unique_filename                

                unique_filename = os.path.join(runtime_dir, os.path.basename(self.options.training_path)+"_foreground_"+unique_id)
                with open(unique_filename, "w") as f:
                    print "Writing positive to" , unique_filename
                    for header,line in zip(self.header_positive_dataset, self.re_positive_dataset):
                        #f.write(header+"\n")
                        f.write(line+"\n")                                   
                self.options.training_path=unique_filename


        #if options.fastamask:
        #    """
        #    Mask sequences based on a fasta mask format 
        #    """
            
        
        """
        if not self.ME.bgseq and self.re_negative_dataset:
            #If background was created with dinucle shuffle then add it to the matrix evaluator
            self.ME.setBgSeq(self.re_negative_dataset)
        print "Background ME loaded with", len(set(self.ME.bgcoords.values())), "sequences"
        """
                
    def multimutate(self, ind, pset=None):
        select = numpy.random.rand()
        if select < 0.5:
            return deap.gp.mutShrink(ind)
        else:
            return deap.gp.mutInsert(ind, pset=pset)

    def set_experiment(self, experiment, options, mstats, logbook, pset, toolbox, population,  invalid_ind, hof, START_GEN=None, POP_SIZE=None,  SEED=42):        

        if options:
            self.options = options

        if not self.re_positive_dataset:
            print "Loading datasets"
            self.load_dataset(self.options)

        if experiment:
            self.experiment = experiment


        if hof:
            self.hof = hof
        else:
            self.hof = deap.tools.ParetoFront(self.pareto_similar)

        if START_GEN:
            self.START_GEN = START_GEN
        if POP_SIZE:
            self.POP_SIZE = POP_SIZE

        self.CXPB = options.cxpb           
        self.MUTPB = options.mutpb
            

        if mstats:
            self.mstats = mstats
        else:
            stats = {}

            def lambda_factory(idx):
                return lambda ind: ind.fitness.values[idx]
                
            for tag in self.FITNESS_TAGS:
                s = deap.tools.Statistics( key=lambda_factory(
                    self.FITNESS_TAGS.index(tag) 
                ))
                stats[tag] = s
    
            #stats["Discrimination"] = deap.tools.Statistics(key=lambda ind: ind.fitness.values[0])
            #stats["Support"] = deap.tools.Statistics(key=lambda ind: ind.fitness.values[1])

            #stats_accuracy = deap.tools.Statistics(key=lambda ind: ind.fitness.values[0])
            #stats_coverage = deap.tools.Statistics(key=lambda ind: ind.fitness.values[1])

            self.mstats = deap.tools.MultiStatistics( **stats  ) #accuracy=stats_accuracy, coverage=stats_coverage )

            self.mstats.register("avg", numpy.mean, axis=0)
            self.mstats.register("std", numpy.std , axis=0)
            self.mstats.register("min", numpy.min , axis=0)
            self.mstats.register("max", numpy.max , axis=0)

        if logbook:
            self.logbook = logbook
        else:
            self.logbook = deap.tools.Logbook()
            #self.logbook.header = "gen", "evals", "hypervolume", "memoize"
            self.logbook.header = "gen", "evals", "memoize"
            self.logbook.header += tuple((self.FITNESS_TAGS)) #tuple(["accuracy","coverage"]) #

        if pset:
            self.pset = pset

        if population:
            self.population = population

        if invalid_ind:
            self.invalid_ind = invalid_ind

        #creator.create("FitnessMax", base.Fitness, weights=self.FITNESS_WEIGHTS)
        #creator.create("Individual", gp.PrimitiveTree, fitness=creator.FitnessMax, pset=self.pset)

        if toolbox:
            self.toolbox = toolbox
        else:
            self.toolbox = deap.base.Toolbox()
	    if self.options.ncpu and ( self.options.ncpu == "auto" or int(self.options.ncpu) > 1):
                """
                Check if PATHOS was loaded
                """
                global PATHOS
                if PATHOS:
                    if self.options.ncpu == "auto":
		        n_cpu = cpu_count()
                    else:
                        n_cpu = int(self.options.ncpu)
		    print "Using pathos.multiprocessing with", n_cpu, "cpus."
                    pool = Pool(n_cpu)
                    self.toolbox.register("map", pool.map)
                else:
                    print "WARNING: Failed to load pathos; using a single processor."

            """ Basic toolbox functions"""
            self.toolbox.register("expr", deap.gp.genHalfAndHalf, pset=self.pset, type_=self.pset.ret, min_=1, max_=4)
            self.toolbox.register("individual", deap.tools.initIterate, deap.creator.Individual, self.toolbox.expr)
            self.toolbox.register("population", deap.tools.initRepeat, list, self.toolbox.individual)
            self.toolbox.register("compile", deap.gp.compile, pset=self.pset)
            
            # Set evaluation 
            if self.options.grammar == "pssm":
                self.toolbox.register("evaluate", self.eval_pssm_match)
            else:
                try:
                #self.toolbox.register("evaluate", self.memoize_pfm_match)
                #self.toolbox.register("evaluate", self.eval_regex_match)
                #self.toolbox.register("evaluate", self.memoize_regex_match)
                #self.toolbox.register("evaluate", self.memoize_python_matcher)
                    if self.options.matcher == "grep":
                        print "Loading grep matcher"
                        self.test_grep()
                        self.toolbox.register("evaluate", self.memoize_grep_match)
                    elif self.options.matcher == "python" :
                        print "Loading Python matcher"
                        self.toolbox.register("evaluate", self.memoize_python_matcher)
                    else:
                        raise Exception("Unknown matching mode" + str(self.options.matcher))
                    print "Loaded."
                except Exception as e:
                    print "Exception occured while setting matcher: ({0}): {1}".format( e.message, str(e.args))
                    print "Use --matcher=python as argument."
                    exit(-1)
                    #print "Falling back to python matcher"                    
                    #self.toolbox.register("evaluate", self.memoize_python_matcher)
                        
                    

                # Add penalties
                for penalty_feasible in self.FITNESS_PENALTIES:
                    if penalty_feasible is not None:
                        self.toolbox.decorate("evaluate", 
                                              deap.tools.DeltaPenality(penalty_feasible, (0.0,)*len(self.FITNESS_TAGS)))
                
                #self.toolbox.register("evaluate", self.eval_pfm_match)
                #self.toolbox.register("evaluate", self.eval_pfm_match)
                #self.toolbox.register("evaluate", self.eval_precomputed_pfm_matches)
            
            #self.toolbox.register("mate", deap.gp.cxOnePoint)
            self.toolbox.register("mate", deap.gp.cxOnePointLeafBiased, termpb=0.1 )
            self.toolbox.register("expr_mut", deap.gp.genGrow, min_=1, max_=4, pset=self.pset)
            self.toolbox.register("mutate", deap.gp.mutUniform, expr=self.toolbox.expr_mut, pset=self.pset)
            #self.toolbox.register("mutate", self.multimutate, pset=self.pset)


            """ Static Limit for GP Tree """
            MAX_HEIGHT = 17 # Koza 
            #MAX_HEIGHT = 8 
            #MAX_HEIGHT = 89
            #MAX_HEIGHT = 25
            self.toolbox.decorate("mate", deap.gp.staticLimit(operator.attrgetter('height'), MAX_HEIGHT))
            self.toolbox.decorate("mutate", deap.gp.staticLimit(operator.attrgetter('height'), MAX_HEIGHT))

            """ Seeding """
            self.toolbox.register("init_seeds", self.initSeeds, deap.creator.Individual)
            self.toolbox.register("init_seeded_population", self.initSeededPopulation, list, self.toolbox.init_seeds)                                                
            
            self.population = self.toolbox.population(self.POP_SIZE) #Random pop for algorithms that need a initial pop
            self.toolbox.register("memoizecount", lambda: self.memoize_count)

            if self.options.moo == "SPEA2":
                self.toolbox.register("select", deap.tools.selSPEA2)
            elif self.options.moo == "NSGA2":
                self.toolbox.register("select", deap.tools.selNSGA2)
            elif self.options.moo == "MOEAD":
                pass
            else:
                if self.options.moo == "ESCMA":
                    from deap import cma
                    """
                    # We start with a centroid in [-4, 4]**D
                    sigma = 2 * 10**(-2 * numpy.random.rand())
                    strategy = cma.StrategyMultiObjective(self.population, 
                                                          #centroid=numpy.random.uniform(-4, 4, len(self.population)), 
                                                          sigma=sigma, 
                                                          lambda_=len(self.population)
                                                     )

                    #self.toolbox.register("generate", strategy.generate, deap.creator.Individual)
                    self.toolbox.register("update", strategy.update)
                    for ind in self.population:
                        ind.fitness.values = self.toolbox.evaluate(ind)
                    """

                    strategy = cma.StrategyMultiObjective(self.population, sigma=1.0, mu=self.POP_SIZE, lambda_=self.POP_SIZE)
                    self.toolbox.register("generate", strategy.generate, deap.creator.Individual)
                    self.toolbox.register("update", strategy.update)
                    pass
                    
                else: #DEFAULT
                    if self.options.moo != "NSGAR":
                        print "Unknown algorithm", self.options.moo, ". Defaulting to NSGAR"
                    self.toolbox.register("preselect", nsgafortin.selTournamentFitnessDCD)
                    self.toolbox.register("select", nsgafortin.selNSGA2)
                    # Evaluate the individuals with an invalid fitness
                    self.population = self.toolbox.population(self.POP_SIZE) #Random pop
                    invalid_ind = [ind for ind in self.population if not ind.fitness.valid]
                    fitnesses = self.toolbox.map(self.toolbox.evaluate, invalid_ind)
                    for ind, fit in zip(invalid_ind, fitnesses):
                        ind.fitness.values = fit
                    self.population = self.toolbox.select(self.population, self.POP_SIZE)
                
            # Set REVCOMP
            self.REVCOMP = options.revcomp


        if not self.options.no_seed:

            seed_options = copy.deepcopy(self.options)
            seed_options.no_seed=True
            seed_options.grammar="alpha"
            #seed_options.moo="SPEA2" #FIXME There's currently a bug preventing seeds from being produced by the FORTIN selection algorithm # Nope, doesn't work at all, just fix it already.
            seedengine = Engine(None)
            #seedengine
            #seedengine.boot(self.options)
            #print self.population
            #raw_input()
            seedengine.set_experiment("regex", seed_options,None, None, None, self.toolbox, self.toolbox.population(self.POP_SIZE), None, None, 0)
            ENGINE_GEN = 50
            print "Seeded for", ENGINE_GEN," generations using grammar:", seed_options.grammar

            self_seeds, self.logbook = seedengine.run(ENGINE_GEN)
            self_seeds = [ind for ind in seedengine.hof]
            #print self_seeds
            #raw_input()
            self.population = self_seeds #TODO: Problem with small seed counts
            pass
        else:
            if len(self.initial_seeds) > 0:
                #self.initial_seeds = ["add,C,A", "C", "add,A,add,C,add,T,G" ]
                print "Seeded population initated"
                self.population = self.toolbox.init_seeded_population()                 
            elif not self.population:
                self.population = self.toolbox.population(self.POP_SIZE)
                
        print "Initial pop size:", len(self.population)
                
    def seed_population(self):
        #print "#SEEDS: ", len(max_rpt_seeds)
        self.max_rpt_seeds()
        self.write_seeds(self.options.training_path, self.initial_seeds)

    def max_rpt_seeds(self):
        print "Using max rpt seeds"
        RELEVANCE_THRESHOLD = 0.001

        positive_seeds = [ [x[1],x[0]] for x in reseed.produce_seeds(self.re_positive_dataset, 2, 6)        ]
        negative_seeds = [  [x[1],x[0]] for x in reseed.produce_seeds(self.re_negative_dataset, 2, 6)        ]

        ns_dict = {}
        for x in negative_seeds:
            ns_dict[x[0]] = x[1]
        
        ps_dict = {}
        for x in positive_seeds:
            if x[0] in ns_dict.keys():
                val = x[1]-ns_dict[x[0]]
            else:
                val = x[1]
            if val > 0:
                ps_dict[x[0]] = val


        self.initial_seeds = [[ps_dict[x], x] for x in ps_dict.keys()]
        #self.initial_seeds = [i for i in self.initial_seeds if int(i[0])/len(self.re_positive_dataset)>RELEVANCE_THRESHOLD]
        self.initial_seeds = sorted( self.initial_seeds, key=lambda x: x[0], reverse=True)
        
        #print self.initial_seeds
        #print len(positive_seeds), len(negative_seeds), len(self.initial_seeds)
        #raw_input()
        #Convert to tree-strings for GP
        #return max_rpt_seeds

    def write_seeds(self, path, seeds):
        f = open(path+".seeds",'w')
        for s in seeds:
            f.write(s[1] + "\n")
        f.close()        

    def GP_seeds(self):
        pass

    def initSeeds(self, icls, _content):
        #print "c",_content
        return icls.from_string(_content, pset=self.pset)

    def initSeededPopulation(self, pcls, init_ind):
        seeds = [self.seed_to_treestring(t[1]) for t in self.initial_seeds]
        pre = pcls(init_ind((s)) for s in seeds)
        return pre 

    def seed_to_treestring(self, string):
        """
        For the seeds to be usable by the toolbox and redable by according to the pset,
        a tree string must be created with the format "add, x, y" where x,y are either
        a char or an instance of another string add, w, z where w,z are recursively the
        same as x,y.
        eg. CCGT's treestring could be "add, C, add, C, add, G, T"
        Could be an equivalent tree.
        """
        stack = list(string)
        tree = ""
        while len(stack) > 1:
            c = stack.pop(0)
            tree += "add,"+str(c)+","
            pass
        tree += stack.pop(0)
        return tree    
