#!/usr/bin/python
"""
MotifGP - Multiobjective Motif Discovery using Genetic Programming

@author Manuel Belmadani

"""
import os
import sys
sys.path.append("./fortin2013/")

#import re2 as re
import re
import math
import random
import numpy
import timeit
import ast

#MotifGP-specific imports
from utils import Utils
from inspector import Inspector
from engine import Engine
from instancer import Instancer

import postprocessing
import statcsv
import evoalgo
import hammingregex
import platform

from optparse import OptionParser
def check_parameter_constraints(options):
    # Inspector and checkpoint
    """
    if options.inspector and not options.checkpoint_path: #DISABLED
        print "-i/--inspector requires -c/--checkpoint parameter"
        sys.exit(-1)
    elif options.inspector:
        options.no_seed=True
        options.moo=None
    """
    # Hamming and pos.file
    if not options.training_path:
        print "Input sequences are required. Use -p <FILE> to specify the path to a FASTA format file of sequences."
        sys.exit(-1)
    if float(options.cxpb + options.mutpb) != 1.0:
        print "EXITING: Crossover and mutation probabilities should add up to 1.0"
        sys.exit(-1)

if __name__ == "__main__":
    parser = OptionParser()
    I = Instancer("config/objectives.map")
    objectives_help = ",".join([a+"="+b for a,b, in I.config.iteritems()])
    
    
    parser.add_option("-p", "--training", dest="training_path",
                      help="Fasta file to use for training (input) sequence data", default=None)

    parser.add_option("-b", "--background", dest="background_path",
                      help="[Optional] Fasta file to use for background (control) sequence data. If not provided, a the generated control sequences will be written to runtime_tmp/", default=None)

    parser.add_option("-m", "--moo", dest="moo",
                      help="Multi-objective optimization [SPEA2, NSGA2, NSGAR, MOEAD]. NSGAR is the NSGA-II_R (NSGA-II Revised) algorithm improvement of NSGA2.", default="NSGAR")

    parser.add_option("-f", "--fitness", dest="fitness", type="str",
                      help="Objective fitness function. Available objectives: "+objectives_help+". Each single character in the string represents an objective. Objectives are mapped by the configuration file at config/objectives. Default is 'DF' for [Discrimination,Fisher] (2-objectives).", default="DF")

    parser.add_option("--cxpb", type="float",
                      dest="cxpb", default=0.7,
                      help="Probability [0.0 to 1.0] for a crossover during variation. Requires --mutpb to be set to (1.0-cxpb). Default is 0.7.")

    parser.add_option("--mutpb", type="float",
                      dest="mutpb", default=0.3,
                      help="Probability [0.0 to 1.0] for a mutation during variation. Requires --cxpb to be set to (1.0-mutpb). Default is 0.3.")

    parser.add_option("--short", type="int",
                      dest="short", default=0,
                      help="Stops reading in after <SHORT> input sequences.")

    parser.add_option("--popsize", type="int",
                      dest="popsize", default=1000,
                      help="Size of the population.")

    parser.add_option("--revcomp", action="store_true",
                      dest="revcomp", default=False,
                      help="Compile regex with reverse complement")

    parser.add_option("--random-seed", type="int",
                  dest="random_seed", default=42,
                      help="Random seed value to set for execution")

    parser.add_option("-n", "--num-gen", type="int",
                  dest="ngen", default=None,
                      help="Generation where runtime stops (even in the case of resumed checkpoints)")

    parser.add_option("--timelimit", type="int",
                  dest="timelimit", default=None,
                      help="Time limit on the GP loop execution.")

    if platform.system() == 'Windows':
        parser.add_option("--matcher", type="string",
                dest="matcher", default="python",
                help="Use a different matcher. Options: 'python', 'grep'. 'python' is a pure python version in case the system doesn't support grep. Install a grep as a system tool for compatibility.")
    else:
        parser.add_option("--matcher", type="string",
              dest="matcher", default="python",
              help="Use a different matcher. Options: 'grep', 'python'. 'grep' is faster on large datasets, while 'python' is a pure python version in case the system doesn't support grep.")

    parser.add_option("-o", "--output", type="string",
                  dest="output_path", default="./OUT/",
                      help="Output directory. Default is ./OUT/ ")

    parser.add_option("-t", "--tag", type="string",
                      dest="tag", default="default",
                      help="A tag for the output subdirectory. Use to describes the run and saves it in the tag's subdirectory in the output directory. default is 'default'.")

    parser.add_option("-i", "--inspector", action="store_true",
                      dest="inspector", default=False,
                      help="Don't print any files. Can be useful with python -i (interactive mode).")

    parser.add_option("--hardmask", action="store_true",
                      dest="hardmask", default=False,
                      help="Replace tandem repeats (lower-case typed nucleotides) by N")

    parser.add_option("-g", "--grammar", dest="grammar", type="str",
                      help="Grammar for the STGP [min, iupac, full, ne]. Default is iupac. 'min' only uses nucleotides. 'iupac' is a network expression grammar. 'full' is a network expression grammar with additional regular expression tokens. 'ne' is like iupac, but built with string primitives instead of booleans.", default="iupac")

    parser.add_option("-e", "--erase", dest="erase",
                      help="Input .nef(t) file to delete from the dataset prior to execution. Used for sequential coverage.", default=None)

    parser.add_option("--backpad", action="store_true",
                      dest="backpad", default=False,
                      help="Pads background sequences with consecutive nucleotides (ie. "+ ",".join([x*postprocessing.PADLENGTH for x in postprocessing.ALPHABET]) +") of length "+ str(postprocessing.PADLENGTH) +" every set of "+str(len(postprocessing.ALPHABET))+" sequences.")

    parser.add_option("--bg-algo", dest="bg_algo",
                      help="Shuffling algorithm for background. Default is 'dinuclShuffle', if no background dataset it provided. Currently, dinuclShuffle is the only implemented method.", default="dinuclShuffle")
    
    parser.add_option("--ncpu", dest="ncpu",
                      help="Number of CPUs to use when mapping evaluation of solutions. Use an integer, \"auto\" to automatically dertmine the maximum number. Default is no parallelism.",
                      default=None)

    parser.add_option("--termination", dest="termination",
                      help="Use automatic termination algorithm. User 'auto' to used the automatic termination algorithm for MOEAs.",
                      default=None)

    parser.add_option("--steps", dest="steps",
                      help="Grammar scafolding steps.",
                      default=None)



    # Currently not used

# Currently disabled, for perfomance
#    parser.add_option("--hypervolume", action="store_true",
#                    dest="hypervolume", default=False,
#                    help="[Experimental] Output a trace of the hypervolume.")

    parser.add_option("--hamming", action="store_true",
                    dest="hamming", default=False,
                    help="[Experimental] Generates statistics on the hamming distance from a template regex and hof candidates.")

    parser.add_option("--seeded-population", action="store_false",
                      dest="no_seed", default=True,
                      help="[Experimental] Use population seeds")

    parser.add_option("-c", "--checkpoint", type="string",
                      dest="checkpoint_path", default=None,
                      help="[Temporarily disabled] Load a checkpoint at path.")

    parser.add_option("-q", "--quiet",
                      action="store_false", dest="verbose", default=True,
                      help="[Unimplemented] don't print status messages to stdout")


    (options, args) = parser.parse_args()
    check_parameter_constraints(options)

    if options.inspector:
        OUTPUT_PATH = None #Don't print file
    else:
        OUTPUT_PATH = os.path.join(options.output_path, options.tag, options.moo, str(options.random_seed))
        if not os.path.exists(OUTPUT_PATH):
            os.makedirs(OUTPUT_PATH)

    print "Options:"
    print ast.literal_eval(options.__str__())
    #print [ x.get_opt_string() for x in parser._get_all_options()[1:]]

    random.seed(options.random_seed)
    numpy.random.seed(options.random_seed)
    toolbox = None
    population = None
    pareto_front = None
    hof = None
    mstats = None
    invalid_ind = None
    logbook = None
    pset = None
    rndstate = None
    numpystate= None

    EXPERIMENT="regex"

    if options.timelimit is None and \
       options.ngen is None:
        options.ngen=10
    elif options.timelimit is not None and\
         options.ngen is None:
        options.ngen=1 # TODO: Edit engine to support infinite generations


    NGEN=options.ngen



    start_gen = 0

    utils = Utils()

    engine = Engine(OUTPUT_PATH=OUTPUT_PATH, fitness=options.fitness)

    if options.checkpoint_path:
        engine.dna_pset(options.grammar) # FIXME: This gets overridden by boot()
        population, toolbox,start_gen,hof,logbook, rndstate, numpystate = evoalgo.explode_checkpoint(options.checkpoint_path, engine.pset)
        print "Checkpoint loaded"

    # BEGIN MAIN EXECUTION
    engine.boot(options)
    engine.set_experiment(EXPERIMENT, options, mstats, logbook, pset, toolbox, population, invalid_ind, hof, START_GEN=start_gen, POP_SIZE=options.popsize)

    if options.erase :
        print "Erasing motifs from", options.erase
        motifs = utils.get_motifs_from_nef(options.erase)
        print "Erasing", len(motifs), "from the dataset."
        engine.re_positive_dataset = utils.motif_eraser(engine.re_positive_dataset,
                                                        motifs)
    if options.hardmask:
        regex = ["[a-z]+"]
        engine.re_positive_dataset = utils.motif_eraser(engine.re_positive_dataset,
                                                        regex)

    if options.backpad:
        print "Padding the background"
        engine.re_negative_dataset = postprocessing.backpad(engine.re_negative_dataset)

    positive_count , negative_count = len(engine.re_positive_dataset), len(engine.re_negative_dataset)

    print "Input sequences#", positive_count
    print "Control sequences#", negative_count
    if positive_count != negative_count:
        print "WARNING: Uneven number of input and control sequences."

    print "Running algorithm:", options.moo
    print "Crossover/Mutations", options.cxpb, options.mutpb
    start_time = timeit.default_timer()

    if rndstate and numpystate:
        random.setstate(rndstate)
        numpy.random.set_state(numpystate)

    offspring, log = engine.run(options.ngen)
    elapsed = timeit.default_timer() - start_time
    # END MAIN EXECUTION

    if options.verbose:
        print "Overall runtime: ", elapsed, "time"
        of = None
        if (options.inspector):
            of = Inspector() # Write to STDOUT
            of.write_nef(engine)

        # Write Network Expression Front
        if OUTPUT_PATH:
            statcsv.write_nef(OUTPUT_PATH, engine)
            statcsv.write_log(OUTPUT_PATH, log)
            

    #if options.hypervolume:
    #   #statcsv.write_hypervolume(OUTPUT_PATH,log)
    #   #pass


    def template_sequence(headers, sequences):
        """
        Creates a list of pairs that are (template,sequence), from hamming benchmarking

        Assumes headers format is

             >$TEMPLATE_REGEX, $otherstuf1, $otherstuff2

        TODO: move away from encoding that stuff in the header.
        """
        template_sequence_list = []
        for h, s in zip(headers, sequences):
            header = ",".join(str(h[1:]).split(",")[:-2])
            #print header, s
            #raw_input()
            template_sequence_list.append( tuple((header, s))  )
        return template_sequence_list


    if options.hamming:
        hamming_benchmark = hammingregex.HammingBenchmark()

        fitnesses = [x[0] for x in hof_fitness_candidate_tuple]
        regexs = [x[1] for x in hof_fitness_candidate_tuple]


        sequence_tuples = template_sequence(engine.header_positive_dataset, engine.re_positive_dataset)
        hamming_benchmark.compile(regexs, sequence_tuples)

        print hamming_benchmark

        hammingout = os.path.split(options.checkpoint_path)[0]
        hammingout = hammingout.split("/")
        hammingout, idx = hammingout[:-1], hammingout[-1]
        hammingout = "/".join(hammingout)
        hammingout = os.path.join(hammingout, "hammingpoints")
        print "Writing hamming at", hammingout

        hamming_benchmark.flush_data_points(regexs, fitnesses, hammingout, idx)

    if options.inspector and options.verbose:
        print "Use interactive mode for items: hof, population, offspring, log"
