import deap
from deap import tools
from deap.algorithms import varOr
from deap.algorithms import varAnd

from deap import creator
from deap import base
from deap import gp
from deap import cma
from deap import benchmarks

from collections import deque

#import dill as pickle
#import pickle 
import math
import json
import random
import numpy
import time

from termination.automoea import MOEATerminationDetection
from termination.timer import TimeTerminator
from termination.generation import GenerationTerminator

DEFAULT_CHECKPOINTFREQ = 9000000

def eaFortin(population, 
             toolbox, 
             mu, 
             cxpb,
             mutpb,
             ngen,
             stats,
             halloffame,
             logbook,
             invalid_ind,
             FREQ=DEFAULT_CHECKPOINTFREQ,
             start_gen=1,
             checkpoint=None,
             verbose=__debug__,
             CPOUT=None,
             pset=None,
             timelimit=None,
             termination=None
             ):
    """
    EA based of the sample in the fortin2013 repository
    """
    if timelimit:
        terminator = TimeTerminator(limit=timelimit)
    else:
        terminator = GenerationTerminator(limit=ngen)
    
    if checkpoint:
        # Removed feature
        print "Resuming from generation", start_gen    
        print logbook.stream
        
    else:
        logbook = tools.Logbook()
        logbook.header = ['gen', 'evals', 'memoize', "maxdepth"] + (stats.fields if stats else [])

        # Evaluate the individuals with an invalid fitness
        invalid_ind = [ind for ind in population if not ind.fitness.valid]
        fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)
        for ind, fit in zip(invalid_ind, fitnesses):
            ind.fitness.values = fit

        # This is just to assign the crowding distance to the individuals
        # no actual selection is done
        population = toolbox.select(population, len(population))

        if halloffame is not None:
            halloffame.update(population)
        #hypervolume_score =  benchmarks.tools.hypervolume(halloffame)

        record = stats.compile(population) if stats is not None else {}
        logbook.record(gen=start_gen, evals=len(invalid_ind), memoize=0, maxdepth=numpy.max([x.height for x in population]), **record)

    if verbose:
        print logbook.stream

    if termination == "auto":
        #TODO: Place somewhere better
        terminator = MOEATerminationDetection()
        
    gen = start_gen    
    while not terminator.done():
        # Vary the population
        offspring = toolbox.preselect(population, len(population))
        offspring = [toolbox.clone(ind) for ind in offspring]        
        offspring = varAnd(offspring, toolbox, cxpb, mutpb)
        
        # Evaluate the individuals with an invalid fitness
        invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
        fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)
        for ind, fit in zip(invalid_ind, fitnesses):
            ind.fitness.values = fit
        # Select the next generation population
        population = toolbox.select(population + offspring, mu)
        #print "\n".join([str(toolbox.compile(p)) + "|" + str(p) for p in population])
        #raw_input()
        

        gen += 1
        # Update termination conditions
        if termination == "auto":
            P = [list(y) for y in [x.fitness.values for x in halloffame[:]]]
            halloffame.update(population)
            Q = [list(y) for y in [x.fitness.values for x in halloffame[:]]]
            terminator.update(P,Q, gen-1) # Update diversity from last generation
        else:
            halloffame.update(population)
            terminator.update(gen=gen)
        
        record = stats.compile(population)
        logbook.record(gen=gen, evals=len(invalid_ind), memoize=toolbox.memoizecount(), maxdepth=numpy.max([x.height for x in population]), **record ) 

        if verbose: print(logbook.stream)

        if CPOUT:
            continue # Checkpoints disabled
            this_checkpoint = CPOUT+"/"+str(gen)+".checkpoint"
            check_store_checkpoint(this_checkpoint, 
                                   population, 
                                   gen, 
                                   halloffame,
                                   logbook,
                                   toolbox,
                                   rndstate=random.getstate(),
                                   numpystate=numpy.random.get_state(),
                                   pset=None                                   
                                   )

    if CPOUT:
        this_checkpoint = CPOUT+"/"+str(ngen)+".checkpoint"
        store_checkpoint(this_checkpoint, 
                         population, 
                         ngen, 
                         halloffame,
                         logbook,
                         toolbox,
                         rndstate=random.getstate(),
                         numpystate=numpy.random.get_state(),
                         pset=None
        )
    return population, logbook



def checkpoint_eaMuPlusLambda(population, 
                              toolbox, 
                              mu, 
                              lambda_, 
                              cxpb, 
                              mutpb, 
                              ngen,
                              logbook,
                              start_gen=0,
                              stats=None, 
                              halloffame=None, 
                              FREQ=DEFAULT_CHECKPOINTFREQ,
                              checkpoint=None,
                              CPOUT=None,
                              verbose=__debug__,
                              timelimit=None,
                              termination=None):
    """This is the :math:`(\mu + \lambda)` evolutionary algorithm.
    
    :param population: A list of individuals.
    :param toolbox: A :class:`~deap.base.Toolbox` that contains the evolution
                    operators.
    :param mu: The number of individuals to select for the next generation.
    :param lambda\_: The number of children to produce at each generation.
    :param cxpb: The probability that an offspring is produced by crossover.
    :param mutpb: The probability that an offspring is produced by mutation.
    :param ngen: The number of generation.
    :param stats: A :class:`~deap.tools.Statistics` object that is updated
                  inplace, optional.
    :param halloffame: A :class:`~deap.tools.HallOfFame` object that will
                       contain the best individuals, optional.
    :param verbose: Whether or not to log the statistics.
    :returns: The final population
    :returns: A class:`~deap.tools.Logbook` with the statistics of the
              evolution.
    
    The algorithm takes in a population and evolves it in place using the
    :func:`varOr` function. It returns the optimized population and a
    :class:`~deap.tools.Logbook` with the statistics of the evolution. The
    logbook will contain the generation number, the number of evalutions for
    each generation and the statistics if a :class:`~deap.tools.Statistics` is
    given as argument. The *cxpb* and *mutpb* arguments are passed to the
    :func:`varOr` function. The pseudocode goes as follow ::

        evaluate(population)
        for g in range(ngen):
            offspring = varOr(population, toolbox, lambda_, cxpb, mutpb)
            evaluate(offspring)
            population = select(population + offspring, mu)

    First, the individuals having an invalid fitness are evaluated. Second,
    the evolutionary loop begins by producing *lambda_* offspring from the
    population, the offspring are generated by the :func:`varOr` function. The
    offspring are then evaluated and the next generation population is
    selected from both the offspring **and** the population. Finally, when
    *ngen* generations are done, the algorithm returns a tuple with the final
    population and a :class:`~deap.tools.Logbook` of the evolution.

    This function expects :meth:`toolbox.mate`, :meth:`toolbox.mutate`,
    :meth:`toolbox.select` and :meth:`toolbox.evaluate` aliases to be
    registered in the toolbox. This algorithm uses the :func:`varOr`
    variation.
    """

    if timelimit:
        terminator = TimeTerminator(limit=timelimit)
    else:
        terminator = GenerationTerminator(limit=ngen)
        
    if checkpoint:
        # Removed feature
        print "Resuming from generation", start_gen
        print logbook.stream
    else:
        
        logbook = tools.Logbook()
        logbook.header = ['gen', 'evals', 'memoize', "maxdepth"] + (stats.fields if stats else [])

        # Evaluate the individuals with an invalid fitness
        invalid_ind = [ind for ind in population if not ind.fitness.valid]
        fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)
        
        for ind, fit in zip(invalid_ind, fitnesses):
            ind.fitness.values = fit

        if halloffame is not None:
            halloffame.update(population)

        record = stats.compile(population) if stats is not None else {}        
        logbook.record(gen=start_gen, evals=len(invalid_ind), memoize=0, maxdepth=numpy.max([x.height for x in population]), **record)
        
    if verbose:
        print logbook.stream

    if len(population) == 1:
        population.append(population[0])

    
    if termination == "auto":
        #TODO: Place somewhere better
        terminator = MOEATerminationDetection()

    # Begin the generational process
    gen = start_gen
    while not terminator.done():
        # Vary the population        
        offspring = varOr(population, toolbox, lambda_, cxpb, mutpb)
        
        # Evaluate the individuals with an invalid fitness
        invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
        fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)
        for ind, fit in zip(invalid_ind, fitnesses):
            ind.fitness.values = fit
            
        # Update the hall of fame with the generated individuals
        # Update termination conditions
        gen+=1
        if termination == "auto":
            P = [list(y) for y in [x.fitness.values for x in halloffame[:]]]
            halloffame.update(offspring)
            Q = [list(y) for y in [x.fitness.values for x in halloffame[:]]]
            terminator.update(P,Q, gen-1)
        else:
            halloffame.update(offspring)
            terminator.update(gen=gen)
        
        # Select the next generation population
        population[:] = toolbox.select(population + offspring, mu)
            
        record = stats.compile(population) if stats is not None else {}

        logbook.record(gen=gen, evals=len(invalid_ind), memoize=toolbox.memoizecount(), maxdepth=numpy.max([x.height for x in population]), **record ) 

        if verbose: print logbook.stream

        # Store checkpoint if gen is valid
        if CPOUT:
            continue # Checkpoints disabled
            this_checkpoint = CPOUT+"/"+str(gen)+".checkpoint"
            check_store_checkpoint(this_checkpoint, 
                             population, 
                             gen, 
                             halloffame,
                             logbook,
                             toolbox,
                             rndstate=random.getstate(),
                             numpystate=numpy.random.get_state(),
                             pset=None
            )

                
    if CPOUT :
        this_checkpoint = CPOUT+"/"+str(ngen)+".checkpoint"
        store_checkpoint(this_checkpoint, 
                             population, 
                             ngen, 
                             halloffame,
                             logbook,
                             toolbox,
                             rndstate=random.getstate(),
                             numpystate=numpy.random.get_state(),
                             pset=None
            )
    return population, logbook





def check_store_checkpoint(path,
                           population,
                           generation,
                           halloffame,
                           logbook,
                           toolbox,
                           rndstate=random.getstate(),
                           numpystate=numpy.random.get_state(),
                           pset=None
                           ):
    """
#    Wrapper method to check if checkpoint should be stored and does it if it's valid
    """
    if path != None and generation%DEFAULT_CHECKPOINTFREQ == 0:
        store_checkpoint(path, 
                         population, 
                         generation, 
                         halloffame, 
                         logbook,
                         toolbox, 
                         rndstate=rndstate,
                         numpystate=numpystate,
                         pset=pset)

def store_checkpoint(path,
                     population,
                     generation,
                     halloffame,
                     logbook,
                     toolbox,
                     rndstate,
                     numpystate,
                     pset=None
                     ):
    """
#    Stores the curent population and stats along with the 'random' state
    """
    cp = dict(population=population,
              toolbox=toolbox,
              generation=int(generation),
              halloffame=halloffame,
              logbook=logbook,
              rndstate=rndstate,
              numpystate=numpystate
              )

    if True:
        return
    
    print "Writing gen", generation, "at path", path
    with open(path, 'w') as f:
        json.dump(cp, f)
        pickle.dump(cp, f)

def explode_checkpoint(path, pset):
    #print "run"
    #creator.create("FitnessMax", base.Fitness, weights=w)
    #creator.create("Individual", gp.PrimitiveTree, fitness=creator.FitnessMax, pset=pset)
    #print "time"

    try:
        deap.creator.Individual
    except:
        deap.creator.create("Individual", deap.gp.PrimitiveTree, fitness=deap.creator.FitnessMax, pset=pset)

    with open(path, "r") as f:
        #cp = pickle.load(f)
        cp = json.load(f)
        
    #return cp['population'], cp['generation'], cp['halloffame'], cp['logbook'], cp['toolbox'], cp['rndstate'], cp['numpystate']
    return cp['population'], cp['toolbox'], cp['generation'], cp['halloffame'], cp['logbook'], cp['rndstate'], cp['numpystate']


"""
def check_store_checkpoint(path,
                           population,
                           generation,
                           halloffame,
                           logbook,
                           toolbox,
                           rndstate=random.getstate(),
                           numpystate=numpy.random.get_state(),
                           pset=None
                           ):
"""
#    Wrapper method to check if checkpoint should be stored and does it if it's valid
"""
    if path != None and generation%DEFAULT_CHECKPOINTFREQ == 0:
        store_checkpoint(path, 
                         population, 
                         generation, 
                         halloffame, 
                         logbook,
                         toolbox, 
                         rndstate=rndstate,
                         numpystate=numpystate,
                         pset=pset)

def store_checkpoint(path,
                     population,
                     generation,
                     halloffame,
                     logbook,
                     toolbox,
                     rndstate,
                     numpystate,
                     pset=None
                     ):
"""
#    Stores the curent population and stats along with the 'random' state
"""
    cp = dict(population=population,
              toolbox=toolbox,
              generation=int(generation),
              halloffame=halloffame,
              logbook=logbook,
              rndstate=rndstate,
              numpystate=numpystate
              )

    print "Writing gen", generation, "at path", path
    if True:
        return
    with open(path, 'w') as f:
        json.dump(cp, f)
        pickle.dump(cp, f)

def explode_checkpoint(path, pset):
    #w = (1.0, 1.0, ) #TODO: Wrong
    #print "run"
    #creator.create("FitnessMax", base.Fitness, weights=w)
    #creator.create("Individual", gp.PrimitiveTree, fitness=creator.FitnessMax, pset=pset)
    #print "time"

    try:
        deap.creator.Individual
    except:
        deap.creator.create("Individual", deap.gp.PrimitiveTree, fitness=deap.creator.FitnessMax, pset=pset)

    with open(path, "r") as f:
        #cp = pickle.load(f)
        cp = json.load(f)
        
    #return cp['population'], cp['generation'], cp['halloffame'], cp['logbook'], cp['toolbox'], cp['rndstate'], cp['numpystate']
    return cp['population'], cp['toolbox'], cp['generation'], cp['halloffame'], cp['logbook'], cp['rndstate'], cp['numpystate']
"""
