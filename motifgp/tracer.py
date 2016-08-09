from deap.benchmarks import tools

class Tracer():
    """
    Traces various aspect of the GP useful for benchmarking
    """
    def __init__(self, old_tracer=None):
        if old_tracer == None:
            self.convergence = []
            self.diversity = []
            self.generation = 0
        else:
            self.convergence = old_tracer.convergence
            self.diversity = old_tracer.diversity
            self.generation = old_tracer.generation

    def trace(self, old, new):
        """
        Takes parent hall of fame and offspring, evaluates metrics to trace.
        """
        #self.convergence.append(tools.convergence(old,new))
        self.diversity.append(tools.diversity(new, new[0], new[-1]))
        self.generation += 1
