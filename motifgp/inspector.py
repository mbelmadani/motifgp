class Inspector():
    def __init__(self):
        pass

    def write(self, string):
        print string.rstrip()

    def write_nef(self, engine):
        self.write("Pareto Front: \n")
        hof_count = 0
        
        for ind in engine.hof:
            pattern = engine.toolbox.compile(expr=ind)
            hof_count+=1
            line = "\t".join([str(hof_count), pattern, str(ind.fitness)])
            self.write(line)
