import csv
import sys
import utils

def write_hypervolume(path, log):
    logs = log.select("gen", "hypervolume")
    with open( path + "_hypervolume.csv" , 'w') as csv_out:
        statswriter = csv.writer(csv_out, delimiter='\t')
        statswriter.writerow(["gen", "hypervolume"])
        for i in range(len(logs[0])):
            statswriter.writerow([logs[0][i], logs[1][i]])
        
def write_nef(path, engine=None):
    tags = engine.FITNESS_TAGS
    header = ["Rank", "NE"] + tags
    u = utils.Utils()
    nef_path = path + ".nef"
    with open( nef_path, 'w') as nef_out:
        writer = csv.writer(nef_out, delimiter='\t')
        writer.writerow(header)


        count = 0
        _front = engine.hof
        front = []
        counter=0
        uniquepatterns = []
        for ind in _front:
            pattern = str(engine.toolbox.compile(expr=ind))
            if pattern in uniquepatterns:
                continue
            else:
                uniquepatterns.append(pattern)            
            score =  [str(float(x)) for x in ind.fitness.getValues()]
            candidate = [str(count), pattern] + score            
            front += [candidate]            
            count += 1

        for candidate in front:            
            writer.writerow(candidate)
            
    print ".nef file written at", nef_path


    if engine.REVCOMP:
        nefrc_path = nef_path+"rc"
        with open( nefrc_path, 'w') as nef_out:
            writer_rc = csv.writer(nef_out, delimiter='\t')
            writer_rc.writerow(header)
            count = 0

            front = engine.hof
            for compiled in uniquepatterns:
                if type(compiled) != str: #FIXME
                    pattern = compiled.spacer.insertInto(compiled.NE)
                else:
                    pattern  = compiled
                    
                if len(pattern) > 1:
                    pattern = u.add_reverse_complement(pattern)
                
                score =  [float(x) for x in ind.fitness.getValues()]
                candidate = [count, pattern] + score
                count += 1
                writer_rc.writerow(candidate)

        print ".nef file w/ revcomp written at", nefrc_path

def write_log(path, log):
    stats_path = path + ".stats"
    stdout = sys.stdout
    sys.stdout=open(stats_path,"w")
    print log
    sys.stdout = stdout
    print "Statistics written at", stats_path
    
