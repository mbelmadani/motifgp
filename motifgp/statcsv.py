import csv
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
        front = engine.hof
        for ind in front:
            pattern = engine.toolbox.compile(expr=ind)
            if len(u.motif_str_to_list(pattern)) < 5:
                continue # Skip motifs shorter than length 5

            score =  [float(x) for x in ind.fitness.getValues()]
            candidate = [count, pattern] + score
            count += 1
            writer.writerow(candidate)
    print "NEF File written at", nef_path
