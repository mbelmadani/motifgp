#!/usr/bin/python
import glob
import os
from utils import Utils

utils = Utils() 

def compute_sequence_pssm(sequence_as_tokens):
    """
    Takes a list of sequences as tokens and computes the weight of the PSSM for that single sequence
    """
    stack = []
    for token in sequence_as_tokens:
        weighted_slice = [0,0,0,0]
        divisor = float(len(token))
        def v(d=divisor):
            return 1.0/d

        if 'A' in token: weighted_slice[0] = v()
        if 'C' in token: weighted_slice[1] = v()
        if 'G' in token: weighted_slice[2] = v()
        if 'T' in token: weighted_slice[3] = v()
        if 'U' in token: weighted_slice[3] = v()
        #if sum(weighted_slice) > 0:
        #    continue
        """ Convert the IUPAC alphabet     A    C   G   T  """
        if 'R' in token: weighted_slice = [v(2), 0, v(2), 0] #[AG]
        if 'Y' in token: weighted_slice = [0, v(2), 0, v(2)] #[CT]
        if 'S' in token: weighted_slice = [0, v(2), v(2), 0] #[CG]
        if 'W' in token: weighted_slice = [v(2), 0, 0, v(2)] #[AT]
        if 'K' in token: weighted_slice = [0, 0, v(2), v(2)] #[GT]
        if 'M' in token: weighted_slice = [v(2), v(2), 0, 0] #[AC]
        
        if 'B' in token: weighted_slice = [0, v(3), v(3), v(3)] #[CGT]
        if 'D' in token: weighted_slice = [v(3), 0, v(3), v(3)] #[AGT]
        if 'H' in token: weighted_slice = [v(3), v(3), 0, v(3)] #[ACT]
        if 'V' in token: weighted_slice = [v(3), v(3), v(3), 0] #[ACG]

        if 'N' in token: weighted_slice = [v(4), v(4), v(4), v(4)] #[ACGT]

        stack.append(weighted_slice)
    return stack

def compute_multiple_sequence_pssm(a_list):
    """
    Expects a list of sequences, returns a 2D list representation of a PSSM
    """
    print "INPUT", a_list
    divisor = float(len(a_list))
    stack = None
    for e in a_list:
        """ Process each sequence """
        seq_pssm = compute_sequence_pssm(e)
        if stack ==None: 
            stack = seq_pssm
        else:
            for idx in range(len(seq_pssm)):
                temp = [stack[idx], seq_pssm[idx]]
                stack[idx] = [sum(i) for i in zip(*temp)]
            pass
    for idx in range(len(stack)):
        """ Divide by the number of input sequences """
        stack[idx] = [x/divisor for x in stack[idx]]        
    return stack

def write_MEME_PSSM(list_motif, output_path=None):
    PSSM = compute_multiple_sequence_pssm([list_motif])
    OUTPUT = ""
    OUTPUT += "MEME version 4.4"+"\n"
    OUTPUT += "\n"
    OUTPUT += "ALPHABET= ACGT"+"\n"
    OUTPUT += "\n"
    OUTPUT += "strands: + -"+"\n"
    OUTPUT += "\n"
    OUTPUT += "Background letter frequencies (from web form):"+"\n"
    OUTPUT += "A 0.25000 C 0.25000 G 0.25000 T 0.25000"+"\n"
    OUTPUT += "\n"
    OUTPUT += "MOTIF " + '['+str("][".join(list_motif))+']'  + "\n"
    OUTPUT += "\n"
    OUTPUT += "letter-probability matrix: alength= 4 w= "+str(len(list_motif))+" nsites= 20 E= 0"
    OUTPUT += "\n"
    for position in PSSM:
        OUTPUT += "\t  ".join([str(i) for i in position])
        OUTPUT += "\n"
    if output_path != None:
        print "Writing to", output_path
        with open(output_path, "w") as text_file:
            text_file.write("{0}".format(OUTPUT))
        return output_path
    if output_path == None:
        print OUTPUT


def concat_tsv(a, b, delimiter='\t'):
    a_s = a.split('\n')
    b_s = b.split('\n')
    r_s = []
    for i in range(len(a_s)):
        # Merge a_s and b_s to r_s
        r_s.append( 
            str(
                '\t'.join(a[i].split('\t') + b[i].split('\t')) 
                ) 
            )
    r = '\n'.join(r_s)
    return r


def NEF_to_MEME(nef):
    """ 
    Converts a Network Expressions for a Network Expression Pareto Front file (.nef)
    to a individual MEME files
    """
    # CREATE OUTPUT DIRECTORIES
    root_dir = os.path.dirname(os.path.abspath(nef))
    filename = os.path.basename(os.path.abspath(nef))
    bin = filename.split(".")[0] # The front ID
    output_path = os.path.join(root_dir, bin)
    if not os.path.isdir(output_path):
        os.mkdir(output_path)

    # Create MEME files for each motifs in the NEF file
    # Stored in each "bin" directory (0.nef MEMEs are in ./0/*.MEME)
    results = []
            
    # Clean older files
    oldfiles = glob.glob(os.path.join(output_path, "*.MEME"))
    print "Deleting", len(oldfiles), "stale MEME files"
    for oldfile in oldfiles:
        os.remove(oldfile)

    with open(nef) as f:
        next(f)
        for tsv_row in f:
            row = tsv_row.split("\t")
            counter = row[0]
            motif_str = row[1]
            # Individual output files for MEMEs
            MEME_output_path = os.path.join(output_path, (str(counter)+".MEME"))

            motif = utils.motif_str_to_list(motif_str)
            result = write_MEME_PSSM(motif, output_path=MEME_output_path)
            results.append(result)
    return results


DEFAULT_FIT = [0.5, 0.0]
def fasta_to_MEME(fasta, minFitnesses=DEFAULT_FIT):
    """
    WARNING: DEPRECATED; motifs are now in .nef format
    Converts a Fasta file to a individual MEME files
    """

    root_dir = os.path.dirname(os.path.abspath(fasta))
    filename = os.path.basename(os.path.abspath(fasta))
    bin = filename[0] # The front ID
    output_path = os.path.join(root_dir, bin)
    if not os.path.isdir(output_path):
        os.mkdir(output_path)

    counter = 0
    while os.path.isfile(os.path.join(root_dir, (str(counter)+".MEME"))):
        counter+=1
    with open(fasta, "r") as sequences:
        skipFlag = -1 # -1: Before optimal pareto. 0: Inside optimal pareto 1: Passed optimal pareto
        for line in sequences:
            if line[0] == ">":
                # Parse (**** [ACCURACY, COVERAGE] ******) from the output, converts to float, strips spaces and quotes.
                accuracy, coverage = [float(x.strip().replace("'","")) for x in line.split("[")[1].split("]")[0].split(",")]
                if minFitnesses == None:
                    skipFlag = 0 # Will read everthing
                else:
                    def isFit(value, minfit):
                        if True:
                            return True #FIXME: Currently we actually want to see how the hypergeometric test performs, so let's return True for all candidates.
                        if value >= minfit:
                            return True
                        return False
                    if isFit(accuracy, minFitnesses[0]) and isFit(coverage, minFitnesses[1]):
                        skipFlag = 0
                    else:
                        skipFlag = -1
            elif skipFlag == 1:
                break
            elif skipFlag == 0:
                # Parse sequence according to skipFlag
                current_output_path = os.path.join(output_path, (str(counter)+".MEME"))
                counter += 1
                motif = utils.motif_str_to_list(line)
                #print line
                #print motif
                #raw_input()
                write_MEME_PSSM(motif, output_path=current_output_path)
                skipFlag = -1

