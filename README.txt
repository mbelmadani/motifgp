===============
= MotifGP 0.1 =
===============
MotifGP is a de novo motif discovery tool for discriminatory network expression identification in ChIP-seq datasets.

Original author: Manuel Belmadani
	mbelm006@uottawa.ca

Acknowledgements:
MotifGP is using source code from these tools:
-hypergeometric.py from the MEME Suite (License and copyright in source file).
-altschulEriksonDinuclShuffle.py from Peter Clote - CLOTE Computational Biology LAB, http://clavius.bc.edu/~clotelab/RNAdinucleotideShuffle/
This software was also made using the DEAP - Fortin, F.-A., De Rainville, F.-M., Gardner, M.-A. G., Parizeau, M. & Gagné, C. DEAP: Evolutionary Algorithms Made Easy. J. Mach. Learn. Res. 13, 2171–2175 (2012).

=======================================================================================
License: (see LICENSE.txt)
=======================================================================================
Installation: (see INSTALL.txt)
=======================================================================================
Examples: (see EXAMPLES.txt)
=======================================================================================
Usage: motifgp.py [options]

Options:
  -h, --help            show this help message and exit
  -p TRAINING_PATH, --training=TRAINING_PATH
                        Fasta file to use for training (input) sequence data
  -b BACKGROUND_PATH, --background=BACKGROUND_PATH
                        [Optional] Fasta file to use for background (control)
                        sequence data. If not provided, a the generated
                        control sequences will be written to runtime_tmp/
  -m MOO, --moo=MOO     Multi-objective optimization [SPEA2, NSGA2, NSGAR].
                        NSGAR is the NSGA-II_R (NSGA-II Revised) algorithm
                        improvement of NSGA2.
  -f FITNESS, --fitness=FITNESS
                        Objective fitness function. Available objectives: D=Di
                        scrimination,F=Fisher,I=ScipyFisher,O=OddsRatio,Q=Fals
                        eDiscoveryRate,S=Support,R=ScipyOddsRatio. Each single
                        character in the string represents an objective.
                        Objectives are mapped by the configuration file at
                        config/objectives. Default is 'DF' for
                        [Discrimination,Fisher] (2-objectives).
  --cxpb=CXPB           Probability [0.0 to 1.0] for a crossover during
                        variation. Requires --mutpb to be set to (1.0-cxpb).
                        Default is 0.7.
  --mutpb=MUTPB         Probability [0.0 to 1.0] for a mutation during
                        variation. Requires --cxpb to be set to (1.0-mutpb).
                        Default is 0.3.
  --short=SHORT         Stops reading in after <SHORT> input sequences.
  --popsize=POPSIZE     Size of the population.
  --revcomp             Compile regex with reverse complement
  --random-seed=RANDOM_SEED
                        Random seed value to set for execution
  -n NGEN, --num-gen=NGEN
                        Generation where runtime stops (even in the case of
                        resumed checkpoints)
  --timelimit=TIMELIMIT
                        Time limit on the GP loop execution.
  --matcher=MATCHER     Use a different matcher. Options: 'grep', 'python'.
                        'grep' is faster on large datasets, while 'python' is
                        a pure python version in case the system doesn't
                        support grep.
  -o OUTPUT_PATH, --output=OUTPUT_PATH
                        Output directory. Default is ./OUT/
  -t TAG, --tag=TAG     A tag for the output subdirectory. Use to describes
                        the run and saves it in the tag's subdirectory in the
                        output directory. default is 'default'.
  -i, --inspector       Don't print any files. Can be useful with python -i
                        (interactive mode).
  --hardmask            Replace tandem repeats (lower-case typed nucleotides)
                        by N
  -g GRAMMAR, --grammar=GRAMMAR
                        Grammar for the STGP [min, iupac, full, ne]. Default
                        is iupac. 'min' only uses nucleotides. 'iupac' is a
                        network expression grammar. 'full' is a network
                        expression grammar with additional regular expression
                        tokens. 'ne' is like iupac, but built with string
                        primitives instead of booleans.
  -e ERASE, --erase=ERASE
                        Input .nef(t) file to delete from the dataset prior to
                        execution. Used for sequential coverage.
  --bg-algo=BG_ALGO     Shuffling algorithm for background. Default is
                        'dinuclShuffle', if no background dataset it provided.
                        Currently, dinuclShuffle is the only implemented
                        method.
  --backpad             Pads background sequences with consecutive nucleotides
                        (ie. AAAAAAAA,CCCCCCCC,GGGGGGGG,TTTTTTTT) of length 8
                        every set of 4 sequences.
  --hamming             [Experimental] Generates statistics on the hamming
                        distance from a template regex and hof candidates.
  --seeded-population   [Experimental] Use population seeds
  -c CHECKPOINT_PATH, --checkpoint=CHECKPOINT_PATH
                        [Temporarily disabled] Load a checkpoint at path.
  -q, --quiet           [Unimplemented] don't print status messages to stdout

Also consider looking at EXAMPLES.txt for basic examples of MotifGP usage.
