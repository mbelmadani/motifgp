from __future__ import print_function

from nose.tools import assert_equals
from nose.tools import assert_almost_equals
from nose.tools import assert_is_none
from nose.tools import assert_is_not_none

import deap.creator
import deap.base

from stgpfitness import STGPFitness
from hypergeometric import *
from utils import Utils

class Object(object):
    pass

    
def skipTestScipy(func):
    try:
        import scipy        
    except:
        #self.skipTest('skipped test due to skip_tests_flag')
        import sys
        msg = "Skipping Scipy dependent test: "+ func.__name__        
        print(msg, file=sys.stderr)
        return
    func()

@skipTestScipy
def test_matchers():
    objectives = "FDSI"
    stgp = STGPFitness(objectives)

    regex = "AAAACGTAAAA"
    regex = "A"
    p_strings = ["TTTTTTTAAAACGTAAAATTTTTTT"]
    n_strings = ["TTTTTTTTTTTTTTTTTTTTTTTTT"]

    stgp.re_positive_dataset = p_strings
    stgp.re_negative_dataset = n_strings

    options = Object()
    options.training_path = "tests/resources/test_positive"
    options.background_path = "tests/resources/test_negative"

    stgp.options = options
    
    python_score = stgp.memoize_or_python_match(regex)
    python_score = [float(sprint_logx(python_score[0], 3, "%6.3fe%-5.0f"))] + [x for x in python_score[1:]]
    
    grep_score = stgp.memoize_or_grep(regex)
    grep_score = [float(sprint_logx(grep_score[0], 3, "%6.3fe%-5.0f"))] + [x for x in grep_score[1:]]
    assert_equals(python_score, [0.5, 1.0, 1.0, 0.5] )
    assert_equals(grep_score, [0.5, 1.0, 1.0, 0.5] )
    assert_equals(python_score, grep_score)

@skipTestScipy
def test_objectives():
    # TEST 1: Odds ratio and Fisher's exact test should work like the scipy version
    
    objectives = "FIOR"
    stgp = STGPFitness(objectives)

    regexs = ["AAAACGTAAAA",
              "AACTGAA",
              "AAA",
              "TTT",
              "GGG",
              "CCC",
              "AA[ACGT]CC[ACGT]"]

    
    options = Object()
    options.training_path = "tests/resources/sample_ctcf.fasta"
    options.background_path = None

    p_strings,  n_strings, _,_ = Utils().load_datasets(options.training_path, None, 'dinuclShuffle', 0, HEADERS=True)

    stgp.re_positive_dataset = p_strings
    stgp.re_negative_dataset = n_strings
    stgp.options = options

    for regex in regexs:
        python_score = stgp.memoize_or_python_match(regex) #, DEBUGGING=True)
        #print python_score
        #print( float(sprint_logx(python_score[0], 10, "%6.10fe%-5.0f")),
        #       python_score[1])


        assert_almost_equals(float(sprint_logx(python_score[0], 10, "%6.10fe%-5.0f")),
                                 python_score[1])
        assert_almost_equals(python_score[2], python_score[3])
