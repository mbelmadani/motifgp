from nose.tools import assert_equals
from nose.tools import assert_is_none
from nose.tools import assert_is_not_none

from utils import Utils


def test_alphabet():
    u = Utils()
    assert_equals('A', u.complement['T'])
    assert_equals('C', u.complement['G'])
    assert_equals('G', u.complement['C'])
    assert_equals('T', u.complement['A'])
    assert_equals('[', u.complement[']'])
    assert_equals(']', u.complement['['])

def test_readfasta():
    """
    Tests that readfasta returns 1 header + 1 sequence continuously
    Test that HARDMASK leaves lower case while if True, else convert to upper()
    """
    
    filename = "tests/resources/sequences1.fasta"
    sequences = ["ACGTTTTTTGCA",
                 "ACGTAAAAGGGG",
                 "ACCCCCCAA",
                 "AAAAttttAAAA"
    ]

    results=[]
    r_headers=[]
    r_sequences=[]
    
    for fasta in Utils().readfasta(filename, HARDMASK=True):
        results.append(fasta)

    r_headers = results[::2]
    r_sequences = results[1::2]
        
    for s1, s2 in zip(r_sequences, sequences):
        """
        Test hardmasked
        """
        assert_equals(s1, s2)

    results=[]
    for fasta in Utils().readfasta(filename, HARDMASK=False):
        results.append(fasta)


    r_sequences = results[1::2]
    for s1, s2 in zip(r_sequences, sequences):
        """
        Test non-hardmasked
        """
        assert_equals(s1, s2.upper())
    

# UNTESTED
# def test_produce_sequences():

def test_list_to_regex():
    u = Utils()

    REGEX = "AAA[AC]TT[A]"
    STRING1 = "AAACTTA" #Pass
    STRING2 = "AAAZTTA" #Fail
    STRING3 = "AAAATTA" #Alternative

    # Test str regex to regex object
    re = u.list_to_regex_obj(REGEX)
    assert_is_not_none(re.match(STRING1))
    assert_is_none(re.match(STRING2))
    assert_is_not_none(re.match(STRING3))

    # Test string regex to token list
    re_list = u.motif_str_to_list(REGEX)
    expected_list = ["A", "A", "A", "AC", "T", "T", "A"]
    for token, expected in zip(re_list, expected_list):
        assert_equals(token, expected)
    
def test_reverse_completement():
    """
    Test reverse compelement functions
    """
    u = Utils()
    REGEX = "AAA[AC]TT[A]"
    RC = "[T]AA[GT]TTT"

    # Test conversion from regex to reverse complement
    assert_equals(u.reverse_complement(REGEX), RC)

    # Test appending reverse complement to regex
    assert_equals(u.add_reverse_complement(REGEX), str(REGEX + "|" + RC))
    

#UNTESTED
# def test_get_motifs_from_nef
