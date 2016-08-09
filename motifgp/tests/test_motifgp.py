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
