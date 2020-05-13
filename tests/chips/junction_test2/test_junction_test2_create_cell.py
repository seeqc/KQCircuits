from tests.chips.chip_test_helpers import errors_test, base_refpoint_existence_test

from kqcircuits.chips.junction_test2 import JunctionTest2


def test_errors(capfd):
    errors_test(capfd, JunctionTest2)


def test_base_refpoint_existence():
    base_refpoint_existence_test(JunctionTest2)
