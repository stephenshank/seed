import unittest

from .epipolymorphism import get_quartets, load_cpgs


class MockRead:
    def __init__(self, reference_start, xm_tag, cigartuples, reference_name,
            reference_end):
        self.reference_start = reference_start
        self.reference_end = reference_end
        self.xm_tag = xm_tag
        self.cigartuples = cigartuples
        self.reference_name = reference_name

    def get_tag(self, tag='XM'):
        return self.xm_tag


class TestLoadCpgs(unittest.TestCase):
    def test_load(self):
        input_cpgs = 'scratch/bismark_CpG.txt'
        output_cpgs = 'scratch/cpgs.json'
        load_cpgs(input_cpgs, output_cpgs)


class TestGetQuartets(unittest.TestCase):
    def test_simple(self):
        # D00796:32:C8BP4ANXX:6:2214:2710:17885_1:N:0:GGCTAC
        read = MockRead(
            reference_start=10588,
            reference_end=10639,
            xm_tag='Z....x....x.........Z.h..h.xZ.xZ....h.....Z.Z.xZ.Z.',
            cigartuples=[(0, 51)],
            reference_name='chr1'
        )
        cpgs = {
            'chr1': [
                5, 10588, 10608, 10616, 10619, 10630, 10632, 10635, 10637, 16243
            ]
        }
        self.assertEqual(len(read.get_tag('XM')), 51)

