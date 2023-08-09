import unittest
import csv
from ast import literal_eval

from .epipolymorphism import get_quartets, load_cpgs


with open('scratch/test.tsv') as f:
    mock_read_table = csv.DictReader(f, delimiter='\t')
    mock_read_hash = {}
    for row in mock_read_table:
        mock_read_hash[row['key']] = {
            'query_name': row['query_name'],
            'reference_start': int(row['reference_start']),
            'reference_end': int(row['reference_end']),
            'cigartuples': literal_eval(row['cigartuples']),
            'xm_tag': row['xm_tag'],
            'is_forward': row['is_forward'] == 'F',
            'reference_name': row['reference_name']
        }

with open('scratch/desired.tsv') as f:
    mock_quartet_table = csv.DictReader(f, delimiter='\t')
    mock_quartet_hash = {}
    for row in mock_quartet_table:
        mock_quartet = tuple(
            int(row['cpg%d_0' % (i+1)]) for i in range(4)
        ) + tuple(
            bool(int(row['meth%d' % (i+1)])) for i in range(4)
        )
        if row['key'] in mock_quartet_hash:
            mock_quartet_hash[row['key']].append(mock_quartet)
        else:
            mock_quartet_hash[row['key']] = [mock_quartet]
    mock_quartet_hash['no_quartets'] = []


class MockRead:
    def __init__(self, key):
        mock_read = mock_read_hash[key]
        self.reference_start = mock_read['reference_start']
        self.reference_end = mock_read['reference_end']
        self.xm_tag = mock_read['xm_tag']
        self.cigartuples = mock_read['cigartuples']
        self.reference_name = mock_read['reference_name']
        self.is_forward = mock_read['is_forward']

    def get_tag(self, tag='XM'):
        return self.xm_tag


class TestLoadCpgs(unittest.TestCase):
    def test_load(self):
        input_cpgs = 'scratch/bismark_CpG.txt'
        output_cpgs = 'scratch/cpgs.json'
        load_cpgs(input_cpgs, output_cpgs)


class TestGetQuartets(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cpgs = load_cpgs('scratch/bismark_CpG.txt')

    def generic_compare(self, key):
        read = MockRead(key)        
        quartets = get_quartets(read, self.cpgs)
        expected_quartets = mock_quartet_hash[key]
        self.assertEqual(quartets, expected_quartets)

    def test_simple(self):
        self.generic_compare('simple')

    def test_insertion(self):
        self.generic_compare('insertion')

    def test_deletion(self):
        self.generic_compare('deletion')

    def test_no_quartets(self):
        self.generic_compare('no_quartets')

    def test_cpg_off_by_one(self):
        self.generic_compare('cpg_off_by_one')

