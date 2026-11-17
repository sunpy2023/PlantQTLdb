import copy
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from validate_benchmark import index_rows, read_tsv, validate_scores, TOTAL, SCORE_COLUMNS


class ValidationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.answers = index_rows(read_tsv(ROOT / 'responses/primary_responses.tsv'),
                                 ('Response ID',), 'test responses')
        cls.score = read_tsv(ROOT / 'scores/primary_scores.tsv')[0]

    def test_valid_record(self):
        validate_scores([self.score], self.answers, 'test')

    def test_reject_incorrect_total(self):
        bad = copy.deepcopy(self.score)
        bad[TOTAL] = str(int(bad[TOTAL]) + 1)
        with self.assertRaisesRegex(ValueError, 'sum of five'):
            validate_scores([bad], self.answers, 'test')

    def test_reject_score_out_of_range(self):
        bad = copy.deepcopy(self.score)
        bad[SCORE_COLUMNS[0]] = '21'
        with self.assertRaisesRegex(ValueError, 'invalid dimension'):
            validate_scores([bad], self.answers, 'test')

    def test_reject_response_metadata_mismatch(self):
        bad = copy.deepcopy(self.score)
        bad['Question ID'] = 'unrelated-question'
        with self.assertRaisesRegex(ValueError, 'metadata mismatch'):
            validate_scores([bad], self.answers, 'test')

    def test_reject_duplicate_evaluations(self):
        with self.assertRaisesRegex(ValueError, 'duplicate key'):
            validate_scores([self.score, self.score], self.answers, 'test')

    def test_reject_unknown_response(self):
        bad = copy.deepcopy(self.score)
        bad['Response ID'] = 'missing-answer'
        with self.assertRaisesRegex(ValueError, 'unknown response'):
            validate_scores([bad], self.answers, 'test')


if __name__ == '__main__':
    unittest.main()
