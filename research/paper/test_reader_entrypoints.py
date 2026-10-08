"""Small adversarial checks for the reader interfaces, not more research cases."""
from pathlib import Path
from copy import deepcopy
import json
import os
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
from research.paper import check_lift_coefficients as algebra
from research.paper import reproduce as reader
from research.five_synthesis.check import check as check_plan
from research.five_feedback.check import replay

SOURCE=Path(__file__).with_name('representability.tex').read_text(encoding='utf-8')


def optimized_reader_process(script, cwd):
    """Use one explicit encoding on both sides of the subprocess boundary."""
    env = dict(os.environ, PYTHONIOENCODING='utf-8')
    return subprocess.run([sys.executable, '-O', '-B', '-S', str(script)],
                          capture_output=True, text=True, encoding='utf-8',
                          errors='strict', cwd=cwd, env=env)


class ReaderChecks(unittest.TestCase):
    def test_actual_pencils(self):
        result=algebra.check_pencils(SOURCE)
        self.assertEqual(result['functional_rank'],19)
        self.assertEqual(result['matrix_entries_checked'],67)

    def test_bad_pencil_rejected(self):
        bad=SOURCE.replace(r'180(\theta_{5}-\theta_{6}-\theta_{11})',
                           r'181(\theta_{5}-\theta_{6}-\theta_{11})')
        self.assertNotEqual(bad,SOURCE)
        with self.assertRaises(AssertionError):algebra.check_pencils(bad)

    def test_missing_pencil_rejected(self):
        bad=SOURCE.replace(r'M_{7}(\theta)&=',r'M_{8}(\theta)&=')
        with self.assertRaises(AssertionError):algebra.parse_pencils(bad)

    def test_non_arithmetic_latex_rejected(self):
        bad=SOURCE.replace(r'\gamma_{1}&=\theta_{3}',r'\gamma_{1}&=\input{secret}')
        with self.assertRaises(AssertionError):algebra.parse_pencils(bad)

    def test_remainder(self):
        self.assertTrue(algebra.check_remainder())

    def test_archive_anchor_mutation_rejected(self):
        with patch.dict(reader.ANCHORS,{'five_synthesis':'0'*64}):
            with self.assertRaisesRegex(ValueError,'anchor mismatch'):
                reader.verified_archive('five_synthesis')

    def test_path_escape_rejected(self):
        with self.assertRaisesRegex(ValueError,'outside'):
            reader.inside(reader.RESEARCH,r'..\outside.json')

    def test_windows_path_normalization(self):
        self.assertEqual(reader.inside(reader.RESEARCH,r'five_synthesis\check.py'),
                         (reader.RESEARCH/'five_synthesis/check.py').resolve())

    def test_exact_counts_not_tolerant(self):
        with self.assertRaises(ValueError):reader.same_summary({'count':16},{'count':15})
        with self.assertRaises(ValueError):reader.same_summary(True,1)

    def test_nonfinite_display_rejected(self):
        with self.assertRaises(ValueError):reader.same_summary(float('nan'),1.0)

    def test_bad_primal_certificate_rejected(self):
        record=reader.read(reader.RESEARCH/'five_synthesis/results/trials/00.json')
        plan=deepcopy(record['recovery']['certificate']);plan['upper']='1'
        with self.assertRaises(AssertionError):check_plan(plan)

    def test_late_activation_mutation_rejected(self):
        episode=reader.read(reader.RESEARCH/'five_feedback/results/episodes/03.json')
        row=episode['steps'][0]
        self.assertEqual(row['decision'],'late_discard')
        row['request']['accepted']=True
        with self.assertRaises(AssertionError):replay(episode)

    def test_replay_refuses_optimized_assertions(self):
        result=optimized_reader_process(Path(reader.__file__), reader.ROOT)
        self.assertNotEqual(result.returncode,0)
        self.assertIn('Do not run replay with -O',result.stderr)

    def test_replay_refusal_with_unicode_path_and_inherited_encodings(self):
        with tempfile.TemporaryDirectory(prefix='reader-encoding-') as temporary:
            directory = Path(temporary) / '中文路径'
            directory.mkdir()
            script = directory / '中文入口.py'
            script.write_text('import runpy\nrunpy.run_path(' +
                              repr(str(Path(reader.__file__))) +
                              ", run_name='__main__')\n", encoding='utf-8')
            for inherited in ('utf-8', 'gbk'):
                with self.subTest(inherited=inherited), patch.dict(
                        os.environ, {'PYTHONIOENCODING': inherited, 'PYTHONUTF8': '0'}):
                    result = optimized_reader_process(script, directory)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn('Do not run replay with -O', result.stderr)
                self.assertIn('中文入口.py', result.stderr)
                self.assertNotIn('\ufffd', result.stderr)


if __name__=='__main__':unittest.main()
