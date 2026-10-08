"""Small exact and mutation tests, independent of production numerical outcomes."""
import unittest
from fractions import Fraction as F
from .experiment import (endpoints, path_record, check_witness, direct_witness,
                         membership_controls, recovered_moments, interval_point,
                         raw_synthesis)
from research.five_synthesis.model import Problem, Waypoint


class TheoryAlignmentTests(unittest.TestCase):
    def test_endpoints_and_small_mixtures(self):
        for name,end in endpoints():
            for eps in (F(0),F(1,10**10)):
                row=path_record(name,end,eps)
                self.assertTrue(check_witness(tuple(map(F,row['point'])),
                    [tuple(map(F,v['values'])) for v in row['sides']]))

    def test_bad_output_witness_rejected(self):
        point=endpoints()[0][1]
        values=direct_witness(point)
        values[0]=list(values[0]);values[0][5]=F(1)
        with self.assertRaises(AssertionError):
            check_witness(point,values)

    def test_bad_psd_witness_rejected(self):
        point=endpoints()[0][1]
        values=direct_witness(point)
        values[0]=list(values[0]);values[0][6]=F(-1)
        with self.assertRaises(AssertionError):
            check_witness(point,values)

    def test_exact_outside_separation(self):
        rows=list(membership_controls())
        self.assertEqual(len(rows),18)
        for row in rows:
            if row['kind']=='outside':
                self.assertEqual(F(row['exact_violation']),F(row['eta']))

    def test_control_moment_roundtrip(self):
        case=Problem('constant',F(1),(F(0),)*5,
                     (Waypoint(F(1),(F(0),)*4+(F(1),),F(2)),))
        record=dict(mesh=['0','1/3','1'],controls=['2','2'],upper='2')
        got=recovered_moments(case,record)
        self.assertEqual(got,[[F(2,j+1) for j in range(5)]])
        point=tuple(v/2 for v in got[0])
        self.assertEqual(point,interval_point(((0,1),)))
        self.assertTrue(check_witness(point,direct_witness(point)))

    def test_direct_convolution_detects_bad_raw_moments(self):
        case=Problem('zero',F(1),(F(0),)*5,
                     (Waypoint(F(1),(F(0),)*4+(F(1),),F(0)),))
        proposal=dict(value=0.,moments=[[0.]*5],auxiliaries=[[0.]*34])
        self.assertTrue(raw_synthesis(case,proposal)['diagnostic_pass'])
        proposal['moments'][0][0]=.1
        self.assertGreater(raw_synthesis(case,proposal)['maximum_scaled_row_error'],.19)
        self.assertFalse(raw_synthesis(case,proposal)['diagnostic_pass'])


if __name__=='__main__':
    unittest.main()
