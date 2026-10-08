"""Exact regression checks for the October 2 proof clarifications.

These are internal algebra/transcription checks, not a formal proof of the
analytic theorems or global novelty.
No optimization or modification of archived records is performed.
"""
from pathlib import Path
import hashlib
import json
import re
import unittest

import sympy as s

from research.paper.check_lift_coefficients import polynomial_data, parse_pencils
from research.four_lift_prior_audit import algebra as generic
from research.moment_four_lift.exact import matrices as small_matrices

ROOT = Path(__file__).resolve().parents[2]
SOURCE = Path(__file__).with_name("representability.tex").read_text(encoding="utf-8")


def interval_data():
    x, delta, n, p, f, basis, features = polynomial_data()
    m, z, q, r = x
    alpha = s.Symbol("alpha")
    moments = tuple(s.expand(((alpha+m)**(j+1)-alpha**(j+1))/(j+1))
                    for j in range(5))
    return x, dict(zip(x, moments[:4])), moments[4]


def source_witness(source):
    raw = source.split(r"\theta^*={}&(", 1)[1].split(").", 1)[0]
    raw = raw.replace(r"\hspace{24mm}", "").replace(r"\\", "").replace("&", "")
    raw = raw.replace("\\ ", " ")
    x, _, _, p, f, _, _ = polynomial_data()
    m, z, q, r = x
    h = q*q/m+m*m*q/12+m*z*z/4-m**5/720
    names = {str(v): v for v in x}
    names["H"] = h
    return tuple(s.sympify(re.sub(r"(?<=[mzqrH])(?=[mzqrH])", "*",
                                 term.strip().replace("^", "**")), locals=names)
                 for term in raw.split(","))


def interval_residuals(witness):
    x, substitutions, _ = interval_data()
    m = x[0]
    _, _, _, _, _, _, features = polynomial_data()
    theta, pencils = parse_pencils(SOURCE)
    assignment = dict(zip(theta, (s.sympify(v).xreplace(substitutions) for v in witness)))
    residuals = []
    for index, (matrix, (_, vector)) in enumerate(zip(pencils, features)):
        b = s.Matrix(vector).xreplace(substitutions)
        target = (s.zeros(matrix.rows) if index < 2 else
                  b*b.T/(m if index in (2, 3, 4) else 1))
        residuals.extend(s.cancel(v) for v in matrix.xreplace(assignment)-target)
    return residuals


class ReviewClarificationChecks(unittest.TestCase):
    def test_normalized_basis_coordinate_identity(self):
        (m, z, q, r), delta, n, p, f, basis, _ = polynomial_data()
        expected = (1, m, z, q, r, f, n*n/p, n/m, delta/m, p, delta,
                    q*q/m, m*q, q*z/m, m*m, m**3, m*z, m*z*z, z*z)
        self.assertEqual(len(expected), 19)
        self.assertTrue(all(s.cancel(b/p-e) == 0
                            for b, e in zip(basis, expected)))
        # Check the newly displayed transcription as well as the algebra.
        compact = re.sub(r'\s+', '', SOURCE)
        for fragment in (r'&(1,m,z,q,r,F,\N^2/p,N/m,\Delta/m,p,\Delta,q^2/m,',
                         r'&\hspace{14mm}mq,qz/m,m^2,m^3,mz,mz^2,z^2),'):
            self.assertIn(fragment, compact)

    def test_center_value_and_concavity_minor(self):
        (m, z, q, r), delta, *_ = polynomial_data()
        kappa = delta/m
        self.assertEqual(s.cancel(kappa.subs(
            {m:s.Rational(1,2), z:s.Rational(1,4), q:s.Rational(1,6)})-s.Rational(1,32)), 0)
        hessian = s.hessian(kappa, (m,z))
        self.assertEqual(s.cancel(hessian.det()-1), 0)
        self.assertEqual(s.cancel(hessian[0,0]+2*z*z/m**3+m/2), 0)

    def test_interval_boundary_formulas(self):
        x, substitutions, fifth = interval_data()
        _, delta, n, p, f, _, _ = polynomial_data()
        h = s.cancel(f-n*n/p)
        for residual in (delta.xreplace(substitutions), n.xreplace(substitutions),
                         h.xreplace(substitutions)-fifth):
            self.assertEqual(s.cancel(residual), 0)

    def test_all_67_entries_of_witness_parsed_from_manuscript(self):
        witness = source_witness(SOURCE)
        self.assertEqual(len(witness), 19)
        residuals = interval_residuals(witness)
        self.assertEqual(len(residuals), 67)
        self.assertTrue(all(v == 0 for v in residuals))

    def test_corrupted_witness_is_detected(self):
        witness = list(source_witness(SOURCE))
        witness[11] += 1
        self.assertTrue(any(v != 0 for v in interval_residuals(witness)))

    def test_all_zero_path_limits_and_apex_blocks(self):
        x, delta, n, p, f, basis, _ = polynomial_data()
        eps = s.Symbol("eps", positive=True)
        center = [s.Rational(1,2),s.Rational(1,4),s.Rational(1,6),s.Rational(1,8)]
        path = dict(zip(x, (eps*c for c in center)))
        for residual in (delta.xreplace(path)-eps**2*(4-eps**2)/192,
                         n.xreplace(path)-delta.xreplace(path),
                         p.xreplace(path)-eps**3*(4-eps**2)/384):
            self.assertEqual(s.cancel(residual), 0)
        values = [s.limit(s.cancel(v.xreplace(path)/p.xreplace(path)), eps, 0, dir="+")
                  for v in basis]
        self.assertEqual(values, [1]+[0]*18)
        theta, pencils = parse_pencils(SOURCE)
        for i, matrix in enumerate(pencils):
            desired = s.zeros(matrix.rows) if i < 5 else s.diag(0,0,1)
            self.assertEqual(matrix.xreplace(dict(zip(theta,values))), desired)

    def test_generic_restrictions_remain_exact_at_symbolic_peak(self):
        values = generic.moment_data()
        blocks = generic.generic_matrices(values)
        self.assertEqual([a.rows for a in blocks], [6,3,3,2,2])
        variables = (generic.s,generic.t,generic.q,generic.r)+tuple(
            values[e] for e in generic.EXPONENTS if e not in ((0,0),(1,0),(0,1)))
        self.assertEqual(len(variables), 16)
        peak = s.Symbol("J")
        def homogenize(mat):
            return mat.applyfunc(lambda e: s.expand(
                e+(peak-1)*e.subs(dict.fromkeys(variables,s.S.Zero))))
        hblocks = list(map(homogenize, blocks))
        changes = [s.Matrix([[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,0,1,0,0]])]
        changes += [s.Matrix([[1,0,0,0,0,0],[0,1,0,0,0,0],
                             [0,0,1,s.Rational(sign,6),0,0]]) for sign in (1,-1)]
        select = s.Matrix([[1,0,0],[0,1,0]])
        derived = [c*hblocks[0]*c.T for c in changes]
        derived += [select*b*select.T for b in hblocks[1:3]]+hblocks[3:]
        aux = [values[e] for e in ((2,0),(3,0),(4,0),(1,1),(2,1),(0,2))]
        target = [homogenize(s.Matrix(b)) for b in
                  small_matrices(variables[:4],aux)]
        for a,b in zip(derived,target):
            self.assertEqual((a-b).applyfunc(s.expand),s.zeros(a.rows))

    def test_generic_definition_and_certification_scope_retained(self):
        for phrase in (r"\Lambda(s^it^j)=\mu_{ij}", r"\mu_{00}=1",
                       "already linearized pencils", "specified exactly as rationals",
                       "not the original real-data equalities",
                       "Rational input data do not guarantee"):
            self.assertIn(phrase, SOURCE)

    def test_received_script_and_replayed_results_match(self):
        folder = ROOT/"research"/"math_derivation_audit"/"20261002"
        files = [folder/part/"verify_derivations.py" for part in ("received","rerun")]
        hashes = [hashlib.sha256(p.read_bytes()).hexdigest() for p in files]
        self.assertEqual(hashes, ["0c3cbb19793058438e058d36580dce7b093387f34fd968515b32db6ef17091d5"]*2)
        reports = [json.loads((folder/part/"symbolic_audit_results.json").read_text(
            encoding="utf-8")) for part in ("received","rerun")]
        self.assertEqual(reports[0],reports[1])
        self.assertEqual(reports[1]["number_of_named_checks"],49)
        self.assertTrue(reports[1]["all_checks_passed"])


if __name__ == "__main__":
    unittest.main()
