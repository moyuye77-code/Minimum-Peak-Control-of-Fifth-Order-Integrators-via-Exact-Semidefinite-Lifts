"""Regression tests for the format checker; no scientific correctness claim."""
from pathlib import Path
import json
import re
import unittest
from research.paper.check_tex_structure import check_source, braced_argument

SOURCE = Path(__file__).with_name('representability.tex').read_text(encoding='utf-8')


class AmoFormatChecks(unittest.TestCase):
    def test_versioned_repository_archive(self):
        record = json.loads(Path(__file__).with_name('ONLINE_RESOURCE_1.json').read_text(encoding='utf-8'))
        self.assertEqual(SOURCE.count(record['archive_url']), 2)
        self.assertIn(r'\url{' + record['repository_url'] + '}', SOURCE)
        self.assertEqual(SOURCE.count(r'\path{' + record['filename'] + '}'), 3)
        self.assertEqual(SOURCE.count('version ' + record['archive_version']), 2)
        self.assertIn(record['release_tag'], SOURCE)
        self.assertNotIn('Online Resource 1', SOURCE)
        self.assertIn(record['caption'], ' '.join(SOURCE.split()))
        self.assertIn(record['email'], SOURCE)
        self.assertIn('reproducibility guide maps the claims', SOURCE)

    def test_disclosure_status_is_reported_separately(self):
        result = check_source(SOURCE)
        self.assertEqual(result['ai_use_disclosure_present'],
                         r'\bmhead{Use of artificial intelligence}' in SOURCE)
        self.assertFalse(result['publication_compliance_verified'])

    def test_missing_disclosure_not_misreported_as_compliant(self):
        changed = SOURCE.replace(r'\bmhead{Use of artificial intelligence}',
                                 r'\bmhead{Other}')
        result = check_source(changed)
        self.assertFalse(result['ai_use_disclosure_present'])
        self.assertFalse(result['publication_compliance_verified'])

    def test_external_input_rejected(self):
        with self.assertRaisesRegex(AssertionError, 'external manuscript input'):
            check_source(SOURCE + r'\input{private-notes.tex}')

    def test_internal_path_rejected(self):
        with self.assertRaisesRegex(AssertionError, 'internal research path'):
            check_source(SOURCE + r'\path{research/private-log.md}')

    def test_internal_editorial_flags_absent(self):
        for phrase in ('Author-review version', 'Author versions inspected',
                       'publication declarations remain incomplete',
                       'Author contributions and approval', 'serialization run'):
            self.assertNotIn(phrase, SOURCE)

    def test_current_source(self):
        result = check_source(SOURCE)
        self.assertEqual(result['bibliography_keys'], 11)
        self.assertEqual(result['attributed_preprints_in_text'], 2)
        self.assertEqual(result['keywords'], 5)
        self.assertFalse(result['compilation_verified'])
        self.assertFalse(result['references_last'])
        self.assertTrue(result['declarations_after_references'])
        self.assertTrue(result['source_layout_valid'])
        self.assertEqual(result['check_scope'],
                         'Source structure and bibliography layout only')

    def test_balanced_macro(self):
        self.assertEqual(braced_argument(r'\abstract{A {nested} argument.}', r'\abstract'),
                         'A {nested} argument.')

    def test_duplicate_label_rejected(self):
        with self.assertRaisesRegex(AssertionError, 'duplicate labels'):
            check_source(SOURCE + r'\label{sec:control-problem}')

    def test_reordered_bibliography_rejected(self):
        pattern = r'(\\bibitem\{averkov\}.*?)(\\bibitem\{gr\}.*?)(?=\\bibitem)'
        changed, count = re.subn(pattern, lambda match: match[2] + match[1], SOURCE, flags=re.S)
        self.assertEqual(count, 1)
        with self.assertRaisesRegex(AssertionError, 'alphabetic'):
            check_source(changed)

    def test_lost_preprint_rejected(self):
        with self.assertRaisesRegex(AssertionError, 'preprint attribution'):
            check_source(SOURCE.replace('math/0607463v2', 'math/0000000'))

    def test_missing_declaration_rejected(self):
        with self.assertRaisesRegex(AssertionError, 'declaration'):
            check_source(SOURCE.replace(r'\bmhead{Funding}', r'\bmhead{Other}'))

    def test_missing_citation_rejected(self):
        with self.assertRaisesRegex(AssertionError, 'bibliography entry'):
            check_source(SOURCE + r'\cite{nonexistent}')

    def test_wrong_main_heading_rejected(self):
        with self.assertRaisesRegex(AssertionError, 'section order'):
            check_source(SOURCE.replace(r'\section{System modeling}',
                                        r'\section{Unrelated narrative}'))

    def test_post_reference_prose_rejected(self):
        heading = r'\section*{Statements and Declarations}'
        with self.assertRaisesRegex(AssertionError, 'unexpected prose'):
            check_source(SOURCE.replace(heading, 'Extra conclusion.\n' + heading))

    def test_declarations_before_references_rejected(self):
        heading = r'\section*{Statements and Declarations}'
        prefix, declarations = SOURCE.split(heading, 1)
        declarations = heading + declarations.split(r'\end{document}', 1)[0]
        changed = prefix.replace(r'\backmatter', declarations + '\n' + r'\backmatter')
        with self.assertRaisesRegex(AssertionError, 'declarations must follow references'):
            check_source(changed + r'\end{document}')

    def test_duplicate_declarations_heading_rejected(self):
        heading = r'\section*{Statements and Declarations}'
        with self.assertRaisesRegex(AssertionError, 'duplicate declarations'):
            check_source(SOURCE.replace(heading, heading + '\n' + heading))

    def test_proof_outside_theory_rejected(self):
        changed = SOURCE.replace(r'\backmatter',
                                 r'\begin{proof}Misplaced.\end{proof}' + '\n' +
                                 r'\backmatter')
        with self.assertRaisesRegex(AssertionError, 'proof outside'):
            check_source(changed)

    def test_validation_outside_experiments_rejected(self):
        label = r'\label{sec:support-validation}'
        changed = SOURCE.replace(label, '').replace(
            r'\section{System modeling}', label + '\n' + r'\section{System modeling}')
        with self.assertRaisesRegex(AssertionError, 'validation outside'):
            check_source(changed)

    def test_unintegrated_appendix_rejected(self):
        with self.assertRaisesRegex(AssertionError, 'appendices must be integrated'):
            check_source(SOURCE.replace(r'\backmatter', r'\appendix' + '\n' + r'\backmatter'))


if __name__ == '__main__':
    unittest.main()
