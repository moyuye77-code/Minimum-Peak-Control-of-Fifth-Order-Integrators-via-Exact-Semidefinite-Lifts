"""Adversarial checks for bibliography transcription, not scientific validity."""
from pathlib import Path
import copy
import json
import unittest

from research.paper.check_references import check

HERE = Path(__file__).resolve().parent
SOURCE = (HERE/"representability.tex").read_text(encoding="utf-8")
METADATA = {p.stem: json.loads(p.read_text(encoding="utf-8"))["message"]
            for p in (HERE/"reference_metadata").glob("*.json")}


class ReferenceChecks(unittest.TestCase):
    def test_current_fields(self):
        result = check(SOURCE, METADATA)
        self.assertEqual(result["references_checked"], 11)
        self.assertFalse(result["citation_relevance_automatically_verified"])

    def test_wrong_doi_rejected(self):
        with self.assertRaisesRegex(ValueError, "DOI mismatch"):
            check(SOURCE.replace("18M1201342", "18M1201343"), METADATA)

    def test_online_year_cannot_replace_volume_year(self):
        with self.assertRaisesRegex(ValueError, "volume year mismatch"):
            check(SOURCE.replace("109--118 (2019)", "109--118 (2018)"), METADATA)

    def test_wrong_pages_rejected(self):
        with self.assertRaisesRegex(ValueError, "pages mismatch"):
            check(SOURCE.replace("775--780", "775--781"), METADATA)

    def test_wrong_title_rejected(self):
        with self.assertRaisesRegex(ValueError, "title mismatch"):
            check(SOURCE.replace("Perfect spline solutions", "Perfect spline conjectures"),
                  METADATA)

    def test_incomplete_inventory_rejected(self):
        records = copy.deepcopy(METADATA)
        del records["fawzi"]
        with self.assertRaisesRegex(ValueError, "reference inventory"):
            check(SOURCE, records)

    def test_wrong_author_rejected(self):
        with self.assertRaisesRegex(ValueError, "author mismatch"):
            check(SOURCE.replace("Fawzi, H.", "Unrelated, H."), METADATA)


if __name__ == "__main__":
    unittest.main()
