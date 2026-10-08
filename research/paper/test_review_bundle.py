"""Small adversarial checks for the review inventory, not scientific proofs."""
import hashlib
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from research.paper.review_checks import verify_inventory


class BundleChecks(unittest.TestCase):
    def fixture(self, root, name="data.txt"):
        (root/"data.txt").write_bytes(b"original")
        record = {"schema": 1, "files": {name: {
            "bytes": 8, "sha256": hashlib.sha256(b"original").hexdigest()}}}
        (root/"bundle-manifest.json").write_text(json.dumps(record), encoding="utf-8")

    def test_valid_inventory(self):
        with TemporaryDirectory() as directory:
            root=Path(directory); self.fixture(root)
            self.assertEqual(verify_inventory(root), 1)

    def test_changed_bytes_rejected(self):
        with TemporaryDirectory() as directory:
            root=Path(directory); self.fixture(root)
            (root/"data.txt").write_bytes(b"mutated!")
            with self.assertRaisesRegex(ValueError, "Changed"):
                verify_inventory(root)

    def test_missing_file_rejected(self):
        with TemporaryDirectory() as directory:
            root=Path(directory); self.fixture(root, "missing.txt")
            with self.assertRaisesRegex(ValueError, "Missing"):
                verify_inventory(root)

    def test_escape_rejected(self):
        with TemporaryDirectory() as directory:
            root=Path(directory); self.fixture(root, "../data.txt")
            with self.assertRaisesRegex(ValueError, "Unsafe"):
                verify_inventory(root)

    def test_windows_path_rejected(self):
        with TemporaryDirectory() as directory:
            root=Path(directory); self.fixture(root, "C:/data.txt")
            with self.assertRaisesRegex(ValueError, "Unsafe"):
                verify_inventory(root)

    def test_unknown_schema_rejected(self):
        with TemporaryDirectory() as directory:
            root=Path(directory)
            (root/"bundle-manifest.json").write_text('{"schema": 2, "files": {}}')
            with self.assertRaisesRegex(ValueError, "Unsupported"):
                verify_inventory(root)

    def publication_fixture(self, root):
        self.fixture(root)
        path = root / "bundle-manifest.json"
        record = json.loads(path.read_text(encoding="utf-8"))
        identification = json.loads(Path(__file__).with_name("ONLINE_RESOURCE_1.json").read_text(
            encoding="utf-8"))
        record["publication"] = identification
        record["files"]["data.txt"]["publication"] = identification.copy()
        return path, record

    def test_publication_metadata(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            path, record = self.publication_fixture(root)
            path.write_text(json.dumps(record), encoding="utf-8")
            self.assertEqual(verify_inventory(root), 1)

    def test_missing_publication_email_rejected(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            path, record = self.publication_fixture(root)
            del record["publication"]["email"]
            path.write_text(json.dumps(record), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "Incomplete publication"):
                verify_inventory(root)

    def test_wrong_resource_number_rejected(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            path, record = self.publication_fixture(root)
            record["publication"]["filename"] = "ESM_2.zip"
            path.write_text(json.dumps(record), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "Wrong resource"):
                verify_inventory(root)

    def test_file_identification_mismatch_rejected(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            path, record = self.publication_fixture(root)
            record["files"]["data.txt"]["publication"]["email"] = "changed@example.invalid"
            path.write_text(json.dumps(record), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "File publication metadata mismatch"):
                verify_inventory(root)

    def test_inconsistent_release_url_rejected(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            path, record = self.publication_fixture(root)
            record["publication"]["archive_url"] += "-wrong"
            path.write_text(json.dumps(record), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "Release URL mismatch"):
                verify_inventory(root)

    def test_inconsistent_release_version_rejected(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            path, record = self.publication_fixture(root)
            record["publication"]["archive_version"] = "1900-01-01"
            path.write_text(json.dumps(record), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "Release version mismatch"):
                verify_inventory(root)

    def test_missing_file_identification_rejected(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            path, record = self.publication_fixture(root)
            del record["files"]["data.txt"]["publication"]
            path.write_text(json.dumps(record), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "File publication metadata mismatch"):
                verify_inventory(root)


if __name__ == "__main__":
    unittest.main()
