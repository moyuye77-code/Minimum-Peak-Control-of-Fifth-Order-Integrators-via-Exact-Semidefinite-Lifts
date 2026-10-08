"""Read-only bundle checks. Default uses Python's standard library only."""
import argparse
from contextlib import ExitStack, redirect_stdout
import hashlib
import io
import json
from pathlib import Path, PurePosixPath
import runpy
import sys
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[2]


def require(condition, message):
    if not condition:
        raise ValueError(message)


def verify_inventory(root):
    manifest = json.loads((root / "bundle-manifest.json").read_text(encoding="utf-8"))
    require(manifest.get("schema") == 1, "Unsupported inventory schema")
    publication = manifest.get("publication")
    if publication is not None:
        required = ("online_resource", "filename", "article_title", "journal", "authors",
                    "corresponding_author", "affiliation", "email", "caption",
                    "delivery", "repository_url", "archive_version", "release_tag", "archive_url")
        require(all(publication.get(key) for key in required), "Incomplete publication metadata")
        require(publication["online_resource"] == "Reproducibility archive"
                and publication["filename"] == "ESM_1.zip", "Wrong resource identification")
        require(publication["archive_url"] == publication["repository_url"]
                + "/releases/tag/" + publication["release_tag"], "Release URL mismatch")
        require(publication["release_tag"] == "reproducibility-" + publication["archive_version"],
                "Release version mismatch")
    for name, record in manifest["files"].items():
        if publication is not None:
            require(record.get("publication") == publication,
                    "File publication metadata mismatch: " + name)
        relative = PurePosixPath(name)
        require(not relative.is_absolute() and ".." not in relative.parts
                and "\\" not in name and ":" not in name, "Unsafe inventory path")
        path = root / relative
        require(path.resolve().is_relative_to(root.resolve()), "Escaped inventory path")
        require(path.is_file() and not path.is_symlink(), "Missing/linked file: " + name)
        data = path.read_bytes()
        require(len(data) == record["bytes"]
                and hashlib.sha256(data).hexdigest() == record["sha256"],
                "Changed bundle file: " + name)
    return len(manifest["files"])


def run_module(name, arguments=()):
    output = io.StringIO()
    with redirect_stdout(output), patch.object(sys, "argv", [name, *arguments]):
        runpy.run_module(name, run_name="__main__")
    return {"module": name, "passed": True, "output": output.getvalue()}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--numerical", action="store_true")
    args = parser.parse_args()
    require(__debug__, "Do not run checks with -O or PYTHONOPTIMIZE")
    sys.path.insert(0, str(ROOT))
    result = {"inventory_files": verify_inventory(ROOT), "checks": [],
              "new_optimizer_calls": 0, "new_timing_measurements": False}
    for name in ("research.paper.check_tex_structure", "research.paper.check_references",
                 "research.paper.reproduce"):
        result["checks"].append(run_module(name))
    if args.numerical:
        import cvxpy
        import scipy.optimize
        import numpy
        import scipy
        import sympy
        import clarabel
        result["versions"] = {m.__name__: m.__version__
                              for m in (numpy, scipy, sympy, cvxpy, clarabel)}
        with ExitStack() as stack:
            stack.enter_context(patch.object(cvxpy.Problem, "solve",
                side_effect=AssertionError("Review replay must not solve")))
            for name in ("linprog", "minimize", "least_squares", "differential_evolution"):
                stack.enter_context(patch.object(scipy.optimize, name,
                    side_effect=AssertionError("Review replay must not optimize")))
            for name in (
                "research.paper.check_lift_coefficients",
                "research.moment_four_lift.verify",
                "research.schur_remainder_recursion.verify",
                "research.moment_lift_obstruction.verify",
                "research.closed_five_lift.verify",
                "research.lift_value_benchmark.experiment",
                "research.peak_synthesis.experiment",
                "research.peak_recovery_ablation",
                "research.peak_feedback.experiment",
            ):
                result["checks"].append(run_module(name))
            result["checks"].append(run_module("research.theory_alignment.verify",
                ("--output", str(ROOT/"research/theory_alignment/results/run-20261002-r2"))))
            modules = [
                "research.paper.test_review_clarifications",
                "research.paper.test_reader_entrypoints", "research.paper.test_amo_format",
                "research.theory_alignment.test_alignment", "research.paper.test_references",
                "research.paper.test_review_bundle",
            ]
            stream = io.StringIO()
            suite = unittest.defaultTestLoader.loadTestsFromNames(modules)
            tested = unittest.TextTestRunner(stream=stream, verbosity=2).run(suite)
            require(tested.wasSuccessful(), stream.getvalue())
            result["tests"] = {"count": tested.testsRun, "passed": True, "log": stream.getvalue()}
    else:
        require(not ({"numpy", "scipy", "sympy", "cvxpy", "clarabel"} & set(sys.modules)),
                "Unexpected numerical dependency in standard-library checks")
    result["inventory_after_files"] = verify_inventory(ROOT)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
