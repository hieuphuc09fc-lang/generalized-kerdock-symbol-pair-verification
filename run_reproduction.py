"""Run unchanged archival scripts in fresh, isolated output directories."""
import argparse
import ast
from datetime import datetime, timezone
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys
import uuid

ROOT = Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text(encoding="utf-8"))


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def preflight():
    """Check syntax and dependency imports without triggering script drivers."""
    files = list(ROOT.glob("*.py")) + list((ROOT / "additional_checks").glob("*.py"))
    for path in files:
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=path.name)
        compile(tree, str(path), "exec")
        # Several archival drivers run exhaustive loops at module scope.
        # Execute their import statements only, not their computation bodies.
        imports = [node for node in ast.walk(tree)
                   if isinstance(node, (ast.Import, ast.ImportFrom))]
        exec(compile(ast.Module(body=imports, type_ignores=[]), str(path), "exec"), {})
    return sorted(str(p.relative_to(ROOT)).replace("\\", "/") for p in files)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("job", choices=["smoke", "binary5", "extra", "consistency",
                                       "weights", "q8", "nonbinary", "all-primitive", "additional"])
    args = parser.parse_args()
    checked = preflight()
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    run = ROOT / "_runs" / (args.job + "_" + stamp + "_" + uuid.uuid4().hex[:8])
    work = run / "verification"
    work.mkdir(parents=True)
    for path in ROOT.glob("*.py"):
        if path.name != Path(__file__).name:
            shutil.copy2(path, work / path.name)
    for name in ("all_primitive.py", "verify_kerdock.py"):
        shutil.copy2(ROOT / "additional_checks" / name, work / name)
    if args.job == "consistency":
        # Consistency checks consume the immutable reference datasets.
        for name in ["results.json", "q8_results.json"]:
            shutil.copy2(ROOT / "results" / name, work / name)
    report = {"job": args.job, "syntax_and_dependency_import_checks": checked,
              "comparisons": []}
    if args.job in ("smoke", "binary5"):
        from verify import Field, audit_weights, binary_distance
        from field_audit import spectrum
        m = 3 if args.job == "smoke" else 5
        coeff = [1, 1, 0] if m == 3 else [1, 0, 1, 0, 0]
        f = Field(2, m, coeff)
        weights = audit_weights(f, 1)
        distance = binary_distance(f)
        independent = spectrum(2, m, coeff)
        saved = read(ROOT / "results/results.json")
        expected_independent = next(d for d in read(ROOT / "results/additional/verification-results.json")
                                    if (d["q"], d["m"]) == (2, m))
        normalize = lambda value: json.loads(json.dumps(value))
        assert normalize(weights) == saved[f"weights_2_{m}_1"]
        assert normalize(distance) == saved[f"distance_2_{m}"]
        assert normalize(independent) == expected_independent
        output = {"weights": weights, "distance": distance, "independent": independent}
        (run / "binary_check.json").write_text(json.dumps(output, indent=2), encoding="utf-8")
        report["comparisons"].append({"q": 2, "m": m, "minimum_weight": weights["min_weight"],
                                      "minimum_distance": distance["min_distance"],
                                      "pairs_per_implementation": distance["pairs"],
                                      "all_recorded_fields_match": True})
    else:
        scripts = {"extra": "extra_checks.py", "consistency": "consistency.py",
                   "weights": "verify.py", "nonbinary": "review2_checks.py",
                   "all-primitive": "all_primitive.py", "additional": "verify_kerdock.py"}
        outputs = {"extra": ["extra_results.json"], "consistency": ["consistency_results.json"],
                   "weights": ["results.json"], "q8": ["q8_results.json"],
                   "nonbinary": ["q4_m3_all_orders.json", "q4_m5_rechecked.json", "q8_m3_rechecked.json"],
                   "all-primitive": ["all_primitive_results.json"], "additional": ["verification-results.json"]}
        if args.job == "q8":
            code = ("import sys,json;from pathlib import Path;sys.path.insert(0,'verification');"
                    "from verify import Field,audit_weights;"
                    "r=audit_weights(Field(8,3,[2,7,2]),7);"
                    "Path('verification/q8_results.json').write_text(json.dumps(r,indent=2),encoding='utf-8')")
            command = [sys.executable, "-B", "-c", code]
        else:
            command = [sys.executable, "-B", str(work / scripts[args.job])]
        subprocess.run(command, cwd=run, check=True)
        for name in outputs[args.job]:
            reference = ROOT / "results" / ("additional" if args.job == "additional" else "") / name
            matches = read(work / name) == read(reference)
            report["comparisons"].append({"file": name, "matches_reference": matches})
            if not matches:
                (run / "check_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
                raise RuntimeError("Result differs from archived reference: " + name)
    (run / "check_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))
    print("Fresh outputs: " + str(run.relative_to(ROOT)))


if __name__ == "__main__":
    main()
