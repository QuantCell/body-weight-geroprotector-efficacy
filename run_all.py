"""Run all analyses, rebuild the figures, check the reported values and compare files with checksums.sha256."""
import argparse
import csv
import hashlib
import io
import math
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).absolute().parent
SCRIPTS = ROOT / "scripts"
MANIFEST = ROOT / "checksums.sha256"

STEPS = [
    ("py", "01_prepare_data.py", []),
    ("R", "02_primary_models.R", []),
    ("R", "03_matched_quartiles.R", []),
    ("R", "04_comparator_models.R", []),
    ("R", "05_developmental_models.R", []),
    ("R", "06_olfactory_cohort.R", []),
    ("py", "07_survival_summaries.py", []),
    ("py", "08_population_context.py", []),
    ("py", "09_experiment_table.py", []),
    ("py", "10_make_figures.py", ["--package", ROOT, "--output", ROOT / "figures"]),
    ("py", "11_verify_reported_results.py", []),
]

WINDOWS_R = [Path(r"C:\Program Files\R") / f"R-{v}" / "bin" / "Rscript.exe"
             for v in ("4.5.1", "4.5.0", "4.4.2", "4.4.1")]

# rebuilt on every run, or local only
NOT_LISTED = {".git", "__pycache__", "derived", "run.log", MANIFEST.name}


def find_rscript(given):
    if given:
        return given
    found = shutil.which("Rscript")
    if found:
        return found
    for candidate in WINDOWS_R:
        if candidate.exists():
            return str(candidate)
    sys.exit("Rscript not found; pass --rscript with its path")


def current_checksums():
    files = (path for path in ROOT.rglob("*")
             if path.is_file() and not NOT_LISTED & set(path.relative_to(ROOT).parts))
    return {path.relative_to(ROOT).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest() for path in files}


def same_table(old, new):
    """Whether two CSV files agree, with numbers compared to a relative tolerance of 1e-9."""
    a, b = (list(csv.reader(io.StringIO(data.decode("utf-8")))) for data in (old, new))
    if len(a) != len(b) or any(len(x) != len(y) for x, y in zip(a, b)):
        return False
    for x, y in zip((v for row in a for v in row), (v for row in b for v in row)):
        if x != y:
            try:
                if not math.isclose(float(x), float(y), rel_tol=1e-9, abs_tol=1e-300):
                    return False
            except ValueError:
                return False
    return True


def compare_with_manifest(update, shipped=None):
    """Rewrite the manifest, or compare with it and return whether every difference is acceptable."""
    current = current_checksums()
    if update:
        lines = [f"{current[name]}  {name}\n" for name in sorted(current, key=str.lower)]
        MANIFEST.write_text("".join(lines), encoding="utf-8", newline="\n")
        print(f"checksums.sha256 rewritten for {len(current)} files")
        return True
    shipped = shipped or {}
    listed = dict(reversed(line.split("  ", 1)) for line in MANIFEST.read_text(encoding="utf-8").splitlines())
    changed = sorted(name for name in current.keys() & listed.keys() if current[name] != listed[name])
    added, missing = sorted(current.keys() - listed.keys()), sorted(listed.keys() - current.keys())
    failures = [f"not in checksums.sha256: {name}" for name in added]
    failures += [f"missing but listed in checksums.sha256: {name}" for name in missing]
    for name in changed:
        old = shipped.get(name)
        if name.startswith("figures/"):
            print(f"rendered differently from checksums.sha256: {name}")
        elif old is not None and hashlib.sha256(old).hexdigest() == listed[name] \
                and same_table(old, (ROOT / name).read_bytes()):
            print(f"within numerical tolerance of checksums.sha256: {name}")
        else:
            failures.append(f"differs from checksums.sha256: {name}")
    if not (changed or added or missing):
        print(f"All {len(current)} files match checksums.sha256")
    for line in failures:
        print(line)
    return not failures


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--python", default=sys.executable)
    parser.add_argument("--rscript", default=None)
    parser.add_argument("--update-checksums", action="store_true",
                        help="rewrite checksums.sha256 from the regenerated files")
    args = parser.parse_args()
    rscript = find_rscript(args.rscript)
    shipped = {path.relative_to(ROOT).as_posix(): path.read_bytes()
               for path in (ROOT / "data" / "outputs").glob("*.csv")}

    for number, (kind, name, extra) in enumerate(STEPS, start=1):
        executable = args.python if kind == "py" else rscript
        flags = ["-B"] if kind == "py" else []
        command = [executable, *flags, str(SCRIPTS / name)] + [str(a) for a in extra]
        print(f"[{number:>2}/{len(STEPS)}] {name}", flush=True)
        if subprocess.run(command).returncode != 0:
            sys.exit(f"failed at step {number}: {name}")
    if not compare_with_manifest(args.update_checksums, shipped):
        sys.exit("Reproduction check failed for the files listed above")
    print("Done.")


if __name__ == "__main__":
    main()
