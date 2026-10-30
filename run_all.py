"""Reproduce the analyses and check the deposited results, figures and supplement."""
from pathlib import Path
import argparse
import hashlib
import shutil
import subprocess
import sys

ROOT = Path(__file__).absolute().parent
SCRIPTS = ROOT / "scripts/analysis"


def find_rscript(given):
    if given:
        return given
    found = shutil.which("Rscript")
    if found:
        return found
    candidates = sorted(Path(r"C:\Program Files\R").glob("R-*\\bin\\Rscript.exe"), reverse=True)
    if candidates:
        return str(candidates[0])
    raise SystemExit("Rscript not found; pass --rscript with its full path")


def publication_files():
    roots = [ROOT / name for name in ["data/inputs", "data/outputs", "figures", "scripts", "supplementary"]]
    files = [path for root in roots for path in root.rglob("*")
             if path.is_file() and "__pycache__" not in path.parts and path.suffix != ".pyc"]
    files += [ROOT / name for name in [".gitattributes", ".gitignore", "DATA_LICENSE.md", "LICENSE",
                                      "README.md", "requirements.txt", "run_all.py"]]
    return sorted(files, key=lambda path: path.relative_to(ROOT).as_posix())


def refresh_checksums():
    """Write publication-file checksums."""
    manifest = "".join(f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.relative_to(ROOT).as_posix()}\n"
                       for path in publication_files())
    (ROOT / "checksums.sha256").write_text(manifest, encoding="utf-8", newline="\n")


def verify_checksums():
    entries = {}
    for line in (ROOT / "checksums.sha256").read_text(encoding="utf-8").splitlines():
        digest, name = line.split("  ", 1)
        if name in entries:
            raise AssertionError(f"Duplicate checksum entry: {name}")
        path = ROOT / name
        if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise AssertionError(f"Checksum mismatch: {name}")
        entries[name] = digest
    expected = {path.relative_to(ROOT).as_posix() for path in publication_files()}
    if set(entries) != expected:
        raise AssertionError("Checksum manifest does not match the publication file inventory")
    print(f"Verified {len(entries)} release checksums.", flush=True)


def verify_figures(generated):
    from PIL import Image, ImageChops
    expected = {path.name for path in (ROOT / "figures").iterdir() if path.suffix.lower() in {".png", ".tif"}}
    actual = {path.name for path in generated.iterdir() if path.suffix.lower() in {".png", ".tif"}}
    if actual != expected:
        raise AssertionError("Regenerated figure inventory differs from the release")
    for name in sorted(expected):
        with Image.open(ROOT / "figures" / name) as shipped, Image.open(generated / name) as rebuilt:
            if shipped.size != rebuilt.size or ImageChops.difference(shipped.convert("RGB"), rebuilt.convert("RGB")).getbbox():
                raise AssertionError(f"Regenerated figure differs: {name}")
    expected_pdf = {path.stem for path in (ROOT / "figures").glob("*.png")}
    actual_pdf = {path.stem for path in generated.glob("*.pdf")}
    if actual_pdf != expected_pdf:
        raise AssertionError("Each figure must have a regenerated PDF")
    for stem in sorted(expected_pdf):
        if (ROOT / "figures" / (stem + ".pdf")).read_bytes() != (generated / (stem + ".pdf")).read_bytes():
            raise AssertionError(f"Regenerated PDF differs: {stem}")
    print("All regenerated figure pixels and PDFs match the release.", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--python", default=sys.executable)
    parser.add_argument("--rscript")
    parser.add_argument("--update-checksums", action="store_true", help="Update checksums and exit")
    args = parser.parse_args()
    if args.update_checksums:
        refresh_checksums()
        verify_checksums()
        return
    verify_checksums()
    rscript = find_rscript(args.rscript)
    derived = ROOT / "data/derived"
    derived.mkdir(parents=True, exist_ok=True)
    steps = ["01_prepare_data.py", "02_treatment_models.R", "03_quartile_models.R",
             "04_comparator_models.R", "05_developmental_models.R", "06_external_models.R",
             "07_survival_context.py", "08_weight_response.py", "make_figures.py",
             "make_tables.py", "verify_reported_results.py"]
    import tempfile
    with tempfile.TemporaryDirectory(prefix="figure_check_", dir=derived) as figure_directory:
        generated = Path(figure_directory)
        for name in steps:
            command = [rscript] if name.endswith(".R") else [args.python, "-B"]
            command += [str(SCRIPTS / name)]
            if name == "make_figures.py":
                command += ["--package", str(ROOT), "--output", str(generated)]
            if name == "make_tables.py":
                command += ["--check"]
            print("\n$ " + " ".join(command), flush=True)
            subprocess.run(command, cwd=ROOT, check=True)
        verify_figures(generated)
    verify_checksums()
    print("All reproduction checks passed.")


if __name__ == "__main__":
    main()
