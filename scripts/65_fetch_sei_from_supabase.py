#!/usr/bin/env python3
"""Stage 65: find a local extracted SEI bundle or fetch the original ZIP.

Archive location was verified in Supabase Storage on 2026-11-16.
Default is strictly local. Accepts the original verified ZIP or
extracted badgers.csv and trajectories_by_draw.csv in the same folder.
Missing local files trigger a directory listing, NOT a network attempt.
Remote download is explicitly opt-in through one of:
  1. SEI_SIGNED_URL (temporary signed Storage object URL)
  2. SEI_STORAGE_TOKEN (JWT/authorized API credential, kept on server)
  3. SEI_FETCH_REMOTE=1 with authenticated Supabase CLI.

No credentials, URLs containing tokens, or user data are printed.
Never writes to Supabase, the raw archive, or the PostgreSQL database.
"""
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from urllib.parse import unquote, urlsplit
from urllib.request import Request, urlopen
from zipfile import ZipFile

PROJECT = "klpvihubzlyzubyyubtw"
BUCKET = "research-model-outputs"
OBJECT = "SEI_only_research_bundle.zip"
EXPECTED_SHA256 = "65a1f1a7c84e2f9bfbbe149a3d92078cebd185e72281d7e5452cf0eadaaf9f98"
EXPECTED_BYTES = 23662164
OUT = Path("data/infection_posteriors/sei_2026/SEI_only_research_bundle.zip")
REQUIRED_MEMBERS = {
    "badgers.csv", "trajectories_by_draw.csv", "p_state_by_quarter.csv",
    "trajectories.csv", "trajectories_infectious.csv",
    "fits/sei_inferred_accuracy_.rds", "README_SEI_ONLY.md",
    "bundle_manifest.json",
}

def digest(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for part in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(part)
    return h.hexdigest()

def verify(path):
    size = path.stat().st_size
    checksum = digest(path)
    if size != EXPECTED_BYTES or checksum != EXPECTED_SHA256:
        raise RuntimeError(
            "SEI file differs from inspected source. Expected size "
            f"{EXPECTED_BYTES}, obtained {size}; expected SHA256 "
            f"{EXPECTED_SHA256}, obtained {checksum}. "
            "STOP: verify the uploaded Storage object before analysing."
        )
    with ZipFile(path) as z:
        actual = set(z.namelist())
        missing = REQUIRED_MEMBERS - actual
        if missing:
            raise RuntimeError(f"Archive missing members: {sorted(missing)}")
        bad_member = z.testzip()
        if bad_member is not None:
            raise RuntimeError(f"ZIP CRC failure: {bad_member}")
    print("Verified source ZIP size:", size)
    print("Verified SHA256:", checksum)
    print("ZIP CRC and required member checks: PASS")

def download_http(tmp, url, bearer):
    # Avoid sending tokens to arbitrary third-party hosts.
    p = urlsplit(url)
    if p.scheme != "https" or p.hostname != PROJECT + ".supabase.co":
        raise RuntimeError("Expected the exact HTTPS URL of the Woodchester Supabase project")
    if unquote(p.path).find("/" + BUCKET + "/" + OBJECT) < 0:
        raise RuntimeError("Unexpected Storage object path in URL")
    headers = {}
    if bearer:
        headers["Authorization"] = "Bearer " + bearer
        # Legacy secret keys can serve as apikey, while JWT user tokens
        # may need a separate publishable key for the gateway.
        api_key = os.environ.get("SEI_STORAGE_APIKEY", bearer)
        headers["apikey"] = api_key
    req = Request(url, headers=headers)
    # Never print exception URLs; they can embed time-limited tokens.
    try:
        with urlopen(req, timeout=180) as response, open(tmp, "wb") as dst:
            shutil.copyfileobj(response, dst, 1024 * 1024)
    except Exception as exc:
        raise RuntimeError(
            "HTTPS download failed (URL and authentication hidden). "
            "Check that the signed URL is still valid or Storage token "
            "has access. Underlying error type: " + type(exc).__name__
        ) from None

def download_cli(tmp):
    if not shutil.which("supabase"):
        raise RuntimeError(
            "No authenticated download method configured. Either install/login "
            "to Supabase CLI or provide SEI_SIGNED_URL or SEI_STORAGE_TOKEN "
            "as a shell environment variable. Do not paste credentials into chat."
        )
    cmd = [
        "supabase", "storage", "cp",
        f"ss:///{BUCKET}/{OBJECT}",
        str(tmp), "--experimental", f"--remote={PROJECT}",
    ]
    try:
        result = subprocess.run(cmd, text=True, capture_output=True, timeout=240)
    except Exception as exc:
        raise RuntimeError("Supabase CLI failed: " + type(exc).__name__) from None
    if result.returncode != 0:
        raise RuntimeError(
            "Supabase CLI could not download the private object. Run "
            "'supabase login' on this server, then retry. "
            "CLI stderr (limited): " + result.stderr[:1300]
        )
    if not tmp.exists():
        raise RuntimeError("Supabase CLI exited without creating the archive")

def find_extracted():
    """Locate a single self-contained extracted SEI folder under sei_2026."""
    candidates = sorted(
        p.parent for p in OUT.parent.rglob("badgers.csv")
        if p.is_file() and (p.parent / "trajectories_by_draw.csv").is_file()
    )
    # Avoid ambiguous version selection: never guess if two copies exist.
    if len(candidates) > 1:
        raise RuntimeError(
            "Multiple extracted SEI folders found: "
            + ", ".join(str(p) for p in candidates)
            + ". Keep one source or provide the original ZIP."
        )
    return candidates[0] if candidates else None


def main():
    OUT.parent.mkdir(parents=True, exist_ok=True)
    extracted = None
    if OUT.exists():
        print("Using existing cached archive:", OUT)
        verify(OUT)
    else:
        extracted = find_extracted()
        if extracted is not None:
            print("Found already-extracted SEI contents:", extracted)
            required = ("badgers.csv", "trajectories_by_draw.csv")
            for member in required:
                path = extracted / member
                if not path.is_file() or path.stat().st_size == 0:
                    raise RuntimeError("Missing/empty extracted source: " + str(path))
            print("Extracted essential tables: PRESENT")
            print("Original ZIP not available: cannot verify the Supabase ZIP SHA256.")
            print("Will fingerprint extracted files and validate rows, IDs, "
                  "draws, and chronology in Stage 66.")
        else:
            with tempfile.TemporaryDirectory(dir=str(OUT.parent)) as tempdir:
                temp = Path(tempdir) / OBJECT
                signed = os.environ.get("SEI_SIGNED_URL", "").strip()
                token = os.environ.get("SEI_STORAGE_TOKEN", "").strip()
                if signed:
                    print("Downloading private Storage object using time-limited URL")
                    download_http(temp, signed, None)
                elif token:
                    print("Downloading private Storage object using authorized bearer token")
                    download_http(
                        temp,
                        f"https://{PROJECT}.supabase.co/storage/v1/object/"
                        f"authenticated/{BUCKET}/{OBJECT}", token
                    )
                elif os.environ.get("SEI_FETCH_REMOTE", "").strip() == "1":
                    print("Explicit remote download requested: using Supabase CLI")
                    download_cli(temp)
                else:
                    root = OUT.parent
                    inventory = sorted(
                        p for p in root.rglob("*") if p.is_file()
                    ) if root.exists() else []
                    print("Local SEI directory:", root.resolve())
                    print("Current working directory:", Path.cwd())
                    print("Files found within local SEI folder:", len(inventory))
                    for candidate in inventory[:60]:
                        print("  ", candidate.relative_to(root),
                              "(", candidate.stat().st_size, "bytes)")
                    if len(inventory) > 60:
                        print("  ... more files omitted")
                    raise RuntimeError(
                        "LOCAL SEI FILES NOT RECOGNISED. Stage 65 needs "
                        "badgers.csv and trajectories_by_draw.csv in the "
                        "same directory anywhere under "
                        "data/infection_posteriors/sei_2026/. "
                        "Check filenames, the extraction folder, and "
                        "uppercase/lowercase spelling on Linux. "
                        "No Supabase access is needed for local files. "
                        "Remote download is now opt-in only."
                    )
                verify(temp)
                os.replace(temp, OUT)

    is_zip = OUT.is_file()
    meta = {
        "project": PROJECT, "bucket": BUCKET, "object": OBJECT,
        "source_mode": "original_zip" if is_zip else "extracted_files",
        "sha256": EXPECTED_SHA256 if is_zip else None,
        "bytes": EXPECTED_BYTES if is_zip else None,
        "local": str(OUT) if is_zip else str(extracted),
        "original_zip_checksum_verified": bool(is_zip),
        "verified": bool(is_zip),
        "original_trajectories": 3214000, "draws": 1000,
        "note": "0 onset means not progressed by last capture, NOT never in future",
    }
    if extracted is not None:
        paths = [extracted / "badgers.csv",
                 extracted / "trajectories_by_draw.csv"]
        meta["extracted_files_sha256"] = {
            p.name: digest(p) for p in paths
        }
        meta["extracted_files_bytes"] = {
            p.name: p.stat().st_size for p in paths
        }
        checksums = OUT.parent / "sei_extracted_files_sha256.tsv"
        with checksums.open("w",encoding="utf-8") as fh:
            fh.write("filename\tsha256\n")
            for p in paths:
                fh.write(f"{p.name}\t{digest(p)}\n")
        print("Wrote extracted-file SHA256 fingerprints:", checksums)
    manifest = OUT.parent / "sei_storage_source_verified.json"
    manifest.write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
    print("SEI source manifest:", manifest)
    print("Archive SHA256 verified:", is_zip)

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print("SEI FETCH FAILED:", str(e), file=sys.stderr)
        sys.exit(1)
