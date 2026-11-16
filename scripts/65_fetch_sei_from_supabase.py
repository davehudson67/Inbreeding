#!/usr/bin/env python3
"""Stage 65: securely fetch the exact private Supabase SEI archive.

Archive location was verified in Supabase Storage on 2026-11-16.
Provides three authentication methods, in order:
  1. SEI_SIGNED_URL (temporary signed Storage object URL)
  2. SEI_STORAGE_TOKEN (JWT/authorized API credential, kept on server)
  3. installed, authenticated Supabase CLI using --remote project-ref

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

def main():
    OUT.parent.mkdir(parents=True, exist_ok=True)
    if OUT.exists():
        print("Using existing cached archive:", OUT)
        verify(OUT)
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
            else:
                print("Downloading private Storage object via Supabase CLI")
                download_cli(temp)
            verify(temp)
            os.replace(temp, OUT)

    meta = {
        "project": PROJECT, "bucket": BUCKET, "object": OBJECT,
        "sha256": EXPECTED_SHA256, "bytes": EXPECTED_BYTES,
        "local": str(OUT), "verified": True,
        "original_trajectories": 3214000, "draws": 1000,
        "note": "0 onset means not progressed by last capture, NOT never in future",
    }
    manifest = OUT.parent / "sei_storage_source_verified.json"
    manifest.write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
    print("Local verified source:", OUT)
    print("Manifest:", manifest)

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print("SEI FETCH FAILED:", str(e), file=sys.stderr)
        sys.exit(1)
