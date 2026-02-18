from django.shortcuts import render
from django.http import JsonResponse, HttpResponse, StreamingHttpResponse
from django.conf import settings
from django.views.decorators.http import require_GET

import os
import datetime
import random
import string
import json
import sys
import subprocess
import threading
import time
import hashlib
from typing import Optional
from urllib import request as urlrequest, parse as urlparse, error as urlerror
import tempfile
import zipfile
import shutil
import io
import json
import xml.etree.ElementTree as ET

# Local helpers for taxonomy validation
from sampling_backend import taxon_lookup
from sampling_backend.classes_ranks_definition import RankType, RANK_ORDER

# --------------------------
# Helper functions
# --------------------------

def generate_job_id():
    timestamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S-")
    random_string = ''.join(random.choices(string.ascii_letters + string.digits, k=10))
    return timestamp + random_string

def create_email_directory(email="guest"):
    safe_email = email.replace("@", "_at_").replace(".", "_dot_")
    email_dir = os.path.join(settings.MEDIA_ROOT, safe_email)
    os.makedirs(email_dir, exist_ok=True)
    return email_dir

def parse_lines(text: str):
    if not text:
        return []
    return [ln.strip() for ln in text.splitlines() if ln.strip()]

def _progress_path(job_dir: str) -> str:
    return os.path.join(job_dir, "progress.json")

def _result_path(job_dir: str) -> str:
    return os.path.join(job_dir, "result.json")

def write_progress(job_dir: str, percent: int, message: str):
    try:
        with open(_progress_path(job_dir), "w", encoding="utf-8") as f:
            json.dump({"percent": int(percent), "message": str(message)}, f)
    except Exception:
        pass  # must never break the run

def write_result(job_dir: str, payload: dict):
    try:
        with open(_result_path(job_dir), "w", encoding="utf-8") as f:
            json.dump(payload, f)
    except Exception:
        pass

def read_result(job_dir: str):
    p = _result_path(job_dir)
    if not os.path.exists(p):
        return None
    try:
        with open(p, "r", encoding="utf-8") as f:
            return json.load(f)
    except Exception:
        return None

def _file_download_url(abs_path: str):
    if not abs_path or not os.path.exists(abs_path):
        return None
    rel = os.path.relpath(abs_path, settings.MEDIA_ROOT).replace("\\", "/")
    return settings.MEDIA_URL + rel

def _cache_root():
    p = os.path.join(settings.MEDIA_ROOT, "_cache")
    os.makedirs(p, exist_ok=True)
    return p

def _cache_key(payload: dict) -> str:
    """
    Stable hash of all inputs that affect sampling.
    """
    s = json.dumps(payload, sort_keys=True, ensure_ascii=True)
    return hashlib.sha256(s.encode("utf-8")).hexdigest()

def _cache_dir_for(key: str) -> str:
    return os.path.join(_cache_root(), key)


# --------------------------
# Output shaping
# --------------------------

def _extract_selected_genomes(text: str) -> str:
    """
    Strategy
    --------
    - Find the first line that starts with "Selected " and contains
      "genomes". Rephrase it to "Selected <N> genomes" (drop timing and
      node counts).
    - After that header, keep only bullet lines that describe genomes
      (lines beginning with a hyphen after optional whitespace).
    - Stop once we leave the bullet list so later sections like "Total Time"
      are omitted.
    """
    if not text:
        return ""

    lines = text.splitlines()
    header = None
    genome_lines: list[str] = []
    capturing = False

    for ln in lines:
        stripped = ln.strip()

        if not capturing:
            if stripped.startswith("Selected ") and "genomes" in stripped:
                # Try to normalize "Selected X genomes" header
                count = None
                try:
                    count = int(stripped.split()[1])
                except Exception:
                    count = None

                if count is not None:
                    plural = "genomes" if count != 1 else "genome"
                    header = f"Selected {count} {plural}"
                else:
                    header = stripped

                capturing = True
            continue

        # Once capturing, keep only genome bullet lines.
        if stripped.startswith("-"):
            genome_lines.append(stripped)
            continue

        # Stop if we reached the end of the bullet block.
        if genome_lines and stripped:
            break

    if not header:
        return ""

    if genome_lines:
        return "\n".join([header, *genome_lines])
    return header


def _extract_accessions(text: str) -> list[str]:
    """
    Best-effort parse of assembly accessions (GCF_/GCA_) from the sampling report.
    """
    if not text:
        return []
    import re
    pattern = re.compile(r"\b(GC[AF]_[0-9]+\.[0-9]+)\b")
    seen = set()
    accs = []
    for match in pattern.finditer(text):
        acc = match.group(1)
        if acc not in seen:
            seen.add(acc)
            accs.append(acc)
    return accs


def _slim_result_payload(payload: Optional[dict]) -> Optional[dict]:
    """Ensure the response only exposes the selected-genome section."""
    if not payload or payload.get("status") != "done":
        return payload

    selected_only = _extract_selected_genomes(
        payload.get("report_text") or payload.get("stdout") or ""
    )

    if not selected_only:
        return payload

    slim = payload.copy()
    slim["report_text"] = selected_only
    slim["stdout"] = selected_only
    return slim


# --------------------------
# Background job runner
# --------------------------

def _run_sampling_job(job_id: str, job_dir: str, cmd: list[str], report_path: str):
    """
    Runs Helene's sampling tool in a background thread.
    Canonical output: plain-text report written via --out.

    Key behaviors:
    - Stream stdout/stderr live (Popen) so we can see where it hangs.
    - Watchdog/heartbeat: if no output for a while, update progress message.
    - Mirror stdout into report_path if report file is empty (some tool versions don't write --out reliably).
    """
    write_progress(job_dir, 20, "Running sampling tool...")

    proc = None
    try:
        # Ensure unbuffered behavior in child
        env = os.environ.copy()
        env["PYTHONUNBUFFERED"] = "1"

        proc = subprocess.Popen(
            cmd,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1,
            env=env,
        )

        # Send a small set of inputs to unblock common interactive prompts.
        # - "0" for "choose index"
        # - "y" for "proceed?"
        try:
            if proc.stdin:
                proc.stdin.write("0\n")
                proc.stdin.write("y\n")
                proc.stdin.flush()
        except Exception:
            pass

        stdout_lines: list[str] = []
        stderr_lines: list[str] = []

        # Shared state for watchdog
        start_time = time.time()
        state_lock = threading.RLock()
        last_output_time = start_time
        last_output_msg = "Sampling started, waiting for output..."
        current_percent = 25  # we will slowly advance this during the run

        def _write_progress(pct: int, msg: str):
            write_progress(job_dir, pct, msg)

        def _bump_progress(amount: int) -> int:
            """
            Increment progress while capped at 85% (finalize step handles 90/100).
            """
            nonlocal current_percent
            with state_lock:
                current_percent = min(85, current_percent + amount)
                return current_percent

        _write_progress(current_percent, last_output_msg)

        def _drain(pipe, sink, prefix: str):
            nonlocal last_output_time, last_output_msg
            if not pipe:
                return
            for line in pipe:
                sink.append(line)
                msg = line.strip()
                if msg:
                    with state_lock:
                        last_output_time = time.time()
                        last_output_msg = (prefix + msg)[:200]
                        pct = _bump_progress(2)
                        msg_to_write = last_output_msg
                    _write_progress(pct, msg_to_write)

        def _watchdog():
            """
            If the subprocess is running but not printing, keep UI alive with a heartbeat.
            """
            nonlocal last_output_time, last_output_msg
            while proc and proc.poll() is None:
                time.sleep(5)
                now = time.time()
                with state_lock:
                    silence = now - last_output_time
                    elapsed = now - start_time
                if silence >= 30:
                    with state_lock:
                        pct = _bump_progress(1)
                        msg = f"Still running... (no output for {int(silence)}s, elapsed {int(elapsed)}s)"
                    _write_progress(pct, msg)

        t_out = threading.Thread(target=_drain, args=(proc.stdout, stdout_lines, ""), daemon=True)
        t_err = threading.Thread(target=_drain, args=(proc.stderr, stderr_lines, "ERR: "), daemon=True)
        t_wd = threading.Thread(target=_watchdog, daemon=True)

        t_out.start()
        t_err.start()
        t_wd.start()

        returncode = proc.wait()

        # Let drain threads finish quickly (pipes should close after wait)
        t_out.join(timeout=2)
        t_err.join(timeout=2)

        stdout_text = "".join(stdout_lines)
        stderr_text = "".join(stderr_lines)

        # Ensure the report file exists: if empty/missing, mirror stdout into it.
        try:
            if report_path:
                needs_fill = (not os.path.exists(report_path)) or (os.path.getsize(report_path) == 0)
                if needs_fill and stdout_text:
                    with open(report_path, "w", encoding="utf-8", errors="replace") as f:
                        f.write(stdout_text)
        except Exception:
            pass

        report_url = _file_download_url(report_path)

        report_text = None
        if report_path and os.path.exists(report_path):
            try:
                with open(report_path, "r", encoding="utf-8", errors="replace") as f:
                    report_text = f.read()
            except Exception:
                report_text = None

        if returncode != 0:
            write_progress(job_dir, 100, "Failed.")
            write_result(job_dir, {
                "status": "failed",
                "error": "Sampling tool failed.",
                "stderr": stderr_text,
                "stdout": stdout_text,
                "cmd": cmd,
                "report_download_url": report_url,
                "report_text": report_text,
                "job_id": job_id,
            })
            return

        # --- FINAL RESULT PAYLOAD (slimmed display; full download) ---
        write_progress(job_dir, 90, "Finalizing report...")

        # Prefer full stdout transcript as canonical report; fallback to file content
        final_report_text = stdout_text if stdout_text else (report_text or "")
        slim_report_text = _extract_selected_genomes(final_report_text) or final_report_text
        selected_accessions = _extract_accessions(final_report_text)

        # Ensure the downloadable report contains the full stdout (not just the slimmed view)
        try:
            if report_path and final_report_text:
                with open(report_path, "w", encoding="utf-8", errors="replace") as f:
                    f.write(final_report_text)
        except Exception:
            pass

        result_payload = {
            "status": "done",
            "message": "Sampling finished successfully.",
            "stderr": stderr_text,
            "stdout": slim_report_text,
            "full_stdout": final_report_text,
            "cmd": cmd,
            "report_download_url": report_url,
            "report_text": slim_report_text,
            "job_id": job_id,
            "selected_accessions": selected_accessions,
        }

        write_result(job_dir, result_payload)
        write_progress(job_dir, 100, "Done.")
        return

    except Exception as exc:
        # If something fails in the wrapper itself, attempt to terminate subprocess
        try:
            if proc and proc.poll() is None:
                proc.kill()
        except Exception:
            pass

        write_progress(job_dir, 100, "Failed.")
        write_result(job_dir, {"status": "failed", "error": str(exc), "job_id": job_id})
        return


# --------------------------
# Views
# --------------------------

def home(request):
    ranked_ranks = [rt.value for rt in RANK_ORDER]
    # Preserve Enum declaration order for ranks that are not part of the canonical ladder
    unranked_ranks = [rt.value for rt in RankType if rt not in RANK_ORDER]

    return render(
        request,
        "sampling_app/home.html",
        {
            "ranked_ranks": ranked_ranks,
            "unranked_ranks": unranked_ranks,
        },
    )

def progress(request, job_id: str):
    email_dir = create_email_directory("guest")
    job_dir = os.path.join(email_dir, job_id)
    p = _progress_path(job_dir)

    if not os.path.exists(p):
        return JsonResponse({"percent": 0, "message": "No progress information yet."})

    try:
        with open(p, "r", encoding="utf-8") as f:
            data = json.load(f)
        return JsonResponse({
            "percent": int(data.get("percent", 0)),
            "message": str(data.get("message", "")),
        })
    except Exception:
        return JsonResponse({"percent": 0, "message": "Progress file unreadable."})

def result(request, job_id: str):
    email_dir = create_email_directory("guest")
    job_dir = os.path.join(email_dir, job_id)

    data = _slim_result_payload(read_result(job_dir))
    if data is None:
        return JsonResponse({"status": "running"}, status=202)

    return JsonResponse(data)


@require_GET
def download_genome_ncbi(request):
    """
    Stream a genome ZIP directly from NCBI Datasets using our API key.
    The key stays server-side; the client only sees a standard file download.
    """
    accession = (request.GET.get("accession") or "").strip()
    if not accession:
        return HttpResponse("Missing required query parameter: accession", status=400)

    include = (request.GET.get("include") or getattr(settings, "NCBI_DATASETS_INCLUDE", "genome,gff3")).strip()
    # sanitize include list: keep comma-separated, drop spaces
    include = ",".join([part for part in include.replace(" ", "").split(",") if part]) or "genome"
    api_key = getattr(settings, "NCBI_API_KEY", None)
    if not api_key:
        return HttpResponse("NCBI_API_KEY not configured on server", status=500)

    base_url = f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{accession}/download"

    def fetch(params):
        full_url = f"{base_url}?{urlparse.urlencode(params)}"
        return urlrequest.urlopen(full_url, timeout=120)

    # First try with requested include; on 400 retry with genome-only to handle assemblies lacking gff3 etc.
    params = {"include": include, "api_key": api_key}
    fallback_params = {"include": "genome", "api_key": api_key}

    try:
        upstream = fetch(params)
    except urlerror.HTTPError as exc:
        if exc.code == 400 and include != "genome":
            try:
                upstream = fetch(fallback_params)
                include = "genome"  # note fallback used
            except Exception:
                return HttpResponse(f"Error contacting NCBI after fallback: {exc}", status=502)
        else:
            return HttpResponse(f"Error contacting NCBI: {exc}", status=502)
    except urlerror.URLError as exc:
        return HttpResponse(f"Error contacting NCBI: {exc}", status=502)

    status = getattr(upstream, "status", 200)
    if status >= 400:
        try:
            detail = upstream.read(500).decode("utf-8", errors="replace")
        except Exception:
            detail = "(no body)"
        return HttpResponse(f"NCBI error {status}: {detail}", status=status)

    filename = f"{accession}.zip"

    def _stream():
        try:
            while True:
                chunk = upstream.read(8192)
                if not chunk:
                    break
                yield chunk
        finally:
            try:
                upstream.close()
            except Exception:
                pass

    resp = StreamingHttpResponse(
        _stream(),
        content_type=upstream.headers.get_content_type() if hasattr(upstream, "headers") else "application/zip",
    )
    resp["Content-Disposition"] = f'attachment; filename="{filename}"'
    if upstream.headers.get("Content-Length"):
        resp["Content-Length"] = upstream.headers["Content-Length"]
    # Expose include used so client can debug
    resp["X-NCBI-Include"] = include
    return resp


def _serve_file(path: str, filename: Optional[str] = None):
    """
    Serve a local file as attachment via StreamingHttpResponse.
    """
    if not os.path.exists(path):
        return HttpResponse("File not found", status=404)
    filename = filename or os.path.basename(path)
    def _stream():
        with open(path, "rb") as f:
            while True:
                chunk = f.read(8192)
                if not chunk:
                    break
                yield chunk
    resp = StreamingHttpResponse(_stream(), content_type="application/zip")
    resp["Content-Disposition"] = f'attachment; filename="{filename}"'
    resp["Content-Length"] = os.path.getsize(path)
    return resp


def _download_single_ncbi_zip(accession: str, include: str, api_key: str, dest_path: str):
    """
    Download a single accession from NCBI Datasets into dest_path (zip file).
    Uses include list, with fallback to genome-only on HTTP 400 or when no FASTA is found.
    Verifies that at least one FASTA (*.fna) is present; otherwise retries with genome-only and
    raises if still missing.
    """
    base_url = f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{accession}/download"
    include = ",".join([part for part in include.replace(" ", "").split(",") if part]) or "genome"

    def fetch(params):
        full_url = f"{base_url}?{urlparse.urlencode(params)}"
        return urlrequest.urlopen(full_url, timeout=300)

    # Force fully hydrated packages (otherwise API may return metadata-only)
    base_params = {"hydrated": "FULLY_HYDRATED", "api_key": api_key}
    params = {**base_params, "include": include}
    fallback_params = {**base_params, "include": "genome,seq-report"}

    def _http_error_detail(exc: urlerror.HTTPError) -> str:
        try:
            body = exc.read(500).decode("utf-8", errors="replace")
        except Exception:
            body = "(no body)"
        return f"{exc.code} {exc.reason}: {body}"

    try:
        upstream = fetch(params)
    except urlerror.HTTPError as exc:
        if exc.code == 400 and include != "genome":
            # retry genome-only
            try:
                upstream = fetch(fallback_params)
                include = "genome"
            except urlerror.HTTPError as exc2:
                raise RuntimeError(f"HTTPError on fallback for {accession}: {_http_error_detail(exc2)}") from exc2
        else:
            raise RuntimeError(f"HTTPError for {accession}: {_http_error_detail(exc)}") from exc
    except urlerror.URLError as exc:
        raise RuntimeError(f"URLError for {accession}: {exc}") from exc

    with open(dest_path, "wb") as f:
        while True:
            chunk = upstream.read(8192)
            if not chunk:
                break
            f.write(chunk)
    try:
        upstream.close()
    except Exception:
        pass

    valid_exts = (".fna", ".fna.gz", ".fa", ".fa.gz", ".fasta", ".fasta.gz")

    def _extract_assembly_ftp(zf: zipfile.ZipFile) -> Optional[str]:
        """
        Try to read assembly_data_report.jsonl to get an FTP/HTTPS path base for manual FASTA fetch.
        """
        try:
            with zf.open("ncbi_dataset/data/assembly_data_report.jsonl") as f:
                for line in f:
                    try:
                        obj = json.loads(line)
                    except Exception:
                        continue
                    ftp = obj.get("refseq_ftp") or obj.get("genbank_ftp") or ""
                    if ftp:
                        return ftp
        except KeyError:
            return None
        except Exception:
            return None
        return None

    def _manual_fetch_fna(ftp_base: str, zf_path: zipfile.ZipFile) -> bool:
        """
        Download genomic FASTA from FTP/HTTPS using paths from assembly report and add to zip.
        Returns True if added.
        """
        # ftp base like ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/...
        if not ftp_base:
            return False
        http_base = ftp_base.replace("ftp://", "https://")
        fn_candidates = [
            ftp_base.split("/")[-1] + "_genomic.fna.gz",
            ftp_base.split("/")[-1] + "_genomic.fna",
        ]
        for fname in fn_candidates:
            url = f"{http_base}/{fname}"
            try:
                with urlrequest.urlopen(url, timeout=600) as resp:
                    data = resp.read()
                if not data:
                    continue
                # Add into zip under ncbi_dataset/data/<accession>/
                arcname = f"ncbi_dataset/data/{accession}/{fname}"
                zf_path.writestr(arcname, data)
                return True
            except Exception:
                continue
        return False

    def _eutils_fetch_ftp(accession: str) -> Optional[str]:
        """
        Use NCBI E-utilities to get the GenBank/RefSeq FTP path.
        """
        try:
            esearch_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=assembly&term={urlparse.quote(accession)}[Assembly+Accession]"
            with urlrequest.urlopen(esearch_url, timeout=30) as resp:
                root = ET.fromstring(resp.read())
            ids = [elem.text for elem in root.findall(".//IdList/Id")]
            if not ids:
                return None
            uid = ids[0]
            esum_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=assembly&id={uid}&retmode=json"
            with urlrequest.urlopen(esum_url, timeout=30) as resp:
                summary = json.load(resp)
            doc = summary.get("result", {}).get(uid, {})
            ftp = doc.get("ftppath_genbank") or doc.get("ftppath_refseq")
            return ftp
        except Exception:
            return None

    def validate_or_retry_for_fna(current_include: str):
        try:
            with zipfile.ZipFile(dest_path, "r") as zf:
                names = zf.namelist()
                has_fna = any(name.lower().endswith(valid_exts) for name in names)
        except zipfile.BadZipFile as exc:
            raise RuntimeError(f"Downloaded file is not a valid zip: {exc}") from exc

        if has_fna:
            return current_include

        # If no FASTA and we weren't already on genome-only, retry with genome-only
        if current_include.split(",")[0] != "genome":
            upstream2 = fetch(fallback_params)
            with open(dest_path, "wb") as f2:
                while True:
                    chunk = upstream2.read(8192)
                    if not chunk:
                        break
                    f2.write(chunk)
            try:
                upstream2.close()
            except Exception:
                pass
            # re-validate; if still none, raise
            with zipfile.ZipFile(dest_path, "a") as zf2:
                names2 = zf2.namelist()
                if any(name.lower().endswith(valid_exts) for name in names2):
                    return "genome"
                # Try manual fetch via FTP links in assembly report
                ftp_base = _extract_assembly_ftp(zf2)
                if ftp_base and _manual_fetch_fna(ftp_base, zf2):
                    return "genome"
                # Try eutils to get FTP path and fetch
                ftp_base2 = _eutils_fetch_ftp(accession)
                if ftp_base2 and _manual_fetch_fna(ftp_base2, zf2):
                    return "genome"

        # As a last resort, attempt manual fetch on the already-downloaded (non-fallback) zip
        try:
            with zipfile.ZipFile(dest_path, "a") as zf3:
                ftp_base = _extract_assembly_ftp(zf3)
                if ftp_base and _manual_fetch_fna(ftp_base, zf3):
                    return current_include
                ftp_base2 = _eutils_fetch_ftp(accession)
                if ftp_base2 and _manual_fetch_fna(ftp_base2, zf3):
                    return current_include
        except Exception:
            pass

        raise RuntimeError("No FASTA (.fna/.fa/.fasta) files found in NCBI package after fallback and FTP fetch.")

    return validate_or_retry_for_fna(include)


@require_GET
def download_selected_genomes(request, job_id: str):
    """
    Bundle all selected genomes for a completed job into a single ZIP.
    Each accession is fetched from NCBI Datasets using the configured API key.
    """
    include = (request.GET.get("include") or getattr(settings, "NCBI_DATASETS_INCLUDE", "genome,gff3")).strip()
    include = ",".join([part for part in include.replace(" ", "").split(",") if part]) or "genome"
    api_key = getattr(settings, "NCBI_API_KEY", None)
    if not api_key:
        return HttpResponse("NCBI_API_KEY not configured on server", status=500)

    # Locate job result to get selected accessions
    email_dir = create_email_directory("guest")
    job_dir = os.path.join(email_dir, job_id)
    data = read_result(job_dir)
    if not data or data.get("status") != "done":
        return HttpResponse("Job not finished or not found", status=404)

    accs = data.get("selected_accessions") or []
    if not accs:
        return HttpResponse("No accessions found in job output.", status=400)

    # Cache: if already built, serve it
    bundle_path = os.path.join(job_dir, "selected_genomes.zip")
    if os.path.exists(bundle_path) and os.path.getsize(bundle_path) > 0:
        return _serve_file(bundle_path, filename=f"{job_id}_selected_genomes.zip")

    tmp_dir = tempfile.mkdtemp(prefix="ncbi_bulk_", dir=job_dir)
    try:
        # Download each accession as {acc}.zip
        for acc in accs:
            dest = os.path.join(tmp_dir, f"{acc}.zip")
            try:
                _download_single_ncbi_zip(acc, include, api_key, dest)
            except Exception as exc:
                return HttpResponse(f"Failed to download {acc}: {exc}", status=502)

        # Bundle all into one zip
        with zipfile.ZipFile(bundle_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
            for fname in os.listdir(tmp_dir):
                full = os.path.join(tmp_dir, fname)
                zf.write(full, arcname=fname)

        return _serve_file(bundle_path, filename=f"{job_id}_selected_genomes.zip")
    finally:
        try:
            shutil.rmtree(tmp_dir)
        except Exception:
            pass

def run_sampling(request):
    """
    Starts a sampling run. Returns immediately with job_id.
    Canonical output is a report file written by Helene's tool via --out.
    """
    if request.method != "POST":
        return HttpResponse("Invalid request", status=405)

    # Required inputs (website)
    taxon = request.POST.get("taxon", "").strip()
    rank = request.POST.get("rank", "").strip()
    genomes = request.POST.get("genomes", "").strip()

    if not taxon:
        return JsonResponse({"error": "Taxon is required"}, status=400)
    if not rank:
        return JsonResponse({"error": "Rank is required"}, status=400)

    # Pre-validate taxon against local taxonomy JSONL to provide fast feedback
    taxonomy_path = getattr(settings, "TAXONOMY_JSON_PATH", None)
    resolved_records = []
    validation_error = None
    try:
        if taxonomy_path and os.path.exists(taxonomy_path):
            resolved_records = taxon_lookup.resolve_to_taxids(taxon, taxonomy_path)
            if not resolved_records:
                validation_error = "Unknown taxon name/ID (not found in taxonomy)."
        else:
            validation_error = "Taxonomy reference file missing on server."
    except Exception as exc:
        validation_error = f"Taxonomy validation error: {exc}"

    if validation_error:
        return JsonResponse({"error": validation_error}, status=400)

    try:
        per_taxon = int(genomes)
        if per_taxon < 1:
            raise ValueError
    except Exception:
        return JsonResponse({"error": "Genomes per taxon must be an integer >= 1."}, status=400)

    # Advanced inputs (mapped to Helene's actual CLI flags)
    tree = request.POST.get("tree", "phantom").strip() or "phantom"
    method = request.POST.get("method", "DFS").strip() or "DFS"
    prefer_reference = request.POST.get("prefer_reference") is not None
    prefer_higher_level = request.POST.get("prefer_higher_level") is not None
    min_assembly_level = request.POST.get("min_assembly_level", "").strip()

    seed_raw = request.POST.get("seed", "").strip()
    seed = None
    if seed_raw:
        try:
            seed = int(seed_raw)
        except Exception:
            return JsonResponse({"error": "Seed must be an integer."}, status=400)

    exclude_names = parse_lines(request.POST.get("exclude_names", ""))
    exclude_taxids_raw = parse_lines(request.POST.get("exclude_taxids", ""))

    exclude_taxids = []
    for ln in exclude_taxids_raw:
        try:
            exclude_taxids.append(int(ln))
        except Exception:
            return JsonResponse({"error": f"Exclude TaxID must be numeric. Invalid line: '{ln}'"}, status=400)

    # Note: exclude_file is not in the UI yet.
    exclude_file_path = None

    # Validate enums (prevents silent wrong calls)
    if tree not in {"basic", "phantom"}:
        return JsonResponse({"error": "Tree must be 'basic' or 'phantom'."}, status=400)

    if method not in {"DFS", "list", "bisect", "sibling"}:
        return JsonResponse({"error": "Method must be one of: DFS, list, bisect, sibling."}, status=400)

    allowed_levels = {"", "COMPLETE GENOME", "CHROMOSOME", "SCAFFOLD", "CONTIG"}
    if min_assembly_level not in allowed_levels:
        return JsonResponse({"error": "min_assembly_level must be one of COMPLETE GENOME, CHROMOSOME, SCAFFOLD, CONTIG."}, status=400)

    # Create job directory
    email_dir = create_email_directory("guest")
    job_id = generate_job_id()
    job_dir = os.path.join(email_dir, job_id)
    os.makedirs(job_dir, exist_ok=True)

    # --------------------------
    # Cache lookup (MVP tree caching)
    # --------------------------
    cache_payload = {
        "taxon": taxon,
        "rank": rank,
        "per_taxon": per_taxon,
        "tree": tree,
        "method": method,
        "prefer_reference": prefer_reference,
        "prefer_higher_level": prefer_higher_level,
        "min_assembly_level": min_assembly_level,
        "seed": seed,
        "exclude_names": exclude_names,
        "exclude_taxids": exclude_taxids,
    }

    cache_key = _cache_key(cache_payload)
    cache_dir = _cache_dir_for(cache_key)
    cache_result = os.path.join(cache_dir, "result.json")

    if os.path.exists(cache_result):
        # Cache hit: return immediately
        try:
            with open(cache_result, "r", encoding="utf-8") as f:
                cached = json.load(f)

            cached = _slim_result_payload(cached)

            # Ensure selected_accessions exists even for older cache entries
            if "selected_accessions" not in cached:
                cached_report = cached.get("full_stdout") or cached.get("report_text") or ""
                cached["selected_accessions"] = _extract_accessions(cached_report)

            write_progress(job_dir, 100, "Done (cached result).")
            write_result(job_dir, cached)

            return JsonResponse({
                "status": "started",
                "message": "Cached result reused.",
                "job_id": job_id,
                "cached": True,
            })
        except Exception:
            pass  # fall through to fresh run if cache is unreadable

    report_path = os.path.join(job_dir, "sampling_results.txt")
    write_progress(job_dir, 5, "Job created. Preparing command...")

    try:
        sampling_py = str(settings.BASE_DIR_BACKEND / "sampling.py")

        cmd = [
            sys.executable,
            "-u",  # IMPORTANT: unbuffered output
            sampling_py,
            "--tree", tree,
            "--out", report_path,
            "sample",
            "--query", taxon,
            "--rank", rank,
            "--per_taxon", str(per_taxon),
            "--method", method,
        ]

        if prefer_reference:
            cmd.append("--prefer_reference")
        if prefer_higher_level:
            cmd.append("--prefer_higher_level")
        if min_assembly_level:
            cmd.extend(["--min_assembly_level", min_assembly_level])
        if seed is not None:
            cmd.extend(["--seed", str(seed)])

        if exclude_names:
            cmd.append("--exclude_name")
            cmd.extend(exclude_names)
        if exclude_taxids:
            cmd.append("--exclude_taxid")
            cmd.extend([str(x) for x in exclude_taxids])

        if exclude_file_path:
            cmd.extend(["--exclude_file", exclude_file_path])

        write_progress(job_dir, 10, "Job queued. Starting...")

        t = threading.Thread(
            target=_run_sampling_job,
            args=(job_id, job_dir, cmd, report_path),
            daemon=True
        )
        t.start()

    except Exception as exc:
        write_progress(job_dir, 100, "Failed.")
        write_result(job_dir, {"status": "failed", "error": str(exc), "job_id": job_id})
        return JsonResponse({"error": str(exc), "job_id": job_id}, status=500)

    return JsonResponse({
        "status": "started",
        "message": "Job started.",
        "job_id": job_id,
    })
