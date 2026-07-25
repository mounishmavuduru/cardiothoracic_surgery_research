"""Anubis-aware downloader for Dryad file_stream endpoints.

Dryad fronts downloads with an Anubis proof-of-work wall (v1.24.0). The agent
proxy passes plain curl/requests through to Dryad but not a headless Chromium
(TLS reset), so we solve the PoW directly: SHA-256(randomData + nonce) must have
`difficulty` leading zero nibbles, then GET pass-challenge to obtain the auth
cookie. One solve yields a domain-wide auth cookie reused for every file.
"""
import hashlib
import json
import os
import re
import sys
import time

import requests

# The original Linux agent container reached the network through a local MITM proxy with
# its own CA bundle. Off that container (e.g. a Windows workstation) there is no proxy and
# the system trust store is correct, so both default to unset and are env-overridable.
PROXY = os.environ.get("ASB_HTTP_PROXY") or None
CA = os.environ.get("ASB_CA_BUNDLE") or None
UA = ("Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 "
      "(KHTML, like Gecko) Chrome/125.0.0.0 Safari/537.36")
REFERER = "https://datadryad.org/dataset/doi:10.5061/dryad.kkwh70sg0"
BASE = "https://datadryad.org"
CHALLENGE_RE = re.compile(r'id="anubis_challenge"[^>]*>(.*?)</script>', re.S)


def _session():
    s = requests.Session()
    if PROXY:
        s.proxies = {"http": PROXY, "https": PROXY}
    if CA:
        s.verify = CA
    s.headers.update({"User-Agent": UA, "Referer": REFERER,
                      "Accept": "*/*", "Accept-Language": "en-US,en;q=0.9"})
    return s


def _solve(challenge_json):
    ch = json.loads(challenge_json)
    rd = ch["challenge"]["randomData"]
    diff = ch["rules"]["difficulty"]
    cid = ch["challenge"]["id"]
    target = "0" * diff
    t0 = time.time()
    nonce = 0
    while True:
        h = hashlib.sha256((rd + str(nonce)).encode()).hexdigest()
        if h.startswith(target):
            return cid, h, nonce, int((time.time() - t0) * 1000) or 1
        nonce += 1


def _is_challenge(resp):
    ctype = resp.headers.get("Content-Type", "")
    return ("text/html" in ctype and b"anubis_challenge" in resp.content[:8000])


def clear_challenge(s, url):
    """GET url; if Anubis-challenged, solve and pass, leaving auth cookie in s."""
    r = s.get(url, allow_redirects=True, timeout=90)
    if not _is_challenge(r):
        return r  # already cleared or direct content
    m = CHALLENGE_RE.search(r.text)
    if not m:
        raise RuntimeError("challenge page without anubis_challenge block")
    cid, resp_hash, nonce, elapsed = _solve(m.group(1))
    print(f"[anubis] solved id={cid[:12]} nonce={nonce} hash={resp_hash[:12]} "
          f"elapsed={elapsed}ms", flush=True)
    pc = s.get(f"{BASE}/.within.website/x/cmd/anubis/api/pass-challenge",
               params={"id": cid, "response": resp_hash, "nonce": nonce,
                       "redir": url, "elapsedTime": elapsed},
               allow_redirects=False, timeout=60)
    print(f"[anubis] pass-challenge status={pc.status_code} "
          f"set-cookie={'techaro' in pc.headers.get('Set-Cookie', '')}", flush=True)
    return s.get(url, allow_redirects=True, timeout=90)


def download(s, url, dest):
    r = s.get(url, stream=True, allow_redirects=True, timeout=180)
    if _is_challenge(r):                       # cookie expired mid-batch -> re-clear
        r.close()
        clear_challenge(s, url)
        r = s.get(url, stream=True, allow_redirects=True, timeout=180)
        if _is_challenge(r):
            raise RuntimeError(f"still challenged for {url}")
    total = 0
    with open(dest, "wb") as fh:
        for chunk in r.iter_content(chunk_size=1 << 20):
            if chunk:
                fh.write(chunk)
                total += len(chunk)
    return r.status_code, total, r.headers.get("Content-Type", "")


def _load_pairs():
    if len(sys.argv) >= 3 and sys.argv[1] == "--from-file":
        pairs = []
        for line in open(sys.argv[2]):
            line = line.strip()
            if line:
                fid, dest = line.split(None, 1)
                pairs.append((fid, dest))
        return pairs
    return list(zip(sys.argv[1::2], sys.argv[2::2]))


if __name__ == "__main__":
    import os
    s = _session()
    pairs = _load_pairs()
    clear_challenge(s, f"{BASE}/downloads/file_stream/{pairs[0][0]}")
    done = 0
    for i, (fid, dest) in enumerate(pairs, 1):
        os.makedirs(os.path.dirname(dest) or ".", exist_ok=True)
        if os.path.exists(dest) and os.path.getsize(dest) > 1024:
            print(f"[skip] {dest} exists ({os.path.getsize(dest):,} B)", flush=True)
            done += 1
            continue
        url = f"{BASE}/downloads/file_stream/{fid}"
        code, n, ct = download(s, url, dest)
        done += 1
        print(f"[dl {i}/{len(pairs)}] {fid} -> {dest}  status={code} "
              f"bytes={n:,} ctype={ct}", flush=True)
    print(f"FETCH_DONE ({done}/{len(pairs)})", flush=True)
