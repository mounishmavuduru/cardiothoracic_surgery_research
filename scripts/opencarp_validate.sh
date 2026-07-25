#!/usr/bin/env bash
# Validate that openCARP -- the PRE-REGISTERED ground-truth solver, substituted in
# section 7.1 because it would not install on the original container -- runs on this
# machine, with the same membrane model as the in-house solver it would replace.
#
# Run from WSL (Ubuntu). Needs no root: the AppImage is extracted rather than mounted,
# because libfuse.so.2 is absent even though fusermount is present.
#
#   bash scripts/opencarp_validate.sh
#
# Established 2026-07-25 (see PRE_REGISTRATION.md amendment 8.6 (g)):
#   * binary runs                     -> openCARP -buildinfo returns cleanly
#   * MitchellSchaeffer is available  -> the same membrane model as asb.labels.monodomain
#   * single cell                     -> APD90 = 246 ms (in-house solver gives 218 ms)
#   * tissue                          -> planar wave sweeps a 10,201-node slab
set -euo pipefail

VER=v19.0
ROOT="${OPENCARP_ROOT:-$HOME/opencarp}"
APP="$ROOT/openCARP-${VER}-x86_64_AppImage"
SQ="$APP/squashfs-root"
URL="https://git.opencarp.org/api/v4/projects/16/packages/generic/opencarp-appimage/${VER}/openCARP-${VER}-x86_64_AppImage.tar.gz"

mkdir -p "$ROOT" && cd "$ROOT"
if [ ! -d "$SQ" ]; then
  [ -f oc.tar.gz ] || curl -sL -o oc.tar.gz "$URL"
  tar -xzf oc.tar.gz
  cd "$APP" && chmod +x "openCARP-${VER}-x86_64.AppImage"
  # FUSE mount fails (libfuse.so.2 missing); extract instead. No root required.
  "./openCARP-${VER}-x86_64.AppImage" --appimage-extract >/dev/null
fi
export LD_LIBRARY_PATH="$SQ/usr/lib:${LD_LIBRARY_PATH:-}"
BIN="$SQ/usr/bin"

echo "== 1. binary =="
"$BIN/openCARP" -buildinfo 2>&1 | head -2

echo "== 2. membrane model available =="
"$BIN/bench" --list-imps 2>&1 | tr ',' '\n' | grep -c MitchellSchaeffer \
  | xargs -I{} echo "   MitchellSchaeffer present: {}"

echo "== 3. single cell =="
rm -rf /tmp/oc_bench && mkdir -p /tmp/oc_bench && cd /tmp/oc_bench
"$BIN/bench" --imp=MitchellSchaeffer --duration=500 --stim-start=10 --stim-dur=2 \
  --stim-curr=60 --dt=0.01 --fout=ms >/dev/null 2>&1
awk 'NR>1{t=$1;v=$2;if(n==0){mx=v}if(v>mx)mx=v;T[++n]=t;V[n]=v}
     END{thr=mx*0.1;for(i=1;i<=n;i++){if(up==0&&V[i]>=thr)up=T[i];
     if(up>0&&V[i]<thr){printf "   APD90 = %.1f ms\n",T[i]-up;exit}}}' ms.txt

echo "== 4. tissue propagation =="
rm -rf /tmp/oc_tissue && mkdir -p /tmp/oc_tissue && cd /tmp/oc_tissue
"$BIN/mesher" -size[0] 1.0 -size[1] 1.0 -size[2] 0.0 -bath[0] 0 -bath[1] 0 -bath[2] 0 \
  -center[0] 0 -center[1] 0 -center[2] 0 -resolution[0] 100 -resolution[1] 100 \
  -resolution[2] 100 -mesh block -Elem3D 0 >/dev/null 2>&1
cat > run.par <<'PAR'
num_phys_regions = 1
phys_region[0].name = "intra"
phys_region[0].ptype = 0
phys_region[0].num_IDs = 1
phys_region[0].ID[0] = 1
num_imp_regions = 1
imp_region[0].im = MitchellSchaeffer
imp_region[0].num_IDs = 1
imp_region[0].ID[0] = 1
num_gregions = 1
gregion[0].num_IDs = 1
gregion[0].ID[0] = 1
gregion[0].g_il = 0.174
gregion[0].g_it = 0.019
gregion[0].g_in = 0.019
num_stim = 1
stimulus[0].name = "S1"
stimulus[0].stimtype = 0
# NB MitchellSchaeffer is a NORMALISED model (Vm ~ 0..1). A physiological
# 250 uA/cm^2 stimulus makes the parabolic solve diverge with NaN; 60 works.
stimulus[0].strength = 60.0
stimulus[0].duration = 2.0
stimulus[0].start = 1.0
stimulus[0].x0 = -600
stimulus[0].xd = 300
stimulus[0].y0 = -600
stimulus[0].yd = 1200
stimulus[0].z0 = -100
stimulus[0].zd = 200
bidomain = 0
tend = 60.0
dt = 10
spacedt = 5
timedt = 5
PAR
"$BIN/openCARP" +F run.par -meshname block -simID out >/dev/null 2>&1
python3 - <<'PY'
import struct, re
d = open("out/vm.igb", "rb").read()
h = d[:1024].decode("latin-1")
x = int(re.search(r"x:(\d+)", h).group(1)); t = int(re.search(r"t:(\d+)", h).group(1))
v = struct.unpack("<%df" % (x * t), d[1024:1024 + x * t * 4])
act = [sum(1 for u in v[i * x:(i + 1) * x] if u > 0.5) for i in range(t)]
print(f"   activated nodes per frame ({x} total): {act}")
print("   VERDICT:", "PROPAGATION CONFIRMED" if max(act) > 100 and max(act) > act[0]
      else "NO PROPAGATION")
PY
