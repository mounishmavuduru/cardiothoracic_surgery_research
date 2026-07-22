#!/bin/bash
# Self-healing supervisor for the 100k neural scale-up.
#
# Runs the (resumable, sharded) builder in a loop: if the builder process dies
# for ANY reason before the run is complete, relaunch it — it skips shards
# already on disk, so it picks up where it left off, losing at most the one
# in-progress 2500-network shard (~20 min). Every launch/exit is timestamped in
# the supervisor log, which doubles as diagnostic evidence: if the builder keeps
# dying at a consistent elapsed time, that reveals a periodic reaper.
#
# Launch fully detached so it is not torn down with the shell that started it:
#   setsid nohup bash scripts/supervise_100k.sh >/dev/null 2>&1 </dev/null &
set -u
cd /home/user/cardiothoracic_surgery_research
LOG=/tmp/run_100k.log
SUP=/tmp/run_100k_sup.log
echo "supervisor START $(date -u +%FT%TZ) pid=$$ sid=$(ps -o sid= -p $$ | tr -d ' ')" >> "$SUP"
for attempt in $(seq 1 500); do
  if grep -q "AGGREGATE_DONE" "$LOG" 2>/dev/null; then
    echo "supervisor: AGGREGATE_DONE seen -> exit $(date -u +%FT%TZ)" >> "$SUP"
    break
  fi
  ns=$(ls outputs/scaled100k/shard_*.json 2>/dev/null | wc -l)
  echo "supervisor: launch #$attempt $(date -u +%FT%TZ) shards_on_disk=$ns" >> "$SUP"
  .venv/bin/python -u scripts/run_gm4_100k.py >> "$LOG" 2>&1
  echo "supervisor: builder exited code=$? $(date -u +%FT%TZ)" >> "$SUP"
  sleep 3
done
echo "supervisor END $(date -u +%FT%TZ)" >> "$SUP"
