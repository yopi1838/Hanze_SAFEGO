#!/usr/bin/env bash
# =====================================================================
# monitor_3dec_cpu.sh
# Log CPU / memory of running processes (default: 3DEC) over time.
#
# Usage:
#   ./monitor_3dec_cpu.sh [pattern] [interval_s] [logfile]
#     pattern    process-name/command match   (default: 3dec)
#     interval_s seconds between samples       (default: 15)
#     logfile    CSV output path               (default: cpu_monitor_<ts>.csv)
#
# CPU% is the TRUE instantaneous usage over each interval, computed from
# /proc (utime+stime deltas), so a process using N cores reads ~N*100%.
# On a multi-core box compare the TOTAL against (cores * 100%).
# Stop with Ctrl-C; the CSV is kept.
#
# Tips:
#   - Find the exact 3DEC binary name first:   pgrep -af 3dec
#     then pass it as the pattern for a tight match, e.g.
#       ./monitor_3dec_cpu.sh 3dec_dp 30 cpu_3runs.csv
#   - Run inside tmux/screen (or with nohup) so it survives an SSH drop.
# =====================================================================
set -u

PATTERN="${1:-3dec}"
INTERVAL="${2:-15}"
LOG="${3:-cpu_monitor_$(date +%Y%m%d_%H%M%S).csv}"

NCORES=$(nproc 2>/dev/null || echo 1)
CLK=$(getconf CLK_TCK 2>/dev/null || echo 100)
PAGE=$(getconf PAGESIZE 2>/dev/null || echo 4096)
SELF=$$

declare -A PREV_J          # previous cpu jiffies per pid
PREV_UP=0

echo "timestamp,pid,cpu_pct,mem_pct,rss_mb,threads,elapsed,command" > "$LOG"
echo "Monitoring '$PATTERN' every ${INTERVAL}s | cores=$NCORES | log=$LOG"
echo "(Ctrl-C to stop)"

cleanup(){ echo; echo "Stopped. Log saved to: $LOG"; exit 0; }
trap cleanup INT TERM

while true; do
    TS=$(date '+%Y-%m-%d %H:%M:%S')
    UP=$(awk '{print $1}' /proc/uptime)
    DT=$(awk -v a="$UP" -v b="$PREV_UP" 'BEGIN{d=a-b; print (d>0)?d:0}')
    LOAD=$(awk '{printf "%s/%s/%s", $1,$2,$3}' /proc/loadavg)

    # PIDs matching pattern, excluding this script itself
    PIDS=$(pgrep -f "$PATTERN" 2>/dev/null | grep -vw "$SELF" || true)

    printf '\n[%s] load(1/5/15)=%s  cores=%s\n' "$TS" "$LOAD" "$NCORES"
    if [ -z "$PIDS" ]; then
        printf '   (no process matches "%s")\n' "$PATTERN"
    else
        printf '   %-7s %8s %6s %9s %7s %11s  %s\n' PID CPU% MEM% RSS_MB THR ELAPSED CMD
        TOTAL=0
        for pid in $PIDS; do
            [ -r "/proc/$pid/stat" ] || continue
            st=$(cat /proc/$pid/stat 2>/dev/null) || continue
            rest=${st#*) }; set -- $rest
            utime=${12}; stime=${13}; nthreads=${18}
            jif=$((utime + stime))
            cpu="..."
            if [ -n "${PREV_J[$pid]:-}" ] && [ "$(awk -v d="$DT" 'BEGIN{print (d>0)?1:0}')" = "1" ]; then
                dj=$(( jif - PREV_J[$pid] ))
                cpu=$(awk -v dj="$dj" -v clk="$CLK" -v dt="$DT" 'BEGIN{printf "%.1f", 100.0*dj/(clk*dt)}')
                TOTAL=$(awk -v t="$TOTAL" -v c="$cpu" 'BEGIN{printf "%.1f", t+c}')
            fi
            PREV_J[$pid]=$jif
            mem=$(ps -o pmem= -p "$pid" 2>/dev/null | tr -d ' '); mem=${mem:-0}
            rss_pages=$(awk '{print $2}' /proc/$pid/statm 2>/dev/null); rss_pages=${rss_pages:-0}
            rss_mb=$(awk -v p="$rss_pages" -v pg="$PAGE" 'BEGIN{printf "%.0f", p*pg/1048576}')
            etime=$(ps -o etime= -p "$pid" 2>/dev/null | tr -d ' '); etime=${etime:-?}
            comm=$(cat /proc/$pid/comm 2>/dev/null); comm=${comm:-?}
            printf '   %-7s %8s %6s %9s %7s %11s  %s\n' "$pid" "$cpu" "$mem" "$rss_mb" "$nthreads" "$etime" "$comm"
            echo "$TS,$pid,$cpu,$mem,$rss_mb,$nthreads,$etime,$comm" >> "$LOG"
        done
        [ "$TOTAL" != "0" ] && printf '   %-7s %8s   (%.0f%% of %s cores)\n' "TOTAL" "$TOTAL" \
            "$(awk -v t="$TOTAL" -v n="$NCORES" 'BEGIN{print 100*t/(n*100)}')" "$NCORES"
    fi
    PREV_UP=$UP
    sleep "$INTERVAL"
done
