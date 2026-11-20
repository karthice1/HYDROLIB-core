#!/bin/bash
# waits for batch 1, then reruns the eight-factor and fixed-factor runs that pass the HFL at +-1.0 m
until [ -f /home/user/ext/rerun_radar/FINISHED ]; do sleep 60; done
cd /home/user/ext/ema_analysis
for d in lhs8 lhs8b fixed; do
  python -I rerun_for_radar.py /home/user/ext/rev2/$d /home/user/ext/rerun_radar/$d /home/user/ext/sar_flood_20180821_corrected.tif $(cat need_$d.txt) --parallel 4 > /home/user/ext/rerun_radar_$d.log 2>&1
done
echo FINISHED > /home/user/ext/rerun_radar/FINISHED2
