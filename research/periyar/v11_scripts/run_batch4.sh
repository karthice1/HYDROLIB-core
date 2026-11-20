#!/bin/bash
# the 294 runs with |HFL error| > 2 m (author decision 18 Nov 2026), for radar-only statistics
until [ -f /home/user/ext/rerun_radar/FINISHED3 ]; do sleep 60; done
cd /home/user/ext/ema_analysis
C=/home/user/ext/sar_flood_20180821_corrected.tif
python -I rerun_for_radar.py /home/user/ext/gsa_v5_lhs200 /home/user/ext/rerun_radar/lhs200 $C $(cat rest_lhs200.txt) > /home/user/ext/rerun_radar_b4.log 2>&1
python -I rerun_for_radar.py /home/user/ext/gsa_v5_lhs400 /home/user/ext/rerun_radar/lhs400 $C $(cat rest_lhs400.txt) >> /home/user/ext/rerun_radar_b4.log 2>&1
echo FINISHED > /home/user/ext/rerun_radar/FINISHED4
