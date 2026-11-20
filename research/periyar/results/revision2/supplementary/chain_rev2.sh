#!/bin/bash
# Second-round revision runs, one stage after another; each stage resumes if restarted.
cd /home/user/hydrolib-core
R="python3 -m research.periyar.sensitivity /home/user/ext/periyar_v5"
C="--threads 1 --dflowfm /home/user/ext/run_dflowfm.sh --sar /home/user/ext/sar_flood_20180821_v3.tif"
E=/home/user/ext/rev2
$R $E/maps --method file --samples-file $E/samples_maps.csv --keep-maps --parallel 2 $C >> $E/maps.log 2>&1
$R $E/multiplicative --method file --samples-file $E/samples_accepted.csv --release-model multiplicative --parallel 4 $C >> $E/multiplicative.log 2>&1
$R $E/fixed --method file --samples-file $E/samples_fixed.csv --parallel 4 $C >> $E/fixed.log 2>&1
$R $E/lhs8 --method lhs --samples 200 --seed 44 --parallel 4 $C >> $E/lhs8.log 2>&1
$R $E/morris6 --method morris --trajectories 20 --factors n_channel n_floodplain channel_depth canopy_factor q_rating tide_offset --parallel 4 $C >> $E/morris6.log 2>&1
echo done > $E/ALL_DONE
