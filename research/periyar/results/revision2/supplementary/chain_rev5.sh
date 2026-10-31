#!/bin/bash
# Round-4 review: second, independently seeded eight-factor ensemble (seed 46). Resumes if restarted.
cd /home/user/hydrolib-core
R="python3 -m research.periyar.sensitivity /home/user/ext/periyar_v5"
C="--threads 1 --dflowfm /home/user/ext/run_dflowfm.sh --sar /home/user/ext/sar_flood_20180821_v3.tif"
E=/home/user/ext/rev2
$R $E/lhs8b --method lhs --samples 200 --seed 46 --parallel 4 $C >> $E/lhs8b.log 2>&1
echo done > $E/ALL_DONE5
