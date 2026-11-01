#!/bin/bash
# Round-5 runs: forcing variants for 14 stratified members of the 50 runs accepted only at +-1.0 m.
cd /home/user/hydrolib-core
E=/home/user/ext/rev5
R="python3 -m research.periyar.sensitivity /home/user/ext/periyar_v5"
C="--threads 1 --dflowfm /home/user/ext/run_dflowfm.sh --sar /home/user/ext/sar_flood_20180821_v3.tif"
$R $E/multiplicative --method file --samples-file $E/samples_central.csv --release-model multiplicative --parallel 4 $C >> $E/multiplicative.log 2>&1
$R $E/lag_short --method file --samples-file $E/samples_central.csv --lags 4 2 --parallel 4 $C >> $E/lag_short.log 2>&1
$R $E/lag_long --method file --samples-file $E/samples_central.csv --lags 12 6 --parallel 4 $C >> $E/lag_long.log 2>&1
echo done > $E/ALL_DONE
