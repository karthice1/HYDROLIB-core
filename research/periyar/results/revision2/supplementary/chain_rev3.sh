#!/bin/bash
# Third-round runs: release travel-time lags for the 23 strict-accepted runs, after chain_rev2.
cd /home/user/hydrolib-core
E=/home/user/ext/rev2
until [ -f $E/ALL_DONE ]; do sleep 60; done
R="python3 -m research.periyar.sensitivity /home/user/ext/periyar_v5"
C="--threads 1 --dflowfm /home/user/ext/run_dflowfm.sh --sar /home/user/ext/sar_flood_20180821_v3.tif"
$R $E/lag_short --method file --samples-file $E/samples_accepted.csv --lags 4 2 --parallel 4 $C >> $E/lag_short.log 2>&1
$R $E/lag_long --method file --samples-file $E/samples_accepted.csv --lags 12 6 --parallel 4 $C >> $E/lag_long.log 2>&1
echo done > $E/ALL_DONE3
