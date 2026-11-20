#!/bin/bash
# Repeat check: three accepted runs whose HFL changed most under the multiplicative form,
# rerun in both forms with the alternative gauge operators recorded. Runs after chain_rev3.
cd /home/user/hydrolib-core
E=/home/user/ext/rev2
until [ -f $E/ALL_DONE3 ]; do sleep 60; done
R="python3 -m research.periyar.sensitivity /home/user/ext/periyar_v5"
C="--threads 1 --dflowfm /home/user/ext/run_dflowfm.sh --sar /home/user/ext/sar_flood_20180821_v3.tif"
$R $E/repeat_additive --method file --samples-file $E/samples_repeat.csv --parallel 3 $C >> $E/repeat.log 2>&1
$R $E/repeat_multiplicative --method file --samples-file $E/samples_repeat.csv --release-model multiplicative --parallel 3 $C >> $E/repeat.log 2>&1
echo done > $E/ALL_DONE4
