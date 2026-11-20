#!/bin/bash
# forcing tests (release formulation, travel times, repeats) of the previously accepted runs; waits for batch 2
until [ -f /home/user/ext/rerun_radar/FINISHED2 ]; do sleep 60; done
cd /home/user/ext/ema_analysis
C=/home/user/ext/sar_flood_20180821_corrected.tif
ids() { python3 -c "import pandas as pd,sys; print(' '.join(map(str,pd.read_csv(sys.argv[1],index_col=0).index)))" $1; }
python -I rerun_for_radar.py /home/user/hydrolib-core/research/periyar/results/revision2/supplementary/multiplicative /home/user/ext/rerun_radar/mult_strict $C $(ids /home/user/hydrolib-core/research/periyar/results/revision2/supplementary/multiplicative/samples.csv) --release-model multiplicative > /home/user/ext/rerun_radar_b3.log 2>&1
python -I rerun_for_radar.py /home/user/hydrolib-core/research/periyar/results/revision2/supplementary/lag_short /home/user/ext/rerun_radar/lagshort_strict $C $(ids /home/user/hydrolib-core/research/periyar/results/revision2/supplementary/lag_short/samples.csv) --lags 4 2 >> /home/user/ext/rerun_radar_b3.log 2>&1
python -I rerun_for_radar.py /home/user/hydrolib-core/research/periyar/results/revision2/supplementary/lag_long /home/user/ext/rerun_radar/laglong_strict $C $(ids /home/user/hydrolib-core/research/periyar/results/revision2/supplementary/lag_long/samples.csv) --lags 12 6 >> /home/user/ext/rerun_radar_b3.log 2>&1
python -I rerun_for_radar.py /home/user/hydrolib-core/research/periyar/results/revision2/central/multiplicative /home/user/ext/rerun_radar/mult_central $C $(ids /home/user/hydrolib-core/research/periyar/results/revision2/central/multiplicative/samples.csv) --release-model multiplicative >> /home/user/ext/rerun_radar_b3.log 2>&1
python -I rerun_for_radar.py /home/user/hydrolib-core/research/periyar/results/revision2/central/lag_short /home/user/ext/rerun_radar/lagshort_central $C $(ids /home/user/hydrolib-core/research/periyar/results/revision2/central/lag_short/samples.csv) --lags 4 2 >> /home/user/ext/rerun_radar_b3.log 2>&1
python -I rerun_for_radar.py /home/user/hydrolib-core/research/periyar/results/revision2/central/lag_long /home/user/ext/rerun_radar/laglong_central $C $(ids /home/user/hydrolib-core/research/periyar/results/revision2/central/lag_long/samples.csv) --lags 12 6 >> /home/user/ext/rerun_radar_b3.log 2>&1
python -I rerun_for_radar.py /home/user/ext/rev2/repeat_additive /home/user/ext/rerun_radar/repeat_add $C $(ids /home/user/ext/rev2/repeat_additive/samples.csv) >> /home/user/ext/rerun_radar_b3.log 2>&1
python -I rerun_for_radar.py /home/user/ext/rev2/repeat_multiplicative /home/user/ext/rerun_radar/repeat_mult $C $(ids /home/user/ext/rev2/repeat_multiplicative/samples.csv) --release-model multiplicative >> /home/user/ext/rerun_radar_b3.log 2>&1
echo FINISHED > /home/user/ext/rerun_radar/FINISHED3
