# Observed August 2018 data (digitised)

Digitised from the vector charts in Central Water Commission, Hydrological Studies
Organisation, *Study Report: Kerala Floods of August 2018* (September 2018), by reading
the plotted coordinates from the PDF and calibrating them against the chart gridlines.

| File | Source figure | Resolution | Check against report text |
|---|---|---|---|
| `neeleswaram_discharge_cwc_fig3.csv` | Fig. 3, Periyar at Neeleswaram G&D site | 6-hourly, 14–20 Aug | Peak 8,803 vs 8,800 m³/s on 16 Aug 14:00 |
| `idukki_cwc_fig4.csv` | Fig. 4, Idukki reservoir inflow, outflow (spill to Periyar) and level | hourly, 9–18 Aug | Spill plateau 1,500 m³/s, FRL reached 17 Aug ~01:00 |
| `idamalayar_cwc_fig5.csv` | Fig. 5, Idamalayar reservoir inflow, outflow and level | daily, 1 Jul–23 Aug | 16 Aug inflow 1,166 vs 1,164; release 964 vs 963; 17 Aug release 1,275 vs 1,272 |

Times are IST. The hourly Idukki inflow is smoothed by hourly averaging (peak 2,467
vs 2,532 m³/s in the text), so use the outflow series, which matches the text closely.
Treat all values as digitised estimates (about ±1 hour, ±2–3 %).
