# Values for the Zenodo "New upload" form

| Field | Value |
|---|---|
| Resource type | Software |
| Title | Periyar 2018 flood: open-data D-Flow FM model, sensitivity and uncertainty workflow (research code and results) |
| Creators | V, Karthik (Karpagam Academy of Higher Education, Coimbatore, India); D, Roopa (Karpagam Academy of Higher Education, Coimbatore, India) |
| Version | v1.0 |
| Publication date | the upload date |
| Licence | MIT License |
| Keywords | flood inundation modelling; global sensitivity analysis; Morris method; uncertainty; D-Flow FM; HYDROLIB-core; Sentinel-1; Kerala 2018 flood |
| Related works | "Is derived from" https://github.com/Deltares/HYDROLIB-core (Software); "Is supplement to" https://github.com/karthice1/HYDROLIB-core (Software) |

Description (paste as is):

Research code, per-run results and manuscript sources for the article "How far can one
gauge mark and one radar image narrow parameter uncertainty in an open-data flood model?
Sensitivity analysis of the 2018 Periyar flood". The study code is in research/periyar and
research/gsa. It builds on the HYDROLIB-core library by Deltares (MIT licence), which the
archive contains unchanged. python -m research.periyar.check_results recomputes the
principal numbers of the paper from the stored results. See README.md and
research/periyar/REPRODUCE.md for installation, reproduction and reuse.
