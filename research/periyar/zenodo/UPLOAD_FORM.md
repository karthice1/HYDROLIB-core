# Uploading to Zenodo

1. Log in at https://zenodo.org (a GitHub or ORCID login works) and click **New upload**.
2. Drag `periyar-flood-workflow-v1.0.zip` into the Files box. Do not unzip it.
3. Under **Digital Object Identifier**, answer "No" to "Do you already have a DOI" and click
   **Get a DOI now!**. This reserves the DOI before publishing.
4. Fill in the form with the values below, then click **Publish**. After publishing, the
   files cannot be changed (only a new version can be added).
5. Copy the DOI (10.5281/zenodo.NNNNNNN) from the record page and send it.

## Form values

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
