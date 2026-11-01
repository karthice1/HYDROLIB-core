# Uploading to Zenodo

1. Log in at https://zenodo.org (a GitHub or ORCID login works) and click **New upload**.
2. Drag `periyar-flood-workflow-v1.0.zip` into the Files box. Do not unzip it. Optionally check
   it first: its SHA-256 must match `periyar-flood-workflow-v1.0.zip.sha256`.
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
| Version | v1.0-periyar |
| Publication date | the upload date |
| Licence | MIT License (the description explains the two copyright holders; see NOTICE in the archive) |
| Keywords | flood inundation modelling; global sensitivity analysis; Morris method; uncertainty; D-Flow FM; HYDROLIB-core; Sentinel-1; Kerala 2018 flood |
| Related works | "Is supplement to" https://github.com/karthice1/HYDROLIB-core/tree/v1.0-periyar (Software); "Is derived from" https://github.com/Deltares/HYDROLIB-core (Software) |

Description (paste as is):

Research code, per-run results and manuscript sources for the article "How far can one gauge mark and one radar image narrow parameter uncertainty in an open-data flood model? Sensitivity analysis of the 2018 Periyar flood". The archive is a snapshot of the GitHub fork karthice1/HYDROLIB-core at tag v1.0-periyar, plus the clean-environment test report of that exact commit. Licensing (see NOTICE): HYDROLIB-core 1.5.0, copyright Deltares, MIT (LICENSE), included unchanged; the research code in research/gsa and research/periyar, copyright Karthik V and Roopa D, MIT (research/LICENSE); the manuscript source in research/periyar/manuscript is not MIT-licensed. python -m research.periyar.check_results recomputes 120 published values from the stored results without running D-Flow FM. See README_FIRST.md and research/periyar/REPRODUCE.md.
