# DAVE1 VKGL VUS Pathogenicity Viewer

This Streamlit app provides an interactive interface for exploring DAVE1 pathogenicity predictions of Variants of Uncertain Significance (VUS) from the VKGL consensus datasharing. It integrates data visualization, protein structure modeling, and molecular feature analysis to support variant interpretation in a clinical genomics context.

Our article is now in preprint! DOI: xxxx.xx.xx.xx

---

## Features

- 🔍 Searchable and filterable VUS table from VKGL dataset
- 📊 Force plot visualization of DAVE1 LP score contributions
- 🧬 Protein structure modeling with FoldX
- 🧪 3D visualization of wild-type and mutant proteins using Py3Dmol

---

## File Requirements

Ensure the following files and directories are present in the working directory:

- `vkgl_apr2024_VUS_pred.csv`: VUS prediction data downloaded from the dave1 resources.
- `mut_wt_structures_vkgl_vus.tar.gz`: AlphaFold wild-type and mutant PDB files: 

---

## Installation

Get the github repo:
```bash
git clone https://github.com/molgenis/dave1_streamlit
```

Install required Python packages (Tip: use a Python virtual environment or Conda environment):

```bash
pip install streamlit matplotlib pandas numpy py3Dmol plotly

```

## Running the streamlit server

```bash

cd /location/of/dave1_streamlit

streamlit run main.py

# Optionally change the port with --server.port

```

## Example NKX2-5 L153P:

![table](images/table.png)

![force_plot](images/explain_plot.png)

![force_plot](images/py3dmol.png)