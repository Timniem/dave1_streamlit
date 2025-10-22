import os
import subprocess
import gzip
import shutil
import glob
import base64
from pathlib import Path


import streamlit as st
import numpy as np
import pandas as pd
import py3Dmol

from explain_plot import explain_plot_plotly

# img_to_bytes and img_to_html inspired from https://pmbaumgartner.github.io/streamlitopedia/sizing-and-images.html
def img_to_bytes(img_path):
    img_bytes = Path(img_path).read_bytes()
    encoded = base64.b64encode(img_bytes).decode()
    return encoded
def img_to_html(img_path, align=None, size=100):
    if align:
        img_html = f"<img src='data:image/png;base64,{img_to_bytes(img_path)}' align={align} style='width:{size}%' class='img-fluid'>"
    else:
        img_html = f"<img src='data:image/png;base64,{img_to_bytes(img_path)}' style='width:{size}%' class='img-fluid'>"
    return img_html

# Define a callback function to disable the button
def disable_button():
    st.session_state.button_disabled = True

###CONST###
FEATURE_NAMES_DICT = {"delta_DNAs_cumu_bin":"DNA binding site (Δ residues) ","delta_RNAs_cumu_bin":"RNA binding site (Δ residues)",
                 "delta_ProtS_cumu_bin":"Protein binding site (Δ residues)","delta_DNAb_binding_affinity_pKd":"DNA binding affinity (pKd)",
                 "delta_RNAb_binding_affinity_pKd":"RNA binding affinity (pKd)","delta_ligand_nr_of_predicted_pockets":"Ligand binding № pockets",
                 "delta_ligand_rank1_sas_points":"Ligand binding top pocket (SAS points)","delta_charge":"Net electric charge at pH 7 (pKa)",
                 "delta_hydrophobicMoment":"Hydrophobic moment","delta_hydrophobicity":"Hydrophobicity (GRAVY)",
                 "delta_isoElecPoint":"Isoelectric point (pH)", "delta_total.energy":"Folding energy ΔG (kJ/mol)"}
###


### init ###
if 'button_disabled' not in st.session_state:
    st.session_state.button_disabled = False
###


features_shap = [f"{key}.sph" for key in FEATURE_NAMES_DICT.keys()]

st.markdown(f"<div style='background-color: #017FFD; padding:10px; display: flex; justify-content: space-between; align-items: center;'> {img_to_html('images/logo_blue.png', size=20)} {img_to_html('images/umcg_logo.png', align='right', size=19)} </div>", unsafe_allow_html=True, width="stretch")
# Title
st.markdown("<h3 style='text-align: center; color: #666;'>DAVE1 scores VKGL Datasharing VUS</h3>", unsafe_allow_html=True)
st.markdown("<p style='text-align: center; font-size: small; color: #B0B0B0;'>Pathogenicity prediction by the DAVE1 model, for more information: </p>", unsafe_allow_html=True)
# Load VUS data
vus_path = "vkgl_apr2024_VUS_pred.csv"
if not os.path.exists(vus_path):
    st.error("CSV file not found. Please check the path.")
    st.stop()


vkgl_consensus_vus = pd.read_csv(vus_path).sort_values(by="LP", ascending=False)
vkgl_consensus_vus["AA change"] = vkgl_consensus_vus["delta_aaSeq"].apply(lambda x: f'{x[0]}{x[2:]}')

c1, c2, c3 = st.columns(3)

with c2:
    search_term = st.text_input(" ").lower()

if search_term:
    mask = vkgl_consensus_vus.astype(str).apply(lambda col: col.str.lower().str.contains(search_term))
    filtered_data = vkgl_consensus_vus[mask.any(axis=1)]
else:
    filtered_data = vkgl_consensus_vus

edited_df = st.dataframe(
    filtered_data[["LP","gene","dna_variant_chrom","dna_variant_pos","dna_variant_ref","dna_variant_alt","AA change","TranscriptID"]].rename(columns={"LP":"DAVE1 LP score",
                                                                                                                                                            "dna_variant_chrom":"chrom",
                                                                                                                                                            "dna_variant_pos":"pos",
                                                                                                                                                            "dna_variant_ref":"ref",
                                                                                                                                                            "dna_variant_alt":"alt",
                                                                                                                                                            "delta_aaSeq":"AAchange"}),
    width="stretch",
    hide_index=True,
    on_select="rerun",
    selection_mode="single-row",
)

selected_rows = edited_df.selection.rows

if len(selected_rows) > 0:
    selected_df = filtered_data.iloc[selected_rows[0]]
    
    selected_prot = selected_df['UniProtID']
    gene_symbol = selected_df['gene']
    prot_change = selected_df['AA change']
    selected_variant = selected_df['delta_aaSeq']
    localization = selected_df['ann_proteinLocalization']
    if type(selected_df['seqFt']) != float:
        feature = selected_df['seqFt'].split("|")[-1].split("~")
    else:
        feature = None

    # Input for FoldX and PDB
    foldx_executable =  "./foldx/foldx5_1Mac/foldx_20251231"
    pdb_file_zipped = f'./AF_pdb/AF-{selected_prot}-F1-model_v4.pdb.gz'

    selected_df = selected_df.rename(FEATURE_NAMES_DICT)

    st.markdown(f"<p style='text-align: center; color: #333;'>Feature contribution: <i>{gene_symbol}</i> {prot_change} </p>", unsafe_allow_html=True)
    # Old matplotlib style: st.pyplot(force_plot(selected_df[features_shap], selected_df[FEATURE_NAMES_DICT.values()], FEATURE_NAMES_DICT.values(), selected_df["LP"]))

    fig = explain_plot_plotly(selected_df[features_shap], selected_df[FEATURE_NAMES_DICT.values()], FEATURE_NAMES_DICT.values(), selected_df["LP"])
    st.plotly_chart(fig, use_container_width=True, config={"displayModeBar": False, "scrollZoom": False,"doubleClick": False })

    st.markdown(f"<p style='text-align: left; color: #555;'>Residue information</p>", unsafe_allow_html=True)
    st.markdown(f"<p style='text-align: left; font-size: small; color: #999;'>Localization: {localization}</p>", unsafe_allow_html=True)
    if feature:
        st.markdown(f"<p style='text-align: left; font-size: small; color: #999;'>Sequence features: {' -> '.join(feature)}</p>", unsafe_allow_html=True)
    

    c1, c2, c3 = st.columns(3)

    with c2:
        viz_button = st.button("Visualize Variant", width="stretch", on_click=disable_button, disabled=st.session_state.button_disabled )

    if viz_button:
        # Generate mutant file content
        mutant_string = f"{selected_variant};"
        mutant_file = "individual_list.txt"
        with open(mutant_file, "w") as f:
            f.write(mutant_string)
        try:
            with gzip.open(pdb_file_zipped, 'rb') as f_in:
                with open(f"{selected_prot}.pdb", 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            #st.success("Unzipped PDB file successfully.")
        except Exception as e:
            st.error(f"Error unzipping file: {e}")
            st.stop()
            
        command = [
            foldx_executable,
            "--command=BuildModel",
            f"--pdb={selected_prot}.pdb",
            f"--mutant-file={mutant_file}"
        ]

        #st.write("Running FoldX command...")
        try:
            result = subprocess.run(command, capture_output=True, text=True)
            #st.success("FoldX BuildModel completed.")
            #st.code(result.stdout)
            if result.stderr:
                st.warning("Warnings or errors:")
                #st.code(result.stderr)
        except Exception as e:
            st.error(f"Error running FoldX: {e}")
            st.stop()

        # Visualize wild-type and mutant structures
        wt_pdb = f"{selected_prot}.pdb"
        mutant_pdb = f"{selected_prot}_1.pdb"  # FoldX default output

        if os.path.exists(wt_pdb) and os.path.exists(mutant_pdb):
            with open(wt_pdb) as f:
                wt_data = f.read()
            with open(mutant_pdb) as f:
                mut_data = f.read()
            # remove mutant file after loading
            fxout_files = glob.glob("*.fxout")

            if fxout_files:
                for file_path in fxout_files:
                    try:
                        os.remove(file_path)
                    except Exception as e:
                        st.error(f"Error deleting {file_path}: {e}")
            try:
                #remove used files
                os.remove(mutant_pdb)
                os.remove(wt_pdb)
                os.remove(f"WT_{mutant_pdb}")
                os.remove(mutant_file)
                #st.info(f"Removed old files")
            except Exception as e:
                st.warning(f"Could not remove one or more files: {e}")

            view = py3Dmol.view(width=700, height=600, linked=True, viewergrid=(2,2))

            residue_number = int(selected_variant[2:-1])
            chain = selected_variant[1]

            view.addModel(wt_data, viewer=(0,0))
            view.addModel(wt_data, viewer=(1,0))
            view.setStyle({"stick": {"color": "#B0B0B0","scale": 0.4}, "cartoon": {'color': '#B0B0B0'}}, viewer=(0,0))
            view.setStyle({'chain': chain, 'resi': residue_number}, {'stick': {'color': '#017FFD'}}, viewer=(0,0))
            view.setStyle({'chain': chain, 'resi': residue_number}, {'stick': {'color': '#017FFD'}}, viewer=(1,0))
            view.addSurface(py3Dmol.VDW, {'opacity':0.8}, viewer=(1,0))

            view.addModel(mut_data, viewer=(0,1))
            view.addModel(mut_data, viewer=(1,1))
            view.setStyle({"stick": {"color": "#B0B0B0","scale": 0.4}, "cartoon": {'color': '#B0B0B0'}}, viewer=(0,1))
            view.setStyle({'chain': chain, 'resi': residue_number}, {'stick': {'color': '#FF0C57'}}, viewer=(0,1))
            view.setStyle({'chain': chain, 'resi': residue_number}, {'stick': {'color': '#FF0C57'}}, viewer=(1,1))
            view.addSurface(py3Dmol.VDW,{'opacity':0.8}, viewer=(1,1))
            
            view.zoomTo({"chain": chain, "resi": residue_number})
            view.setBackgroundColor("white")
            st.components.v1.html(view._make_html(), height=600)
            st.caption("Py3DMol visualization (wild-type AA in blue, mutant AA in red/pink): Top left = wild-type protein, " \
                "Top right = mutant protein, Bottom left: wild-type protein Van der Waals force surface view, " \
                "Bottom right = mutant protein Van der Waals force surface view." \
                " Controls: Rotate using the left mouseclick, zoom using scroll or right mouseclick.", width="stretch")
        
        else:
            st.error("Could not find output PDB files for visualization.")
        st.session_state.button_disabled = False  # Re-enable button
    else:
        st.markdown("<p style='text-align: center; color: #B0B0B0;'>Note: Large proteins may take longer to visualize </p>", unsafe_allow_html=True)
        