import os
import base64
import tarfile
from pathlib import Path


import streamlit as st
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
        img_html = f"<img src='data:image/png;base64,{img_to_bytes(img_path)}' align={align} style='width:{size}px' class='img-fluid'>"
    else:
        img_html = f"<img src='data:image/png;base64,{img_to_bytes(img_path)}' style='width:{size}px' class='img-fluid'>"
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

HEADER_TEXT = "Pathogenicity prediction by the Digital Approximation of Variant Effect, version 1 (DAVE1) model on 11k+ VUS variants from the April 2024 release of the <a href=https://vkgl.molgeniscloud.org target='_blank' rel='noopener noreferrer' > VKGL datasharing</a>. \
                The ideal cutoff point for Likely Pathogenic classification >= 0.286. For more information please refer to: DOI"
###


### init ###
if 'button_disabled' not in st.session_state:
    st.session_state.button_disabled = False
###


features_shap = [f"{key}.sph" for key in FEATURE_NAMES_DICT.keys()]

# Top Bar

st.markdown(f"""
<div style='background-color: #017FFD; padding:10px; display: flex; justify-content: space-between; align-items: center;'>
    <a href='https://molgenis.org' style='display: inline-block;' target='_blank' rel='noopener noreferrer'>
      {img_to_html('images/logo_blue.png', size=150)}
    </a>
    {img_to_html('images/umcg_logo.png', align='right', size=100)}
</div>
""", unsafe_allow_html=True)

# Title
st.markdown("<p style='text-align: center; color: #666; font-size:32px;'>DAVE1 scores VKGL Datasharing VUS</p>", unsafe_allow_html=True)
st.markdown(f"<p style='text-align: left; font-size: small; color: #B0B0B0;'>{HEADER_TEXT}</p>", unsafe_allow_html=True)
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
    st.markdown("---")
    selected_df = filtered_data.iloc[selected_rows[0]]

    if type(selected_df['seqFt']) != float:
        feature = selected_df['seqFt'].split("|")[-1].split("~")
    else:
        feature = None

    selected_df = selected_df.rename(FEATURE_NAMES_DICT)

    st.markdown(f"<p style='text-align: center; color: #333;'>Feature contribution: <i>{selected_df['gene']}</i> {selected_df['AA change']} </p>", unsafe_allow_html=True)

    fig = explain_plot_plotly(selected_df[features_shap], selected_df[FEATURE_NAMES_DICT.values()], FEATURE_NAMES_DICT.values(), selected_df["LP"])
    st.plotly_chart(fig, use_container_width=True, config={"displayModeBar": False, "scrollZoom": False,"doubleClick": False })

    st.markdown(f"<p style='text-align: left; color: #555;'>Residue information</p>", unsafe_allow_html=True)
    st.markdown(f"<p style='text-align: left; font-size: small; color: #999;'>Localization: {selected_df['ann_proteinLocalization']}</p>", unsafe_allow_html=True)
    if feature:
        st.markdown(f"<p style='text-align: left; font-size: small; color: #999;'>Sequence features: {' -> '.join(feature)}</p>", unsafe_allow_html=True)
    
    st.markdown("---")

    c1, c2, c3 = st.columns(3)

    
    with c1:
        left, center, right = st.columns([1, 1, 1])  
        with center:
            toggle_state = st.toggle(" ", disabled=st.session_state.button_disabled)
        st.markdown("<p style='text-align: left; color: #B0B0B0; font-size: 10px;'>Show Solvent Accessible Surface (SAS)<br>Default = Van der Waals forces (VDW)</p>", unsafe_allow_html=True)
        
        
    with c2:
        viz_button = st.button("Visualize Variant", width="stretch", on_click=disable_button, disabled=st.session_state.button_disabled )

    if viz_button:

        try:
            wt_pdb = f"protein_structures/{selected_df['UniProtID']}.pdb"
            mutant_pdb = f"protein_structures/{selected_df['UniProtID']}_{selected_df['delta_aaSeq']}.pdb"
            # Read PDB contents from archive
            with tarfile.open("./mut_wt_structures_vkgl_vus.tar.gz", "r:gz") as tar:
                wt_file = tar.extractfile(wt_pdb)
                mutant_file = tar.extractfile(mutant_pdb)
                wt_data = wt_file.read().decode() if wt_file else ""
                mut_data = mutant_file.read().decode() if mutant_file else ""
        except Exception as e:
            st.error(f"Error opening files: {e}")
            st.session_state.button_disabled = False  # Re-enable button
            st.stop()
        
        view = py3Dmol.view(width=700, height=600, linked=True, viewergrid=(2,2))

        residue_number = int(selected_df['delta_aaSeq'][2:-1])
        chain = selected_df['delta_aaSeq'][1]
        # Add the molecules
        view.addModel(wt_data, viewer=(0,0))
        view.addModel(wt_data, viewer=(1,0))
        view.addModel(mut_data, viewer=(0,1))
        view.addModel(mut_data, viewer=(1,1))
        view.setStyle({})  # Hide all atom/bond styles

        # WT
        view.setStyle({"stick": {"color": "#B0B0B0","scale": 0.4}, "cartoon": {'color': '#B0B0B0'}}, viewer=(0,0))
        view.setStyle({'chain': chain, 'resi': residue_number}, {'stick': {'color': '#017FFD'}}, viewer=(0,0))
        view.setStyle({'chain': chain, 'resi': residue_number}, {'stick': {'color': '#017FFD'}}, viewer=(1,0))
        if toggle_state:
            view.addSurface(py3Dmol.SAS, {'opacity':0.8}, viewer=(1,0))
            surface_shown = "solvent accessible surface"
        else:
            view.addSurface(py3Dmol.VDW, {'opacity':0.8}, viewer=(1,0))
            surface_shown = "Van der Waals force"
        
        # MT
        view.setStyle({"stick": {"color": "#B0B0B0","scale": 0.4}, "cartoon": {'color': '#B0B0B0'}}, viewer=(0,1))
        view.setStyle({'chain': chain, 'resi': residue_number}, {'stick': {'color': '#FF0C57'}}, viewer=(0,1))
        view.setStyle({'chain': chain, 'resi': residue_number}, {'stick': {'color': '#FF0C57'}}, viewer=(1,1))

        if toggle_state:
            view.addSurface(py3Dmol.SAS, {'opacity':0.8}, viewer=(1,1))
            
        else:
            view.addSurface(py3Dmol.VDW, {'opacity':0.8}, viewer=(1,1))
            
        
        view.zoomTo({"chain": chain, "resi": residue_number})
        view.setBackgroundColor("white")
        st.components.v1.html(view._make_html(), height=600)
        st.markdown("<p style='text-align: center; color: #B0B0B0; font-size:smaller;' > <a href=https://github.com/avirshup/py3dmol >Py3DMol</a> visualization (wild-type AA in blue, mutant AA in red/pink): Top left = wild-type protein, " \
            f"Top right = mutant protein, Bottom left: wild-type protein {surface_shown} surface view, " \
            f"Bottom right = mutant protein {surface_shown} surface view." \
            " Controls: Rotate using the left mouseclick, zoom using scroll or right mouseclick.</p>", width="stretch", unsafe_allow_html=True)
            
        st.session_state.button_disabled = False  # Re-enable button
    else:
        st.markdown("<p style='text-align: center; color: #B0B0B0;'>Note: Large proteins may take longer to visualize </p>", unsafe_allow_html=True)
        