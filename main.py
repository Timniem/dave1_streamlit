import os
import subprocess
import gzip
import shutil
import glob

import streamlit as st
import pandas as pd
import py3Dmol


# Title
st.title("VKGL Datasharing VUS Analysis")

# Load VUS data
vus_path = "vkgl_apr2024_VUS_pred.csv"
if not os.path.exists(vus_path):
    st.error("CSV file not found. Please check the path.")
    st.stop()


vkgl_consensus_vus = pd.read_csv(vus_path)


edited_df = st.dataframe(
    vkgl_consensus_vus[["gene","TranscriptID","dna_variant_chrom","dna_variant_pos","dna_variant_ref","dna_variant_alt","delta_aaSeq","LP"]],
    width=800,
    hide_index=True,
    on_select="rerun",
    selection_mode="single-row",
)

selected_rows = edited_df.selection.rows

if len(selected_rows) > 0:
    selected_df = vkgl_consensus_vus.iloc[selected_rows[0]]
    
    selected_gene = selected_df['UniProtID']
    selected_variant = selected_df['delta_aaSeq']
    #st.write(selected_gene)

    # Input for FoldX and PDB
    foldx_executable =  "./foldx/foldx5_1Linux64/foldx_20251231"
    pdb_file_zipped = f'./AF_pdb/AF-{selected_gene}-F1-model_v4.pdb.gz'

    if st.button("Visualize Variant"):
        # Generate mutant file content
        mutant_string = f"{selected_variant};"
        mutant_file = "individual_list.txt"
        with open(mutant_file, "w") as f:
            f.write(mutant_string)

        try:
            with gzip.open(pdb_file_zipped, 'rb') as f_in:
                with open(f"{selected_gene}.pdb", 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            #st.success("Unzipped PDB file successfully.")
        except Exception as e:
            st.error(f"Error unzipping file: {e}")
            st.stop()
            
        command = [
            foldx_executable,
            "--command=BuildModel",
            f"--pdb={selected_gene}.pdb",
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
        wt_pdb = f"{selected_gene}.pdb"
        mutant_pdb = f"{selected_gene}_1.pdb"  # FoldX default output

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

            view.addModel(wt_data, viewer=(0,0))
            view.addModel(wt_data, viewer=(1,0))
            view.setStyle({"stick": {"color": "grey","scale": 0.4}, "cartoon": {'color': 'grey'}}, viewer=(0,0))
            view.addSurface(py3Dmol.VDW, {'opacity':0.7}, viewer=(1,0))

            view.addModel(mut_data, viewer=(0,1))
            view.addModel(mut_data, viewer=(1,1))
            view.setStyle({"stick": {"color": "grey","scale": 0.4}, "cartoon": {'color': 'grey'}}, viewer=(0,1))
            view.addSurface(py3Dmol.VDW,{'opacity':0.7}, viewer=(1,1))

            residue_number = int(selected_variant[2:-1])
            chain = selected_variant[1]
            view.setStyle({'chain': chain, 'resi': residue_number}, {'stick': {'color': 'yellow'}})
            view.zoomTo({"chain": chain, "resi": residue_number})
            view.setBackgroundColor("white")
            
            st.components.v1.html(view._make_html(), height=600)
        else:
            st.error("Could not find output PDB files for visualization.")
