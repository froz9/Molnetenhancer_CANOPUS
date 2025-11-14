import streamlit as st
import pandas as pd
import requests
import zipfile
import io

# --- Configuration & Title ---
st.set_page_config(page_title="GNPS & SIRIUS Consensus Tool", layout="wide")
st.title("🧬 GNPS & SIRIUS Consensus Propagation Tool")
st.markdown("""
**Workflow:**
1. Paste your **GNPS Task ID** (from the URL of your GNPS job).
2. Upload your **SIRIUS/CANOPUS** summary file.
3. Preview and download the propagated consensus classifications.
""")

# --- Helper Functions ---

def get_gnps_network_data(task_id):
    """Downloads the GNPS output zip and extracts the cluster info file."""
    url = f"https://gnps.ucsd.edu/ProteoSAFe/DownloadResult?task={task_id}&view=download_cytoscape_data"
    
    try:
        # Use POST to satisfy GNPS server
        response = requests.post(url, data={}) 
        response.raise_for_status()
    except requests.exceptions.RequestException as e:
        st.error(f"Failed to download from GNPS. Check Task ID. Error: {e}")
        return None

    try:
        with zipfile.ZipFile(io.BytesIO(response.content)) as z:
            candidates = [f for f in z.namelist() if "clusterinfo" in f and "__MACOSX" not in f]
            target_file = None
            
            for f in candidates:
                if "clusterinfo_summary" in f and not f.endswith('/'):
                    target_file = f
                    break
            
            if not target_file and candidates:
                target_file = candidates[0]

            if target_file:
                with z.open(target_file) as f:
                    df = pd.read_csv(f, sep='\t')
                return df
            else:
                st.error("Could not locate a 'clusterinfo' file inside the GNPS Task Zip.")
                return None
    except zipfile.BadZipFile:
        st.error("The downloaded file was not a valid Zip. The Task ID might be incorrect.")
        return None

def calculate_consensus_score(df, class_col, component_col='componentindex'):
    """Calculates the dominant class for each component index."""
    valid = df[
        df[class_col].notna() & 
        (df[class_col] != "") & 
        (df[class_col] != "NA")
    ].copy()

    if valid.empty:
        return pd.DataFrame(columns=[component_col, 'consensus_class', 'score'])

    counts = valid.groupby([component_col, class_col]).size().reset_index(name='count')
    counts = counts.sort_values([component_col, 'count'], ascending=[True, False])
    totals = counts.groupby(component_col)['count'].sum().reset_index(name='total')
    consensus = counts.drop_duplicates(subset=[component_col], keep='first')
    
    result = pd.merge(consensus, totals, on=component_col)
    result['score'] = result['count'] / result['total']
    
    return result[[component_col, class_col, 'score']].rename(columns={class_col: 'consensus_class'})

def process_pipeline(gnps_df, canopus_df):
    # --- Standardize GNPS ---
    gnps_df.columns = [c.replace('cluster index', 'cluster.index') for c in gnps_df.columns]
    
    if 'componentindex' not in gnps_df.columns or 'cluster.index' not in gnps_df.columns:
        st.error("GNPS file is missing 'componentindex' or 'cluster index' columns.")
        return None
        
    net_data = gnps_df[['componentindex', 'cluster.index']]

    # --- Standardize CANOPUS ---
    if 'mappingFeatureId' in canopus_df.columns:
        canopus_df = canopus_df.rename(columns={'mappingFeatureId': 'cluster.index'})

    col_mapping = {
        'NPC#pathway': 'NPC_Pathway',
        'NPC#superclass': 'NPC_Superclass',
        'NPC#class': 'NPC_Class',
        'ClassyFire#superclass': 'ClassyFire_Superclass',
        'ClassyFire#class': 'ClassyFire_Class',
        'ClassyFire#subclass': 'ClassyFire_Subclass'
    }
    
    canopus_df = canopus_df.rename(columns=col_mapping)
    
    target_cols = list(col_mapping.values())
    for col in target_cols:
        if col not in canopus_df.columns:
            canopus_df[col] = pd.NA

    sirius_subset = canopus_df[['cluster.index'] + target_cols].copy()

    # --- Merge & Propagate ---
    merged_df = pd.merge(sirius_subset, net_data, on='cluster.index', how='right')

    propagation_targets = [
        ('NPC_Pathway', 'NPC_Pathway_Consensus', 'NPC_Pathway_Score'),
        ('NPC_Superclass', 'NPC_Superclass_Consensus', 'NPC_Superclass_Score'),
        ('NPC_Class', 'NPC_Class_Consensus', 'NPC_Class_Score'),
        ('ClassyFire_Superclass', 'ClassyFire_Superclass_Consensus', 'ClassyFire_Superclass_Score'),
        ('ClassyFire_Class', 'ClassyFire_Class_Consensus', 'ClassyFire_Class_Score'),
        ('ClassyFire_Subclass', 'ClassyFire_Subclass_Consensus', 'ClassyFire_Subclass_Score'),
    ]

    progress_bar = st.progress(0)
    
    for i, (target_col, cons_col_name, score_col_name) in enumerate(propagation_targets):
        consensus_df = calculate_consensus_score(merged_df, target_col)
        
        if not consensus_df.empty:
            merged_df = pd.merge(merged_df, consensus_df, on='componentindex', how='left')
            merged_df = merged_df.rename(columns={'consensus_class': cons_col_name, 'score': score_col_name})
        else:
            merged_df[cons_col_name] = pd.NA
            merged_df[score_col_name] = pd.NA
        
        # Propagate Logic
        is_network = merged_df['componentindex'] != -1
        has_consensus = merged_df[cons_col_name].notna()
        original_missing = (merged_df[target_col].isna()) | (merged_df[target_col] == "") | (merged_df[target_col] == "NA")
        
        mask = is_network & has_consensus & original_missing
        merged_df.loc[mask, target_col] = merged_df.loc[mask, cons_col_name]
        
        progress_bar.progress((i + 1) / len(propagation_targets))

    progress_bar.empty()
    return merged_df

# --- Main UI Layout ---

with st.sidebar:
    st.header("Settings")
    
    task_id_input = st.text_input(
        "GNPS Task ID", 
        placeholder="e.g., dacf2cdcdb8b4558a62716f4e9ca4fdb"
    )
    
    uploaded_file = st.file_uploader("Upload CANOPUS Summary (tsv/txt)", type=['tsv', 'txt', 'csv'])

    st.markdown("---")
    run_btn = st.button("Run Analysis", type="primary")

# --- Execution Logic ---

# 1. IF button is clicked, run logic and save to session state
if run_btn:
    if uploaded_file and task_id_input:
        st.subheader("Data Processing")
    
    if st.button("Run Analysis"):
        with st.spinner(f"Downloading GNPS Task {task_id_input}..."):
            gnps_data = get_gnps_network_data(task_id_input)
            
        if gnps_data is not None:
            try:
                # Reset file pointer just in case
                uploaded_file.seek(0)
                canopus_data = pd.read_csv(uploaded_file, sep='\t')
                
                with st.spinner("Propagating annotations..."):
                    final_df = process_pipeline(gnps_data, canopus_data)
                
                if final_df is not None:
                    # SAVE TO SESSION STATE
                    st.session_state['processed_data'] = final_df
                    st.success("Processing Complete!")
            except Exception as e:
                st.error(f"Error reading the CANOPUS file: {e}")
    else:
        st.error("Please provide both a Task ID and a CANOPUS file.")

# 2. CHECK if data exists in session state (This persists across re-runs/downloads)
if 'processed_data' in st.session_state:
    final_df = st.session_state['processed_data']

    # --- Prepare Output Dataframes ---
    base_cols = ['cluster.index']
    
    cols_npc = [c for c in final_df.columns if 'NPC' in c and ('Consensus' in c)]
    cols_classy = [c for c in final_df.columns if 'ClassyFire' in c and ('Consensus' in c)]
    
    df_npc = final_df[base_cols + cols_npc]
    df_classy = final_df[base_cols + cols_classy]
    df_all = final_df[base_cols + cols_npc + cols_classy]
    
    # --- PREVIEW ---
    st.subheader("🔍 Result Preview (Combined)")
    st.dataframe(df_all.head(50))
    
    # --- DOWNLOADS ---
    st.subheader("📥 Download Results")
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.download_button(
            label="Download NPC Only",
            data=df_npc.to_csv(index=False).encode('utf-8'),
            file_name="Consensus_NPC.csv",
            mime="text/csv"
        )
    
    with col2:
        st.download_button(
            label="Download ClassyFire Only",
            data=df_classy.to_csv(index=False).encode('utf-8'),
            file_name="Consensus_ClassyFire.csv",
            mime="text/csv"
        )
        
    with col3:
        st.download_button(
            label="Download Combined",
            data=df_all.to_csv(index=False).encode('utf-8'),
            file_name="Consensus_All.csv",
            mime="text/csv"
        )