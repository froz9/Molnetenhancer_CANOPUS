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
1. Downloads network data from a **GNPS Task ID**.
2. Merges it with your uploaded **SIRIUS/CANOPUS** file.
3. Propagates **NPC** and **ClassyFire** classifications across the network (MolNetEnhancer logic).
""")

# --- Helper Functions ---

def get_gnps_network_data(task_id):
    """Downloads the GNPS output zip and extracts the cluster info file."""
    url = f"https://gnps.ucsd.edu/ProteoSAFe/DownloadResult?task={task_id}&view=download_cytoscape_data"
    
    try:
        response = requests.get(url)
        response.raise_for_status()
    except requests.exceptions.RequestException as e:
        st.error(f"Failed to download from GNPS. Check Task ID. Error: {e}")
        return None

    try:
        with zipfile.ZipFile(io.BytesIO(response.content)) as z:
            # Intelligent search for the correct cluster info file within the zip structure
            # We look for files containing 'clusterinfo' but avoid MAC/OS artifacts
            candidates = [f for f in z.namelist() if "clusterinfo" in f and "__MACOSX" not in f]
            
            target_file = None
            # Priority 1: standard clusterinfo_summary
            for f in candidates:
                if "clusterinfo_summary" in f and not f.endswith('/'):
                    target_file = f
                    break
            
            # Priority 2: If not found, look for other variations (common in older/newer GNPS versions)
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
        st.error("The downloaded file was not a valid Zip. The Task ID might be incorrect or the job failed.")
        return None

def calculate_consensus_score(df, class_col, component_col='componentindex'):
    """
    Calculates the dominant class for each component index.
    Returns a DataFrame with [componentindex, consensus_class, score].
    """
    # Filter for valid annotations only
    valid = df[
        df[class_col].notna() & 
        (df[class_col] != "") & 
        (df[class_col] != "NA")
    ].copy()

    if valid.empty:
        return pd.DataFrame(columns=[component_col, 'consensus_class', 'score'])

    # Count occurrences of each class per component
    counts = valid.groupby([component_col, class_col]).size().reset_index(name='count')
    
    # Sort so the highest count is at the top for each component
    counts = counts.sort_values([component_col, 'count'], ascending=[True, False])
    
    # Calculate total annotations per component (for the score)
    totals = counts.groupby(component_col)['count'].sum().reset_index(name='total')
    
    # Pick the top one (Consensus)
    consensus = counts.drop_duplicates(subset=[component_col], keep='first')
    
    # Merge to calculate score
    result = pd.merge(consensus, totals, on=component_col)
    result['score'] = result['count'] / result['total']
    
    return result[[component_col, class_col, 'score']].rename(columns={class_col: 'consensus_class'})

def process_pipeline(gnps_df, canopus_df):
    # --- 1. Standardize GNPS Columns ---
    # We need 'componentindex' and 'cluster.index'
    # GNPS usually gives 'cluster index', we map it to 'cluster.index' to match R style
    gnps_df.columns = [c.replace('cluster index', 'cluster.index') for c in gnps_df.columns]
    
    if 'componentindex' not in gnps_df.columns or 'cluster.index' not in gnps_df.columns:
        st.error("GNPS file is missing 'componentindex' or 'cluster index' columns.")
        return None
        
    net_data = gnps_df[['componentindex', 'cluster.index']]

    # --- 2. Standardize CANOPUS Columns ---
    # Rename SIRIUS ID to match GNPS
    if 'mappingFeatureId' in canopus_df.columns:
        canopus_df = canopus_df.rename(columns={'mappingFeatureId': 'cluster.index'})

    # Define the mapping from ORIGINAL (Raw) names to Clean names
    # Note: We are fixing the '#' vs '.' issue here by renaming them to clean underscores
    col_mapping = {
        'NPC#pathway': 'NPC_Pathway',
        'NPC#superclass': 'NPC_Superclass',
        'NPC#class': 'NPC_Class',
        'ClassyFire#superclass': 'ClassyFire_Superclass',
        'ClassyFire#class': 'ClassyFire_Class',
        'ClassyFire#subclass': 'ClassyFire_Subclass'
    }
    
    # Rename columns if they exist
    canopus_df = canopus_df.rename(columns=col_mapping)
    
    # Ensure all target columns exist (fill with NA if missing)
    target_cols = list(col_mapping.values())
    for col in target_cols:
        if col not in canopus_df.columns:
            canopus_df[col] = pd.NA

    # Filter to relevant columns
    sirius_subset = canopus_df[['cluster.index'] + target_cols].copy()

    # --- 3. Merge ---
    # Right join: We want all nodes from the GNPS network, even if unannotated
    merged_df = pd.merge(sirius_subset, net_data, on='cluster.index', how='right')

    # --- 4. Propagation Loop ---
    
    # List of (Column to analyze, Output Consensus Name, Output Score Name)
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
        # Calculate consensus
        consensus_df = calculate_consensus_score(merged_df, target_col)
        
        # Merge results back
        if not consensus_df.empty:
            merged_df = pd.merge(merged_df, consensus_df, on='componentindex', how='left')
            merged_df = merged_df.rename(columns={'consensus_class': cons_col_name, 'score': score_col_name})
        else:
            merged_df[cons_col_name] = pd.NA
            merged_df[score_col_name] = pd.NA
        
        # --- Apply Logic: Fill Missing Values ---
        # 1. We only propagate if componentindex != -1 (it's actually a network)
        is_network = merged_df['componentindex'] != -1
        
        # 2. We have a consensus value to use
        has_consensus = merged_df[cons_col_name].notna()
        
        # 3. The original value is missing (NA, Empty, or string "NA")
        original_missing = (merged_df[target_col].isna()) | (merged_df[target_col] == "") | (merged_df[target_col] == "NA")
        
        # Apply fill
        mask = is_network & has_consensus & original_missing
        merged_df.loc[mask, target_col] = merged_df.loc[mask, cons_col_name]
        
        # Update progress
        progress_bar.progress((i + 1) / len(propagation_targets))

    progress_bar.empty()
    return merged_df

# --- Main UI Layout ---

with st.sidebar:
    st.header("Settings")
    
    # 1. GNPS Task ID Input
    task_id_input = st.text_input("GNPS Task ID", value="dacf2cdcdb8b4558a62716f4e9ca4fdb")
    
    # 2. File Upload
    uploaded_file = st.file_uploader("Upload CANOPUS Summary (tsv/txt)", type=['tsv', 'txt', 'csv'])

# --- Execution Logic ---

if uploaded_file and task_id_input:
    st.subheader("Data Processing")
    
    if st.button("Run Analysis"):
        # 1. Load GNPS
        with st.spinner(f"Downloading GNPS Task {task_id_input}..."):
            gnps_data = get_gnps_network_data(task_id_input)
            
        # 2. Load CANOPUS
        if gnps_data is not None:
            try:
                canopus_data = pd.read_csv(uploaded_file, sep='\t')
                st.success("Files loaded successfully.")
                
                # 3. Process
                with st.spinner("Propagating annotations..."):
                    final_df = process_pipeline(gnps_data, canopus_data)
                
                if final_df is not None:
                    st.success("Processing Complete!")
                    
                    # Preview
                    st.write("### Result Preview")
                    st.dataframe(final_df.head(20))
                    
                    # Download
                    csv_buffer = final_df.to_csv(index=False).encode('utf-8')
                    st.download_button(
                        label="Download MolNetEnhancer Results (.csv)",
                        data=csv_buffer,
                        file_name="MolNetEnhancer_Propagated.csv",
                        mime="text/csv"
                    )
            
            except Exception as e:
                st.error(f"Error reading the CANOPUS file: {e}")
else:
    st.info("Please provide a GNPS Task ID and upload your CANOPUS summary file to begin.")