import streamlit as st
from pubmed_utils import fetch_pubmed_articles, extract_pmids_from_text
from mistral_chains import answer_user_query  # updated import
import pandas as pd

st.set_page_config(page_title="Evidence Analysis", layout="wide")
st.title("📊 Evidence Analysis")

# Initialize session state variables
for key in ['articles', 'analysis_complete', 'metadata_df', 'summaries', 
            'comparison', 'conclusion', 'overall_summary']:
    if key not in st.session_state:
        st.session_state[key] = None

# Function to process articles
def process_articles(pmids):
    if len(pmids) > 10:
        st.warning("⚠ Only the first 10 articles will be processed.")
        pmids = pmids[:10]
    
    with st.spinner("Processing articles..."):
        # Fetch articles
        articles = fetch_pubmed_articles(pmids)
        
        if not articles:
            st.error("No articles were retrieved.")
            return None
        
        # Store articles in session state
        st.session_state.articles = articles
        
        # Combine abstracts for multi-article prompts
        combined_abstracts = "".join([art['abstract'] for art in articles])
        
        # Generate overall evidence summary (only overall_summary task)
        if st.session_state.overall_summary is None:
            st.session_state.overall_summary = answer_user_query(
                "Generate a 300-word evidence-based summary synthesizing findings across all articles",
                articles,
                combined_abstracts,
                forced_tasks=["overall_summary"]
            )
        
        # Metadata Table
        if st.session_state.metadata_df is None:
            metadata_data = []
            for art in articles:
                metadata_data.append({
                    "Title": art['title'],
                    "Authors": ", ".join(art['authors']),
                    "Journal": art['journal'],
                    "Year": art['date'],
                    "PMID": art['pmid'],
                    "Abstract": art['abstract'][:200] + "..."
                })
            st.session_state.metadata_df = pd.DataFrame(metadata_data)
        
        # Individual Article Summaries (only per_article_summary task)
        if not st.session_state.summaries:
            summaries = {}
            for i, art in enumerate(articles, 1):
                summaries[i] = answer_user_query(
                f"Generate a 100-word structured summary for this article including Background, Methodology, Key Findings, and Conclusion",
                [art],
                art['abstract'],
                forced_tasks=["per_article_summary"]
            )
        st.session_state.summaries = summaries

    # Store in articles_data for protocol use
    if "articles_data" not in st.session_state:
        st.session_state.articles_data = {}

    st.session_state.articles_data["summaries"] = summaries

        
        # Comparative Analysis (only comparison task)
    if st.session_state.comparison is None:
        st.session_state.comparison = answer_user_query(
             "Create a comparative analysis table of the 3 with these parameters: Treatment, Mechanism of Action, Indication, Onset of Action, Healing Rate, Bleeding Risk, CYP Interaction, GI Side Effects, Trial Design, Population",
             articles,
            combined_abstracts,
            forced_tasks=["comparison"]
            )
        
        # Unified Clinical Conclusion (only conclusion task)
        if st.session_state.conclusion is None:
            st.session_state.conclusion = answer_user_query(
                "Generate an overall evidence-based conclusion summarizing comparative advantages and clinical guidance",
                articles,
                combined_abstracts,
                forced_tasks=["conclusion"]
            )
        
        st.session_state.analysis_complete = True

# Create two columns for input methods
col1, col2 = st.columns(2)

# Column 1: File Upload
with col1:
    st.markdown("### 📁 Upload File")
    uploaded_file = st.file_uploader("Upload PubMed Search Results (TXT file)", type=['txt'])
    
    if uploaded_file is not None:
        content = uploaded_file.getvalue().decode()
        pmids = extract_pmids_from_text(content)
        
        if st.button("Process Uploaded File", key="upload_btn"):
            process_articles(pmids)

# Column 2: Text Input
with col2:
    st.markdown("### 📝 Paste Text")
    st.markdown("""You can paste:
    - Full PubMed search results
    - List of PMIDs
    - PubMed URLs
    - Any text containing PMIDs
    """)
    
    text_input = st.text_area("Paste your text here", height=200)
    
    if text_input:
        pmids = extract_pmids_from_text(text_input)
        
        if pmids:
            st.success(f"✅ Found {len(pmids)} PMIDs in the text.")
            
            if st.button("Process Pasted Text", key="paste_btn"):
                process_articles(pmids)
        else:
            st.warning("❗ No PMIDs found in the pasted text.")

# Display analysis results if available
if st.session_state.analysis_complete and st.session_state.articles:
    st.markdown("---")
    
    # Create tabs for different sections
    analysis_tab1, analysis_tab2, analysis_tab3, analysis_tab4, analysis_tab5 = st.tabs([
        "📝 Overall Summary",
        "📋 Metadata",
        "📚 Article Summaries",
        "🔍 Comparison",
        "🎯 Conclusion"
    ])
    
    with analysis_tab1:
        st.markdown("### Overall Evidence Summary")
        st.markdown(st.session_state.overall_summary)
    
    with analysis_tab2:
        st.markdown("### Metadata Table")
        st.dataframe(st.session_state.metadata_df, use_container_width=True)
        
        st.download_button(
            label="📥 Download Metadata",
            data=st.session_state.metadata_df.to_csv(index=False),
            file_name="metadata.csv",
            mime="text/csv"
        )
    
    with analysis_tab3:
        st.markdown("### Individual Article Summaries")
        for i, summary in st.session_state.summaries.items():
            with st.expander(f"Article {i}"):
                st.markdown(summary)
    
    with analysis_tab4:
        st.markdown("### Comparative Analysis")
        st.markdown(st.session_state.comparison)
    
    with analysis_tab5:
        st.markdown("### Unified Clinical Conclusion")
        st.markdown(st.session_state.conclusion)