import streamlit as st

st.set_page_config(
    page_title="Medical Research Article Summaries & Regulatory Content Generation",
    page_icon="🩺",
    layout="wide"
)

st.image("https://s3ktech.ai/wp-content/uploads/2025/03/S3Ktech-Logo.png", width=140)

st.title("🩺 Medical Research Article Summaries & Regulatory Content Generation")

st.markdown("""

This tool helps you analyze medical research articles from PubMed in a structured way.

---
### How to Use

#### Option 1: Manual Upload
1. Go to [PubMed](https://pubmed.ncbi.nlm.nih.gov/)
2. Perform your search and save the results as a text file (use the 'Save' button)
3. Go to the **"Article Summaries"** page
4. Upload your saved results and click **Process Articles**
---

#### Option 2: Direct Search (Recommended)
- Go to the **"PubMed Search & Summary"** page in the sidebar.
- Enter your search query.
- Apply filters (year, article type, language, etc.) as needed.
- Click **Search and Summarize** to instantly get:
    - Evidence summaries
    - Metadata tables
    - Comparative analyses
    - Clinical conclusions
---

### Explore the Results

- The tool will generate:
    - Overall evidence summary
    - Metadata table
    - Comparative analysis
    - Clinical conclusions

---

### Generate Regulatory Content

The **"Regulatory Content Generation"** page provides:
- Clinical Trial Protocol Introduction (ICH M11 format)
- Clinical Study Report Discussion (TransCelerate format)

---

## Getting Started

- For the fastest experience, use **"PubMed Search & Summary"** in the sidebar.
- Or use **"Article Summaries"** for manual upload.

""")

# Add some styling
st.markdown("""
<style>
    .main {
        padding: 2rem;
    }
    .stButton>button {
        width: 100%;
    }
</style>
""", unsafe_allow_html=True)
