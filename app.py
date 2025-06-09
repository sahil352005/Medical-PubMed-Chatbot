import streamlit as st

st.set_page_config(
    page_title="Medical PubMed Analysis",
    page_icon="🩺",
    layout="wide"
)

st.title("🩺 Medical PubMed Analysis Tool")

st.markdown("""
## Welcome to the Medical PubMed Analysis Tool

This tool helps you analyze medical research articles from PubMed in a structured way. Follow these steps:

### Step 1: Search on PubMed
1. Go to [PubMed](https://pubmed.ncbi.nlm.nih.gov/)
2. Perform your search
3. Save the results as a text file (use the 'Save' button on PubMed)

### Step 2: Upload and Process
1. Go to the "Evidence Analysis" page
2. Upload your saved PubMed results
3. Click "Process Articles" to analyze the content

### Step 3: Explore the Results
The tool will generate:
- Overall evidence summary
- Metadata table
- Individual article summaries
- Comparative analysis
- Clinical conclusions

### Step 4: Visualize Data
Visit the "Data Visualization" page to see:
- Healing rates comparison
- Adverse event rates
- Time to symptom relief
- Drug interaction risks

### Step 5: Generate Clinical Documents
The "Clinical Documents" page provides:
- Clinical Trial Protocol Introduction (ICH M11 format)
- Clinical Study Report Discussion (TransCelerate format)

## Getting Started
Click on "Evidence Analysis" in the sidebar to begin!
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