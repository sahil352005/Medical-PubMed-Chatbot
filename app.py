import streamlit as st
from pubmed_utils import fetch_pubmed_articles, extract_pmids_from_urls, search_pubmed_by_keyword
from vector_store import build_faiss_vectorstore
from mistral_chains import answer_user_query  # Import the centralized function

st.set_page_config(page_title="🧠 Medical PubMed Chatbot", layout="wide")
st.title("🩺 Medical Chatbot for PubMed Articles")

# -------- Step 1: Input Mode --------
input_mode = st.selectbox(
    "Choose Input Mode", 
    ["Enter PMID", "Paste URL", "Search by Keyword", "Paste Raw Text"]
)

user_inputs = []
if input_mode == "Enter PMID":
    user_inputs = st.text_area("Enter up to 10 PMIDs (comma-separated)").split(",")

elif input_mode == "Paste URL":
    urls = st.text_area("Paste up to 10 PubMed URLs (one per line)").splitlines()
    user_inputs = extract_pmids_from_urls(urls)

elif input_mode == "Search by Keyword":
    keyword = st.text_input("Enter keyword for PubMed search")

    if "searched_pmids" not in st.session_state:
        st.session_state.searched_pmids = []

    if keyword and st.button("Search"):
        searched_pmids = search_pubmed_by_keyword(keyword)
        st.session_state.searched_pmids = searched_pmids
        st.success(f"Top {len(searched_pmids)} articles fetched.")

    user_inputs = st.session_state.searched_pmids

elif input_mode == "Paste Raw Text":
    raw_text = st.text_area("Paste full reference text with PMIDs", height=300)
    
    import re
    def extract_pmids_from_text(text):
        pmid_pattern = r"PMID: (\d+)"
        return re.findall(pmid_pattern, text)

    if raw_text:
        user_inputs = extract_pmids_from_text(raw_text)
        if not user_inputs:
            st.warning("❗ No valid PMIDs found in the pasted text.")
        else:
            st.success(f"✅ Found {len(user_inputs)} PMIDs.")



# Clean input and limit to 10
user_inputs = [pmid.strip() for pmid in user_inputs if pmid.strip()]
if len(user_inputs) > 10:
    st.error("⚠️ Please enter only up to 10 articles.")
    st.stop()

# -------- Step 2: Fetch Articles --------
if len(user_inputs) > 0:
    if st.button("Fetch Articles"):
            with st.spinner("🔎 Fetching PubMed articles..."):
                articles = fetch_pubmed_articles(user_inputs)
                if not articles:
                    st.warning("⚠️ No articles were retrieved.")
                    st.stop()

        # Store articles in session state
            st.session_state.articles = articles
            st.success("✅ Articles fetched successfully.")

# -------- Step 3: Query Input --------
if "articles" in st.session_state:
    articles = st.session_state.articles

    st.subheader("📄 Retrieved Articles")
    for i, art in enumerate(articles, 1):
        st.markdown(f"**{i}. {art['title']}**")

    # Prepare article text for processing
    articles_text = "\n\n".join(
        f"Title: {a['title']}\nAuthors: {', '.join(a['authors'])}\nAbstract: {a['abstract']}"
        for a in articles
    )

    st.session_state.articles_text = articles_text

    # Input query
    st.subheader("💬 Ask Your Question")
    user_query = st.text_input("e.g., Compare the studies based on side effects")

    if st.button("Run Query"):
        if not user_query:
            st.warning("Please enter a query.")
        else:
            with st.spinner("🧠 Processing your query..."):
                # Pass all three parameters here as per your mistral_chains.py function signature
                response = answer_user_query(user_query, articles, articles_text)

            st.markdown("### 📝 Response")
            st.markdown(response)
