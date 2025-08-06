import streamlit as st
from pubmed_utils import search_pubmed_by_keyword, fetch_pubmed_articles
from mistral_chains import answer_user_query
import pandas as pd
from datetime import date

st.set_page_config(page_title="PubMed Search & Summary", layout="wide")

st.markdown("""
    <style>
    .big-title {
        font-size:2.2rem;
        font-weight:700;
        text-align:center;
        margin-bottom:0.5em;
    }
    .subtitle {
        font-size:1.2rem;
        font-weight:500;
        color:#343a40;
        margin-bottom:0.5em;
    }
    .card {
        background: #f3f4f6;
        border-radius: 12px;
        padding: 1.5em 1.5em 1em 1.5em;
        margin-bottom: 1.5em;
        box-shadow: 0px 2px 8px rgba(0,0,0,0.04);
        border: 1px solid #dee2e6;
    }
    .sample-q {
        background: #e9ecef;
        border-radius: 8px;
        padding: 0.5em 1em;
        margin: 0.2em 0.2em 0.2em 0;
        display: inline-block;
        cursor: pointer;
        color: #343a40;
        font-weight: 500;
        font-size: 1rem;
        border: 1px solid #dee2e6;
    }
    .stButton>button {
        background-color: #343a40;
        color: white;
        border-radius: 8px;
        padding: 0.6em 1.2em;
        font-size: 1rem;
        font-weight: 600;
        border: none;
        box-shadow: 0px 3px 6px rgba(0,0,0,0.08);
    }
    </style>
""", unsafe_allow_html=True)

st.image("https://s3ktech.ai/wp-content/uploads/2025/03/S3Ktech-Logo.png", width=140)
st.markdown('<div class="big-title">🔎 PubMed Search & Summary</div>', unsafe_allow_html=True)

st.markdown("""
<div class="card">
<b>How to use this page:</b>
<ul style="margin-bottom:0;">
  
  <li><span style="color:#343a40"><b>Search PubMed directly below</b></span> and get instant summaries and tables.</li>
</ul>
</div>
""", unsafe_allow_html=True)

with st.expander("ℹ️ About this tool", expanded=False):
    st.write("""
    This tool lets you search PubMed, fetch up to 10 articles, and instantly generate evidence summaries, metadata tables, and comparative analyses using AI.
    """)

st.markdown('<div class="subtitle">Sample Questions</div>', unsafe_allow_html=True)
sample_questions = [
    "What are the latest treatments for type 2 diabetes?",
    "Compare proton pump inhibitors and H2 blockers for peptic ulcer.",
    "Recent advances in immunotherapy for lung cancer.",
    "What is the efficacy of GLP-1 agonists in obesity?",
    "Are SGLT2 inhibitors safe in heart failure patients?"
]
qcols = st.columns(len(sample_questions))
for i, q in enumerate(sample_questions):
    if qcols[i].button(q, key=f"sampleq_{i}"):
        st.session_state['sample_query'] = q

st.markdown('<div class="subtitle">PubMed Search</div>', unsafe_allow_html=True)

with st.form("search_form"):
    col1, col2 = st.columns([3,1])
    with col1:
        query = st.text_input(
            "Enter your search query",
            value=st.session_state.get('sample_query', ""),
            placeholder="e.g., proton pump inhibitors peptic ulcer"
        )
    with col2:
        max_results = st.slider("Articles", 1, 10, 5, key="max_results_slider")
    st.markdown("#### Filters (optional)")
    colf1, colf2, colf3, colf4, colf5 = st.columns(5)
    with colf1:
        year_range = st.date_input(
            "Publication Year Range",
            value=(date(2020,1,1), date.today()),
            min_value=date(1950,1,1),
            max_value=date.today(),
            key="year_range"
        )
        year_from = str(year_range[0].year) if isinstance(year_range, tuple) else ""
        year_to = str(year_range[1].year) if isinstance(year_range, tuple) else ""
    with colf2:
            article_type = st.selectbox("Article Type", ["", "Clinical Trial", "Review", "Meta-Analysis", "Randomized Controlled Trial"])
    with colf3:
            language = st.selectbox("Language", ["", "English", "Spanish", "French", "German", "Chinese"])
            status = None
    with colf4:
        sex = st.selectbox("Sex", ["", "Male", "Female", "Both"])
    with colf5:
        species = st.selectbox("Species", ["", "Human", "Animal", "Both"])
    submitted = st.form_submit_button("🔍 Search and Summarize")

def build_pubmed_query(base_query, year_from, year_to, article_type, language, sex, species):
    q = base_query
    if year_from and year_to:
        q += f" AND ({year_from}:{year_to}[dp])"
    elif year_from:
        q += f" AND ({year_from}[dp])"
    elif year_to:
        q += f" AND ({year_to}[dp])"
    if article_type:
        pt_map = {
            "Clinical Trial": "clinical trial[pt]",
            "Review": "review[pt]",
            "Meta-Analysis": "meta-analysis[pt]",
            "Randomized Controlled Trial": "randomized controlled trial[pt]"
        }
        q += f" AND {pt_map.get(article_type, '')}"
    if language:
        lang_map = {
            "English": "english[lang]",
            "Spanish": "spanish[lang]",
            "French": "french[lang]",
            "German": "german[lang]",
            "Chinese": "chinese[lang]"
        }
        q += f" AND {lang_map.get(language, '')}"
    if sex:
        sex_map = {
            "Male": "male[sex]",
            "Female": "female[sex]",
            "Both": "male[sex] OR female[sex]"
        }
        if sex != "Both":
            q += f" AND {sex_map.get(sex, '')}"
        else:
            q += f" AND (male[sex] OR female[sex])"
    if species:
        species_map = {
            "Human": "humans[mh]",
            "Animal": "animals[mh]",
            "Both": "humans[mh] OR animals[mh]"
        }
        if species != "Both":
            q += f" AND {species_map.get(species, '')}"
        else:
            q += f" AND (humans[mh] OR animals[mh])"
    return q

if submitted:
    if not query.strip():
        st.warning("Please enter a search query.")
    else:
        full_query = build_pubmed_query(
            query,
            year_from,
            year_to,
            article_type if article_type else "",
            language if language else "",
            sex if sex else "",
            species if species else ""
        )
        with st.spinner("🔄 Generating summaries and fetching articles..."):
            pmids = search_pubmed_by_keyword(full_query, max_results=max_results)
            if not pmids:
                st.error("No articles found for this query.")
            else:
                articles = fetch_pubmed_articles(pmids)
                if not articles:
                    st.error("No article details could be retrieved.")
                else:
                    combined_abstracts = "".join([art['abstract'] for art in articles])
                    st.session_state.articles = articles
                    st.session_state.combined_abstracts = combined_abstracts
                    st.session_state.overall_summary = answer_user_query(
                        "Generate a 300-word evidence-based summary synthesizing findings across all articles",
                        articles,
                        combined_abstracts,
                        forced_tasks=["overall_summary"]
                    )
        if (
            "articles" in st.session_state and st.session_state.articles
            and "combined_abstracts" in st.session_state and st.session_state.combined_abstracts
        ):
            st.success("Summaries generated successfully!")

# Always show summaries/tabs/chatbot if articles are in session_state
if (
    "articles" in st.session_state and st.session_state.articles
    and "combined_abstracts" in st.session_state and st.session_state.combined_abstracts
):
    articles = st.session_state.articles
    combined_abstracts = st.session_state.combined_abstracts
    overall_summary = st.session_state.get("overall_summary", "")
    st.success(f"Found {len(articles)} articles for your query.")
    tab1, tab2, tab3, tab4, tab5 = st.tabs([
        "📝 Overall Summary",
        "📋 Metadata Table",
        "📚 Article Summaries",
        "🔍 Comparison Table",
        "🎯 Unified Conclusion"
    ])
    with tab1:
        st.markdown("#### Overall Evidence Summary")
        st.markdown(overall_summary)
    with tab2:
        st.markdown("#### Article Metadata Table")
        metadata_data = []
        for art in articles:
            metadata_data.append({
                "Title": art['title'],
                "Authors": ", ".join(art['authors']),
                "Journal": art['journal'],
                "Year": art['date'],
                "PMID": art['pmid'],
                "URL": art.get('url', "")
            })
        df = pd.DataFrame(metadata_data)
        st.dataframe(df, use_container_width=True)
        st.download_button(
            label="📥 Download Metadata as CSV",
            data=df.to_csv(index=False),
            file_name="pubmed_metadata.csv",
            mime="text/csv"
        )
    with tab3:
        st.markdown("#### Individual Article Summaries")
        for i, art in enumerate(articles, 1):
            summary = answer_user_query(
                "Generate a 100-word structured summary for this article including Background, Methodology, Key Findings, and Conclusion",
                [art],
                art['abstract'],
                forced_tasks=["per_article_summary"]
            )
            with st.expander(f"Article {i}: {art['title']}"):
                st.markdown(summary)
    with tab4:
        st.markdown("#### Comparative Analysis Table")
        st.markdown(answer_user_query(
            "Create a comparative analysis table of the articles with these parameters: Treatment, Mechanism of Action, Indication, Onset of Action, Healing Rate, Bleeding Risk, CYP Interaction, GI Side Effects, Trial Design, Population",
            articles,
            combined_abstracts,
            forced_tasks=["comparison"]
        ))
    with tab5:
        st.markdown("#### Unified Clinical Conclusion")
        st.markdown(answer_user_query(
            "Generate an overall evidence-based conclusion summarizing comparative advantages and clinical guidance",
            articles,
            combined_abstracts,
            forced_tasks=["conclusion"]
        ))