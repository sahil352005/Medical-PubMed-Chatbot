import streamlit as st
import pandas as pd
import plotly.express as px
from mistral_chains import run_comparison

st.set_page_config(page_title="📊 Medical Article Comparison", layout="wide")
st.title("📊 Comparative Visualization of PubMed Articles")

# --- Step 1: Ensure articles are loaded ---
if "articles_text" not in st.session_state:
    st.warning("⚠️ No articles found. Please go to the main page and input articles.")
    st.stop()

articles_text = st.session_state.articles_text

# --- Step 2: Get Comparison Table from Mistral ---
if "comparison_table" not in st.session_state:
    with st.spinner("🔍 Generating comparison table from Mistral..."):
        st.session_state.comparison_table = run_comparison(articles_text)

markdown_table = st.session_state.comparison_table

# --- Step 3: Parse Markdown Table into DataFrame ---
def parse_markdown_table(md_table):
    lines = [line.strip() for line in md_table.split('\n') if line.strip()]
    if len(lines) < 3:
        return pd.DataFrame()
    
    headers = [h.strip() for h in lines[0].strip('|').split('|')]
    rows = []
    for row in lines[2:]:
        values = [v.strip() for v in row.strip('|').split('|')]
        if len(values) == len(headers):
            rows.append(values)

    return pd.DataFrame(rows, columns=headers)

df = parse_markdown_table(markdown_table)

if df.empty:
    st.error("❌ Could not parse the comparison table into data. Try reloading.")
    st.stop()

# --- Step 4: Show Raw Table ---
st.subheader("🧾 Article Comparison Table")
st.dataframe(df, use_container_width=True)

# --- Step 5: Plot Interactive Charts using Plotly ---
def plot_metric_bar_plotly(metric: str, label: str):
    if metric not in df.columns:
        st.warning(f"⚠️ '{metric}' column is missing.")
        return

    try:
        chart_data = df[[df.columns[0], metric]].copy()
        chart_data[metric] = chart_data[metric].str.extract(r"([\d.]+)").astype(float)
        fig = px.bar(
            chart_data,
            x=chart_data.columns[0],
            y=metric,
            title=label,
            labels={chart_data.columns[0]: "Article", metric: label},
            color=chart_data.columns[0],
        )
        st.plotly_chart(fig, use_container_width=True)
    except Exception as e:
        st.warning(f"⚠️ Could not visualize '{metric}': {str(e)}")

# Plot key medical metrics
plot_metric_bar_plotly("Healing Rate", "Healing Rates")
plot_metric_bar_plotly("Bleeding Risk", "Adverse Event Incidence")
plot_metric_bar_plotly("CYP Interactions", "Drug-Drug Interaction Potential")
plot_metric_bar_plotly("Onset of Action", "Time to Symptom Relief")
