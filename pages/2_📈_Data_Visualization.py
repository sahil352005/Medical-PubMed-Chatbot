import streamlit as st
import plotly.express as px
import pandas as pd
import re

st.set_page_config(page_title="Data Visualization", layout="wide")
st.title("📈 Data Visualization")

# Check if articles are available
if 'articles' not in st.session_state or not st.session_state.articles:
    st.warning("⚠ Please process articles in the Evidence Analysis page first.")
    st.stop()

# Initialize visualization data in session state
if 'viz_data' not in st.session_state:
    st.session_state.viz_data = None

def extract_numerical_data(text: str) -> dict:
    """Extract numerical values using regex patterns."""
    data = {
        'healing_rate': [],
        'adverse_event_rate': [],
        'time_to_relief': [],
        'drug_interaction_risk': []
    }
    patterns = {
        'healing_rate': r'(\d+(?:\.\d+)?)\s*%\s*(?:healing|healed|resolution)',
        'adverse_event_rate': r'(\d+(?:\.\d+)?)\s*%\s*(?:adverse|side effects|complications)',
        'time_to_relief': r'(\d+(?:\.\d+)?)\s*(?:days|weeks|months)',
        'drug_interaction_risk': r'(high|moderate|low)\s*risk'
    }

    text = text.lower()

    for key, pattern in patterns.items():
        for match in re.finditer(pattern, text):
            if key == 'drug_interaction_risk':
                data[key].append(match.group(1).capitalize())
            else:
                try:
                    data[key].append(float(match.group(1)))
                except:
                    continue
    return data

def extract_visualization_data(articles: list) -> pd.DataFrame:
    """Extract and organize data for visualization."""
    if st.session_state.viz_data is not None:
        return st.session_state.viz_data

    all_data = {
        'healing_rate': [],
        'adverse_event_rate': [],
        'time_to_relief': [],
        'drug_interaction_risk': []
    }

    for idx, article in enumerate(articles):
        abstract = article.get('abstract', '')
        if not abstract.strip():
            continue

        data = extract_numerical_data(abstract)
        for key in all_data:
            if key == 'drug_interaction_risk':
                all_data[key].extend([(v, f'Article {idx+1}') for v in data[key]])
            else:
                all_data[key].extend([(v, f'Article {idx+1}') for v in data[key]])

    df_rows = []
    for metric, values in all_data.items():
        if metric == 'drug_interaction_risk':
            risk_df = pd.DataFrame(values, columns=['Category', 'Article'])
            risk_summary = risk_df['Category'].value_counts().reset_index()
            risk_summary.columns = ['Category', 'Value']
            risk_summary['Metric'] = metric
            df_rows.append(risk_summary)
        else:
            df_rows.append(pd.DataFrame({
                'Metric': metric,
                'Value': [v[0] for v in values],
                'Category': [v[1] for v in values]
            }))

    final_df = pd.concat(df_rows, ignore_index=True)
    st.session_state.viz_data = final_df
    return final_df

# Tabs for visualization
tabs = st.tabs(["📊 Healing Rates", "⚠ Adverse Events", "⏱ Time to Relief", "💊 Drug Interactions"])

try:
    df = extract_visualization_data(st.session_state.articles)
    if df.empty:
        st.warning("No numerical data found in the articles for visualization.")
        st.stop()

    tab_titles = ['healing_rate', 'adverse_event_rate', 'time_to_relief', 'drug_interaction_risk']
    tab_labels = ['Healing Rate (%)', 'Adverse Event Rate (%)', 'Time (Days)', 'Interaction Count']

    for i, (tab, metric, y_label) in enumerate(zip(tabs, tab_titles, tab_labels)):
        with tab:
            st.markdown(f"### {metric.replace('_', ' ').title()}")
            metric_data = df[df['Metric'] == metric]
            if metric_data.empty:
                st.info(f"No {metric.replace('_', ' ')} data available.")
            else:
                if metric != 'drug_interaction_risk':
                    fig = px.bar(metric_data, x='Category', y='Value',
                                 labels={'Category': 'Article', 'Value': y_label},
                                 title=f'{y_label} Across Articles')
                else:
                    fig = px.pie(metric_data, values='Value', names='Category',
                                 title='Drug Interaction Risk Distribution')
                st.plotly_chart(fig, use_container_width=True)

    # Download button
    st.download_button(
        label="📥 Download Visualization Data",
        data=df.to_csv(index=False),
        file_name="visualization_data.csv",
        mime="text/csv"
    )

except Exception as e:
    st.error(f"❌ Error while visualizing: {str(e)}")
    st.info("Tip: Ensure articles contain numeric data in expected formats (%, days, risks).")

# Tips section
with st.expander("💡 Tips for Better Visualizations"):
    st.markdown("""
    - Ensure abstracts contain measurable data (percentages, time durations, or risk terms).
    - Standard phrases like "80% healing", "5 days to relief", or "high risk" improve extraction.
    - Try to use similar article formats for consistent graphs.
    """)