import streamlit as st
import plotly.express as px
import pandas as pd
import re
import time
import json

st.set_page_config(page_title="Data Visualization", layout="wide")
st.title("📈 Data Visualization")

# Check if articles are available
if not st.session_state.get('articles'):
    st.warning("⚠️ Please process articles in the Evidence Analysis page first.")
    st.stop()

# Initialize visualization data in session state
if 'viz_data' not in st.session_state:
    st.session_state.viz_data = None

def extract_numerical_data(text):
    """Extract numerical data from text using regex patterns."""
    data = {
        'healing_rate': [],
        'adverse_event_rate': [],
        'time_to_relief': [],
        'drug_interaction_risk': []
    }
    
    # Patterns for different types of data
    patterns = {
        'healing_rate': r'(\d+(?:\.\d+)?)\s*%\s*(?:healing|healed|resolution)',
        'adverse_event_rate': r'(\d+(?:\.\d+)?)\s*%\s*(?:adverse|side effects|complications)',
        'time_to_relief': r'(\d+(?:\.\d+)?)\s*(?:days|weeks|months)',
        'drug_interaction_risk': r'(?:high|moderate|low)\s*risk'
    }
    
    # Extract data using patterns
    for key, pattern in patterns.items():
        matches = re.finditer(pattern, text.lower())
        for match in matches:
            if key == 'drug_interaction_risk':
                data[key].append(match.group(0))
            else:
                try:
                    value = float(match.group(1))
                    data[key].append(value)
                except ValueError:
                    continue
    
    return data

def extract_visualization_data(articles):
    """Extract data for visualization from articles."""
    if st.session_state.viz_data is not None:
        return st.session_state.viz_data
    
    all_data = {
        'healing_rate': [],
        'adverse_event_rate': [],
        'time_to_relief': [],
        'drug_interaction_risk': []
    }
    
    for article in articles:
        text = article['abstract']
        data = extract_numerical_data(text)
        
        for key in all_data:
            all_data[key].extend(data[key])
    
    # Convert to DataFrame
    df = pd.DataFrame({
        'Metric': [],
        'Value': [],
        'Category': []
    })
    
    for metric, values in all_data.items():
        if values:
            if metric == 'drug_interaction_risk':
                # Count occurrences of each risk level
                risk_counts = pd.Series(values).value_counts()
                for risk, count in risk_counts.items():
                    df = pd.concat([df, pd.DataFrame({
                        'Metric': [metric],
                        'Value': [count],
                        'Category': [risk]
                    })])
            else:
                # Add numerical values
                df = pd.concat([df, pd.DataFrame({
                    'Metric': [metric] * len(values),
                    'Value': values,
                    'Category': [f'Article {i+1}' for i in range(len(values))]
                })])
    
    st.session_state.viz_data = df
    return df

# Create tabs for different visualizations
viz_tab1, viz_tab2, viz_tab3, viz_tab4 = st.tabs([
    "📊 Healing Rates",
    "⚠️ Adverse Events",
    "⏱️ Time to Relief",
    "💊 Drug Interactions"
])

try:
    df = extract_visualization_data(st.session_state.articles)
    
    if df.empty:
        st.warning("No numerical data found in the articles for visualization.")
        st.stop()
    
    with viz_tab1:
        st.markdown("### Healing Rates")
        healing_data = df[df['Metric'] == 'healing_rate']
        if not healing_data.empty:
            fig = px.bar(healing_data, x='Category', y='Value',
                        title='Healing Rates Across Articles',
                        labels={'Value': 'Healing Rate (%)', 'Category': 'Article'})
            st.plotly_chart(fig, use_container_width=True)
        else:
            st.info("No healing rate data available.")
    
    with viz_tab2:
        st.markdown("### Adverse Event Rates")
        adverse_data = df[df['Metric'] == 'adverse_event_rate']
        if not adverse_data.empty:
            fig = px.bar(adverse_data, x='Category', y='Value',
                        title='Adverse Event Rates Across Articles',
                        labels={'Value': 'Adverse Event Rate (%)', 'Category': 'Article'})
            st.plotly_chart(fig, use_container_width=True)
        else:
            st.info("No adverse event data available.")
    
    with viz_tab3:
        st.markdown("### Time to Symptom Relief")
        time_data = df[df['Metric'] == 'time_to_relief']
        if not time_data.empty:
            fig = px.bar(time_data, x='Category', y='Value',
                        title='Time to Symptom Relief',
                        labels={'Value': 'Time (days)', 'Category': 'Article'})
            st.plotly_chart(fig, use_container_width=True)
        else:
            st.info("No time to relief data available.")
    
    with viz_tab4:
        st.markdown("### Drug Interaction Risk Distribution")
        interaction_data = df[df['Metric'] == 'drug_interaction_risk']
        if not interaction_data.empty:
            fig = px.pie(interaction_data, values='Value', names='Category',
                        title='Distribution of Drug Interaction Risk Levels')
            st.plotly_chart(fig, use_container_width=True)
        else:
            st.info("No drug interaction risk data available.")
    
    # Add download button for visualization data
    st.download_button(
        label="📥 Download Visualization Data",
        data=df.to_csv(index=False),
        file_name="visualization_data.csv",
        mime="text/csv"
    )

except Exception as e:
    st.error(f"An error occurred while creating visualizations: {str(e)}")
    st.info("Tip: Make sure the articles contain numerical data in a format that can be extracted.")

# Add tips for better visualization results
with st.expander("💡 Tips for Better Visualizations"):
    st.markdown("""
    To get better visualization results:
    1. Ensure articles contain numerical data in standard formats
    2. Look for percentages, time periods, and risk levels in the text
    3. Use consistent units across articles
    4. Include clear labels and measurements in the text
    """) 