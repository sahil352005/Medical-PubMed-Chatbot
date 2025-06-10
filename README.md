
# 🩺 Medical Research Article Summarizer & Regulatory Content Generation

This tool helps medical professionals quickly **summarize PubMed research articles** and **generate regulatory content** (like clinical trial protocols and study reports). It uses **Generative AI** powered by the **Mistral API** to provide structured insights and draft compliant documents.

---
## Glipmse

![image](https://github.com/user-attachments/assets/3f913aaa-a397-4b24-812d-192c932ae1b4)

![image](https://github.com/user-attachments/assets/0991936f-445f-4826-bead-d945cc1eccbb)


## How It Works

1.  **Search & Save on PubMed:** Find articles on PubMed and save the results as a text file.
2.  **Upload & Process:** Upload the saved file to the "Article Summaries" page in the tool.
3.  **Get Summaries:** The tool generates an **overall evidence summary, metadata table, comparative analysis, and clinical conclusions**.
4.  **Generate Regulatory Content:** On the "Regulatory Content Generation" page, get drafts for **Clinical Trial Protocol Introductions (ICH M11)** and **Clinical Study Report Discussions (TransCelerate)**.

---

## Key Features

* **Fast Medical Summaries:** Quickly understand research papers.
* **Automated Regulatory Drafts:** Speeds up documentation for trials and reports.
* **AI-Powered:** Uses Mistral API for smart content generation.
* **User-Friendly:** Easy-to-follow steps in a Streamlit interface.

---

## Getting Started

To run this tool locally, follow these concise steps:

1.  **Clone the Repository:**
    ```bash
    git clone https://github.com/sahil352005/Medical-PubMed-Chatbot.git
    cd medical-research-summarizer
    ```
2.  **Set Up Environment:** Create and activate a virtual environment.
    ```bash
    python -m venv venv
    source venv/bin/activate # Use `.\venv\Scripts\activate` on Windows
    ```
3.  **Install Dependencies:** Install required libraries.
    ```bash
    pip install -r requirements.txt
    ```
4.  **Configure API Key:** Create a `.env` file in the project root and add your Mistral API key:
    ```
    MISTRAL_API_KEY="your_mistral_api_key"
    ```
5.  **Run the Tool:**
    ```bash
    streamlit run app.py
    ```

---




