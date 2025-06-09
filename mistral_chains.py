from langchain_mistralai import ChatMistralAI
from langchain.prompts import PromptTemplate
from langchain.chains import LLMChain
import os
from dotenv import load_dotenv

load_dotenv()
mistral_api_key = os.getenv("MISTRAL_API_KEY")

llm = ChatMistralAI(api_key=mistral_api_key, temperature=0.3)

# --- Classify query ---
def classify_query(user_query: str) -> list[str]:
    prompt_template = PromptTemplate.from_template(
        """
You are an expert medical writer with extensive experience in evidence synthesis, clinical protocol development (ICH M11 format), and regulatory-standard documentation (TransCelerate CSR format) and a classifier that analyzes a medical query and returns the relevant task labels from the following options:

1. overall_summary
2. article_metadata
3. per_article_summary
4. comparison
5. conclusion

If the query is broad (e.g., "Give all insights" or "Summarize everything"), return all the above task labels separated by commas.

Return a comma-separated list of relevant task labels (e.g., "overall_summary, comparison"). Do not explain. Only return labels.

Query: "{query}"
Tasks:
"""
    )
    chain = LLMChain(llm=llm, prompt=prompt_template)
    result = chain.run({"query": user_query}).strip().lower()
    tasks = [task.strip() for task in result.split(",") if task.strip()]
    return tasks

# --- Task Functions ---

def run_overall_summary(article_text: str) -> str:
    prompt = f"""
You are a medical expert. Read the following abstracts from multiple articles and write a single integrated summary in about 150–200 words. Do not summarize each article separately. Instead, synthesize the key findings across all studies and present a unified perspective, mentioning similarities, differences, and overall conclusions.

Text:
{article_text}
"""
    response = llm.invoke(prompt)
    return response.content if hasattr(response, "content") else response

def run_per_article_summary(articles: list[dict]) -> str:
    summaries = []
    for i, art in enumerate(articles, 1):
        prompt = f"""
You are a medical reviewer. Generate a ~100-word structured summary for the article below.

Use the following format:

- *Background:* Brief context of the study.
- *Methodology:* Describe the study design or analysis method.
- *Key Findings:* Summarize the main results.
- *Conclusion:* Clinical or scientific implications.

Article:
Title: {art['title']}
Authors: {', '.join(art['authors'])}
Abstract: {art['abstract']}
"""
        response = llm.invoke(prompt)
        result = response.content if hasattr(response, "content") else response
        summaries.append(f"### 📝 Article {i}: {art['title']}\n{result.strip()}")

    return "\n\n".join(summaries)

def run_conclusion(article_text: str) -> str:
    prompt = f"""
You are a medical reviewer. Based on the following articles, generate an overall evidence-based conclusion. Your response should:

- Compare the main interventions or treatments discussed
- Highlight key differences in efficacy or safety
- Mention any clinical advantages or disadvantages
- Offer a brief clinical guidance or recommendation based on the collective evidence

Do not summarize each article separately. Write a unified, professional conclusion suitable for a medical report or systematic review.

Articles:
{article_text}
"""
    response = llm.invoke(prompt)
    return response.content if hasattr(response, "content") else response

def run_comparison(article_text: str) -> str:
    prompt = f"""
You are a medical researcher. Compare the following articles using 10 key medical parameters such as:
- Treatment
- Mechanism of Action
- Indication
- Onset of Action
- Healing Rate
- Bleeding Risk
- CYP Interactions
- GI Side Effects
- Trial Design
- Population

Format your response as a Markdown table with one row per article and one column per parameter.

Articles:
{article_text}
"""
    response = llm.invoke(prompt)
    return response.content if hasattr(response, "content") else response

def run_article_metadata(articles: list[dict]) -> str:
    metadata_str = ""
    for art in articles:
        metadata_str += (
            f"*PMID:* {art.get('pmid', '')}\n"
            f"*Title:* {art.get('title', '')}\n"
            f"*Authors:* {', '.join(art.get('authors', []))}\n"
            f"*Journal:* {art.get('journal', '')}\n"
            f"*Date:* {art.get('date', '')}\n"
            f"*Abstract:* {art.get('abstract', '')[:300]}...\n\n"
        )
    return metadata_str

# --- Master Handler Function ---

def answer_user_query(user_query: str, articles: list[dict], articles_text: str, forced_tasks: list[str] = None) -> str:
    if forced_tasks is not None:
        tasks = forced_tasks
    else:
        tasks = classify_query(user_query)

    if not tasks:
        return "❌ Couldn't classify your query properly. Please rephrase."

    responses = []

    if "overall_summary" in tasks:
        responses.append("### 🟡  Overall Summary\n" + run_overall_summary(articles_text))

    if "article_metadata" in tasks:
        responses.append("### 🔵  Article Metadata\n" + run_article_metadata(articles))

    if "per_article_summary" in tasks:
        responses.append("### 🟢  Per-Article Summary\n" + run_per_article_summary(articles))


    if "comparison" in tasks:
        responses.append("### 🔴  Comparison on 10 Medical Parameters\n" + run_comparison(articles_text))

    if "conclusion" in tasks:
        responses.append("### 🟣  Final Conclusion\n" + run_conclusion(articles_text))

    return "\n\n".join(responses)