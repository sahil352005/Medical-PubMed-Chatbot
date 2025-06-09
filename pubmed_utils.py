import re
from typing import List, Dict
from Bio import Entrez, Medline

# Set your email (required by Entrez)
Entrez.email = "sohambagayatkar7@gmail.com"  # ✅ Replace with your real email if needed

# ----------- Utility Functions -----------

def extract_pmids_from_urls(urls: List[str]) -> List[str]:
    """Extract PMIDs from PubMed URLs"""
    pmid_pattern = re.compile(r"pubmed\.ncbi\.nlm\.nih\.gov/(\d+)/?")
    pmids = [match.group(1) for url in urls if (match := pmid_pattern.search(url))]
    return pmids

def extract_pmids_from_text(text: str) -> List[str]:
    """Extract PMIDs from text content (e.g., saved PubMed results)"""
    patterns = [
        r"PMID: (\d+)",
        r"PMID-(\d+)",
        r"\[PMID: (\d+)\]",
        r"pubmed/(\d+)",
        r"https?://pubmed\.ncbi\.nlm\.nih\.gov/(\d+)/?",
        r"\b(\d{8})\b"  # Match 8-digit standalone PMIDs
    ]
    
    pmids = set()
    for pattern in patterns:
        matches = re.findall(pattern, text)
        pmids.update(matches)
    
    return list(pmids)

def fetch_pubmed_articles(pmids: List[str]) -> List[Dict]:
    """Fetch PubMed articles using PMIDs and return metadata"""
    if not pmids:
        raise ValueError("No PMIDs provided.")
    
    pmids = pmids[:10]  # Ensure only 10 are processed

    with Entrez.efetch(db="pubmed", id=",".join(pmids), rettype="medline", retmode="text") as handle:
        records = Medline.parse(handle)
        articles = []
        for record in records:
            pmid = record.get("PMID", "")
            article = {
                "pmid": pmid,
                "title": record.get("TI", "No title available"),
                "authors": record.get("AU", ["No authors listed"]),
                "journal": record.get("JT", "No journal"),
                "abstract": record.get("AB", "No abstract available"),
                "date": record.get("DP", "No date"),
                "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid else ""
            }
            articles.append(article)
        return articles

def get_articles_from_pmids_or_urls(inputs: List[str]) -> List[Dict]:
    """Process user inputs (PMIDs or URLs) and fetch article data"""
    pmids = [x.strip() for x in inputs if x.strip().isdigit()]
    urls = [x.strip() for x in inputs if "pubmed.ncbi.nlm.nih.gov" in x]
    pmids_from_urls = extract_pmids_from_urls(urls)
    
    all_pmids = list(set(pmids + pmids_from_urls))[:10]

    if not all_pmids:
        raise ValueError("No valid PMIDs or PubMed URLs found in input.")

    return fetch_pubmed_articles(all_pmids)

def search_pubmed_by_keyword(keyword: str, max_results: int = 10) -> List[str]:
    """Search PubMed using a keyword and return a list of PMIDs"""
    handle = Entrez.esearch(db="pubmed", term=keyword, retmax=max_results)
    record = Entrez.read(handle)
    handle.close()
    return record.get("IdList", [])