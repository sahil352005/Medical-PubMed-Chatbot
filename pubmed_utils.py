import re
from typing import List, Dict
from Bio import Entrez, Medline

# Set your email (required by Entrez)
Entrez.email = "sohambagayatkar7@gmail.com"  # ✅ Use your real email

# ----------- Utility Functions -----------

def extract_pmids_from_urls(urls: List[str]) -> List[str]:
    """Extract PMIDs from PubMed URLs"""
    pmid_pattern = re.compile(r"pubmed\.ncbi\.nlm\.nih\.gov/(\d+)/?")
    pmids = [match.group(1) for url in urls if (match := pmid_pattern.search(url))]
    return pmids

def extract_pmids_from_text(text: str) -> List[str]:
    """Extract PMIDs from text content (e.g., saved PubMed results)"""
    # Pattern for PMID in various formats
    patterns = [
        r"PMID: (\d+)",  # Standard PMID format
        r"PMID-(\d+)",   # Alternative format
        r"\[PMID: (\d+)\]",  # PMID in brackets
        r"pubmed/(\d+)",  # URL format
        r"(\d{8})"       # Just the number (8 digits)
    ]
    
    pmids = set()
    for pattern in patterns:
        matches = re.findall(pattern, text)
        pmids.update(matches)
    
    return list(pmids)

def fetch_pubmed_articles(pmids: List[str]) -> List[Dict]:
    """Fetch PubMed articles using PMIDs and return metadata"""
    if len(pmids) > 10:
        raise ValueError("You can only process up to 10 articles at a time.")

    with Entrez.efetch(db="pubmed", id=",".join(pmids), rettype="medline", retmode="text") as handle:
        records = Medline.parse(handle)
        articles = []
        for record in records:
            article = {
                "pmid": record.get("PMID", ""),
                "title": record.get("TI", ""),
                "authors": record.get("AU", []),
                "journal": record.get("JT", ""),
                "abstract": record.get("AB", ""),
                "date": record.get("DP", ""),
            }
            articles.append(article)
        return articles

def get_articles_from_pmids_or_urls(inputs: List[str]) -> List[Dict]:
    """Process user inputs (PMIDs or URLs) and fetch article data"""
    # Separate valid PMIDs and extract PMIDs from URLs
    pmids = [x.strip() for x in inputs if x.strip().isdigit()]
    urls = [x.strip() for x in inputs if "pubmed.ncbi.nlm.nih.gov" in x]
    pmids_from_urls = extract_pmids_from_urls(urls)
    all_pmids = list(set(pmids + pmids_from_urls))[:10]  # max 10 articles

    if not all_pmids:
        raise ValueError("No valid PMIDs or PubMed URLs provided.")

    return fetch_pubmed_articles(all_pmids)

def search_pubmed_by_keyword(keyword: str, max_results: int = 10) -> List[str]:
    """Search PubMed using a keyword and return a list of PMIDs"""
    handle = Entrez.esearch(db="pubmed", term=keyword, retmax=max_results)
    record = Entrez.read(handle)
    handle.close()
    return record["IdList"]
