import os
from typing import List
from dotenv import load_dotenv
from langchain.vectorstores import FAISS
from langchain_mistralai import MistralAIEmbeddings  
from langchain.docstore.document import Document

load_dotenv()
mistral_api_key = os.getenv("MISTRAL_API_KEY")

# ----------- Initialize Embeddings -----------

embedding_model = MistralAIEmbeddings(
    model="mistral-embed",
    api_key=mistral_api_key
)

# ----------- Build Vector Store -----------

def build_faiss_vectorstore(articles: List[dict]) -> FAISS:
    """
    Convert articles to FAISS vector store.
    Each article abstract becomes a document.
    """
    docs = []
    for i, art in enumerate(articles):
        content = f"Title: {art['title']}\nAbstract: {art['abstract']}"
        metadata = {
            "pmid": art.get("pmid", ""),
            "source": f"Article-{i+1}"
        }
        docs.append(Document(page_content=content, metadata=metadata))

    vectorstore = FAISS.from_documents(docs, embedding_model)
    return vectorstore


# ----------- Query Vector Store -----------

def query_vectorstore(vectorstore: FAISS, user_query: str, k: int = 3) -> List[Document]:
    """
    Retrieve top k relevant documents from vector store for user query
    """
    return vectorstore.similarity_search(user_query, k=k)