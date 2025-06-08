from langchain.memory import ConversationBufferMemory
from langchain.chains import LLMChain
from langchain.prompts import PromptTemplate
from langchain_mistralai import ChatMistralAI
import os
from dotenv import load_dotenv

load_dotenv()
mistral_api_key = os.getenv("MISTRAL_API_KEY")

def get_mistral_with_memory():
    llm = ChatMistralAI(
        api_key=mistral_api_key,
        temperature=0.3
    )

    memory = ConversationBufferMemory(memory_key="chat_history", return_messages=True)

    prompt = PromptTemplate.from_template("""
You are a helpful and professional medical chatbot working with scientific articles from PubMed.
Use the user's query and the article content to figure out what task they're referring to
(e.g., summary, comparison, metadata, conclusion), and respond with structured, accurate information.

Chat history: {chat_history}
User input: {input}
""")

    return LLMChain(llm=llm, prompt=prompt, memory=memory)
