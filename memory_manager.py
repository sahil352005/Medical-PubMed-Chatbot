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

    memory = ConversationBufferMemory(
        memory_key="chat_history",
        return_messages=True,
        input_key="human_input"
    )

    prompt = PromptTemplate.from_template("""
You are a helpful and professional medical chatbot working with scientific articles from PubMed.
Use the user's query and the article content to provide accurate, evidence-based responses.

Context from PubMed articles:
{articles_text}

Previous conversation:
{chat_history}

Current user query: {human_input}

Instructions:
1. If the query is about the articles, use the provided article context to answer
2. If the query is a follow-up question, use the chat history to maintain context
3. Always maintain a professional medical tone
4. Cite specific articles when relevant
5. If you're unsure about something, acknowledge the limitations

Your response:
""")

    return LLMChain(
        llm=llm,
        prompt=prompt,
        memory=memory,
        verbose=True
    )