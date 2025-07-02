import re
import time
from typing import Dict, List, Any

import requests

# Vector search and embeddings
try:
    from qdrant_client import QdrantClient
    from qdrant_client.models import Distance, VectorParams, PointStruct, Batch
    import numpy as np
    VECTOR_SEARCH_AVAILABLE = True
except ImportError:
    VECTOR_SEARCH_AVAILABLE = False
    print("Vector search not available - install qdrant-client")

# MixedBread embeddings
try:
    import requests as embedding_requests
    MIXEDBREAD_AVAILABLE = True
except ImportError:
    MIXEDBREAD_AVAILABLE = False

# FastEmbed for local embeddings
try:
    from fastembed import TextEmbedding
    FASTEMBED_AVAILABLE = True
except ImportError:
    FASTEMBED_AVAILABLE = False
    print("FastEmbed not available - install fastembed package for local embeddings")

class MixedBreadEmbeddings:
    """MixedBread embeddings API integration with FastEmbed fallback"""

    def __init__(self, api_key: str = None, model: str = "mixedbread-ai/mxbai-embed-large-v1"):
        self.api_key = api_key
        self.model = model
        self.base_url = "https://api.mixedbread.ai/v1"
        self.local_embedding_model = None
        if not self.api_key and FASTEMBED_AVAILABLE:
            print("Initializing local FastEmbed model...")
            self.local_embedding_model = TextEmbedding()
            print("Local FastEmbed model initialized.")

    def embed_texts(self, texts: List[str]) -> List[List[float]]:
        """Generate embeddings using MixedBread API or FastEmbed (local fallback)"""
        if self.api_key:
            try:
                headers = {
                    "Authorization": f"Bearer {self.api_key}",
                    "Content-Type": "application/json"
                }

                payload = {
                    "model": self.model,
                    "input": texts,
                    "encoding_format": "float"
                }

                response = embedding_requests.post(
                    f"{self.base_url}/embeddings",
                    headers=headers,
                    json=payload,
                    timeout=30
                )

                if response.status_code == 200:
                    data = response.json()
                    embeddings = [item["embedding"] for item in data["data"]]
                    return embeddings
                else:
                    print(f"MixedBread API error: {response.status_code}. Falling back to local embeddings.")
                    if self.local_embedding_model:
                        return self.local_embedding_model.embed(texts).tolist()
                    else:
                        return [[float(i % 1000) / 1000 for i in range(1024)] for _ in texts]

            except Exception as e:
                print(f"MixedBread embedding error: {e}. Falling back to local embeddings.")
                if self.local_embedding_model:
                    return self.local_embedding_model.embed(texts).tolist()
                else:
                    return [[float(i % 1000) / 1000 for i in range(1024)] for _ in texts]
        else:
            if self.local_embedding_model:
                print("Using local FastEmbed model for embeddings.")
                return self.local_embedding_model.embed(texts).tolist()
            else:
                print("No MixedBread API key and FastEmbed not available - using mock embeddings.")
                return [[float(i % 1000) / 1000 for i in range(1024)] for _ in texts]

class EnhancedDocumentationIndexer:
    """Enhanced documentation indexer with MixedBread embeddings"""

    def __init__(self, mixedbread_api_key: str = None):
        self.qdrant_client = None
        self.collection_name = "spatial_docs_enhanced"
        self.chunk_size = 800
        self.chunk_overlap = 100
        self.embeddings_model = MixedBreadEmbeddings(mixedbread_api_key)

        if VECTOR_SEARCH_AVAILABLE:
            self.setup_qdrant()

    def setup_qdrant(self):
        """Initialize Qdrant client (in-memory for demo)"""
        try:
            self.qdrant_client = QdrantClient(":memory:")
            # Create collection for documentation
            self.qdrant_client.create_collection(
                collection_name=self.collection_name,
                vectors_config=VectorParams(size=1024, distance=Distance.COSINE)
            )
            print("✅ Qdrant vector database initialized")
        except Exception as e:
            print(f"Failed to setup Qdrant: {e}")

    def extract_text_from_html(self, html_content: str) -> str:
        """Extract clean text from HTML content"""
        # Remove scripts and style elements
        html_content = re.sub(r'<script[^>]*>.*?</script>', '', html_content, flags=re.DOTALL | re.IGNORECASE)
        html_content = re.sub(r'<style[^>]*>.*?</style>', '', html_content, flags=re.DOTALL | re.IGNORECASE)

        # Remove HTML tags
        text = re.sub(r'<[^>]+>', ' ', html_content)

        # Clean up whitespace
        text = re.sub(r'\s+', ' ', text).strip()

        # Remove common navigation/footer text
        text = re.sub(r'(Copyright|©|\d{4}|Terms|Privacy|Cookie)', '', text, flags=re.IGNORECASE)

        return text[:100000]  # Limit for demo

    def chunk_text_smart(self, text: str, url: str) -> List[Dict]:
        """Intelligently chunk text based on structure"""
        chunks = []

        # Split by common delimiters (headers, paragraphs)
        sections = re.split(r'\n\s*\n|\. (?=[A-Z])', text)

        current_chunk = ""
        chunk_id = 0

        for section in sections:
            section = section.strip()
            if not section:
                continue

            # If adding this section would exceed chunk size, save current chunk
            if len(current_chunk) + len(section) > self.chunk_size and current_chunk:
                chunks.append({
                    "text": current_chunk.strip(),
                    "url": url,
                    "chunk_id": chunk_id,
                    "char_count": len(current_chunk)
                })
                current_chunk = section
                chunk_id += 1
            else:
                current_chunk += f" {section}" if current_chunk else section

        # Add the last chunk
        if current_chunk.strip():
            chunks.append({
                "text": current_chunk.strip(),
                "url": url,
                "chunk_id": chunk_id,
                "char_count": len(current_chunk)
            })

        return chunks

    def scrape_documentation_enhanced(self, urls: List[str]) -> Dict[str, Any]:
        """Enhanced documentation scraping with better text extraction"""
        docs = {}

        headers = {
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36'
        }

        for url in urls:
            try:
                print(f"📄 Scraping: {url}")
                response = requests.get(url, headers=headers, timeout=15)

                if response.status_code == 200:
                    clean_text = self.extract_text_from_html(response.text)

                    # Extract title if possible
                    title_match = re.search(r'<title[^>]*>([^<]+)</title>', response.text, re.IGNORECASE)
                    title = title_match.group(1) if title_match else "Documentation"

                    docs[url] = {
                        "content": clean_text,
                        "title": title.strip(),
                        "timestamp": time.time(),
                        "char_count": len(clean_text),
                        "status": "success"
                    }
                    print(f"✅ {url}: {len(clean_text)} characters extracted")
                else:
                    docs[url] = {
                        "error": f"HTTP {response.status_code}",
                        "timestamp": time.time(),
                        "status": "error"
                    }
                    print(f"❌ {url}: HTTP {response.status_code}")

            except Exception as e:
                docs[url] = {
                    "error": str(e),
                    "timestamp": time.time(),
                    "status": "error"
                }
                print(f"❌ {url}: {str(e)}")

        return docs

    def index_documentation(self, progress_callback=None) -> str:
        """Index documentation with MixedBread embeddings"""
        if not self.qdrant_client:
            return "❌ Vector search not available - install qdrant-client"

        # Official documentation URLs
        urls = [
            "https://viash.io/docs/",
            "https://www.nextflow.io/docs/latest/",
            "https://openproblems.bio/documentation",
            "https://docs.docker.com/"
        ]

        results = []

        if progress_callback:
            progress_callback(0.1, "🌐 Scraping documentation...")

        docs = self.scrape_documentation_enhanced(urls)

        if progress_callback:
            progress_callback(0.3, "📝 Processing and chunking content...")

        all_chunks = []
        for url, doc_data in docs.items():
            if doc_data.get("status") == "success" and "content" in doc_data:
                chunks = self.chunk_text_smart(doc_data["content"], url)
                all_chunks.extend(chunks)
                results.append(f"✅ {doc_data.get('title', url)}: {len(chunks)} chunks ({doc_data['char_count']} chars)")
            else:
                results.append(f"❌ {url}: {doc_data.get('error', 'Failed')}")

        if progress_callback:
            progress_callback(0.6, "🧠 Generating MixedBread embeddings...")

        if all_chunks:
            # Generate embeddings with MixedBread
            texts = [f"{chunk['url']} | {chunk['text']}" for chunk in all_chunks]
            embeddings = self.embeddings_model.embed_texts(texts)

            if progress_callback:
                progress_callback(0.9, "💾 Storing in Qdrant vector database...")

            # Store in Qdrant using batch upsert
            points = []
            for i, (chunk, embedding) in enumerate(zip(all_chunks, embeddings)):
                points.append(PointStruct(
                    id=i,
                    vector=embedding,
                    payload={
                        **chunk,
                        "indexed_at": time.time()
                    }
                ))

            # Batch upsert for efficiency
            batch_size = 100
            for i in range(0, len(points), batch_size):
                batch = points[i:i + batch_size]
                self.qdrant_client.upsert(
                    collection_name=self.collection_name,
                    points=batch
                )

            results.append(f"\n🎉 Successfully indexed {len(all_chunks)} chunks with MixedBread embeddings!")
            results.append(f"📊 Vector database contains {len(embeddings)} vectors")

        if progress_callback:
            progress_callback(1.0, "✅ Documentation indexing complete!")

        return "\n".join(results)

    def search_documentation(self, query: str, limit: int = 5) -> List[Dict]:
        """Search indexed documentation using MixedBread embeddings"""
        if not self.qdrant_client:
            return [{"text": "Vector search not available", "url": "", "score": 0}]

        try:
            # Generate query embedding
            query_embedding = self.embeddings_model.embed_texts([query])[0]

            # Search in Qdrant
            search_results = self.qdrant_client.search(
                collection_name=self.collection_name,
                query_vector=query_embedding,
                limit=limit,
                with_payload=True
            )

            results = []
            for result in search_results:
                results.append({
                    "text": result.payload["text"][:600] + "..." if len(result.payload["text"]) > 600 else result.payload["text"],
                    "url": result.payload["url"],
                    "score": float(result.score),
                    "chunk_id": result.payload.get("chunk_id", 0),
                    "char_count": result.payload.get("char_count", 0)
                })

            return results
        except Exception as e:
            return [{"text": f"Search error: {str(e)}", "url": "", "score": 0}]