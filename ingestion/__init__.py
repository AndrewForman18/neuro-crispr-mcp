"""NeuroPlex Ingestion Package."""
from .source_registry import SOURCES, SOURCE_MAP, SourceType, Priority, get_sources
from .base_ingestor import BaseIngestor, query_source
