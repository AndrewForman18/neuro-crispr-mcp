"""NeuroPlex Source Registry — Central configuration for all data sources.

This module defines every external data source that NeuroPlex can ingest and query.
Sources are grouped by Type (Literature, Expression, Genetics, Pathways, Druggability)
and each entry carries the metadata needed for:
  1. Ingestion pipelines (API endpoints, rate limits, schema mapping)
  2. Tool generation (auto-registers MCP tools per source)
  3. UI rendering (sidebar checkboxes, labels, icons)

Target catalog/schema are loaded from config/neuroplex_env.yml
Table prefix: neuroplex_
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Optional

from config.neuroplex_config import load_config


# ═══════════════════════════════════════════════════════════════════════════════
# Enums
# ═══════════════════════════════════════════════════════════════════════════════

class SourceType(str, Enum):
    """Data source category — maps to sidebar dataset groups."""
    LITERATURE = "literature"
    EXPRESSION = "expression"
    GENETICS = "genetics"
    PATHWAYS = "pathways"
    DRUGGABILITY = "druggability"


class Priority(str, Enum):
    HIGH = "high"
    MEDIUM = "medium"
    LOW = "low"


class IngestMethod(str, Enum):
    """How data is fetched from the source."""
    REST_API = "rest_api"          # Standard REST endpoints
    BULK_DOWNLOAD = "bulk_download"  # FTP/HTTP file downloads (TSV, XML)
    GRAPHQL = "graphql"            # GraphQL API (e.g., OpenTargets)
    SCRAPE = "scrape"              # Structured HTML/XML parsing
    SPARQL = "sparql"              # RDF/SPARQL endpoint (e.g., UniProt)


# ═══════════════════════════════════════════════════════════════════════════════
# Source Definition
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class SourceConfig:
    """Configuration for a single NeuroPlex data source."""
    # Identity
    key: str                        # Snake_case identifier (e.g., "pubmed")
    name: str                       # Display name (e.g., "PubMed")
    source_type: SourceType
    priority: Priority
    description: str                # "Best when you want" column

    # Ingestion
    ingest_method: IngestMethod
    base_url: str                   # Primary API/download URL
    rate_limit_rpm: int = 60        # Requests per minute
    bulk_url: Optional[str] = None  # Alternate bulk download URL
    auth_required: bool = False     # Whether API key/token is needed
    auth_secret_key: Optional[str] = None  # Databricks secret scope key

    # Storage
    target_table: str = ""          # Auto-derived: neuroplex_{key}
    partition_cols: list[str] = field(default_factory=list)
    incremental: bool = True        # Supports incremental ingestion

    # Tool generation
    tool_name: str = ""             # Auto-derived: query_{key}
    tool_description: str = ""      # LLM-facing tool description
    searchable_fields: list[str] = field(default_factory=list)

    # UI
    icon: str = "\U0001F4C4"        # Sidebar emoji
    i18n_key: str = ""              # Key into _I18N dict

    def __post_init__(self):
        if not self.target_table:
            self.target_table = f"neuroplex_{self.key}"
        if not self.tool_name:
            self.tool_name = f"query_{self.key}"
        if not self.i18n_key:
            self.i18n_key = f"ds_{self.key}"


# ═══════════════════════════════════════════════════════════════════════════════
# Source Registry
# ═══════════════════════════════════════════════════════════════════════════════

SOURCES: list[SourceConfig] = [
    # ──────────────────────────────────────────────────────────────────────────
    # LITERATURE
    # ──────────────────────────────────────────────────────────────────────────
    SourceConfig(
        key="pubmed",
        name="PubMed",
        source_type=SourceType.LITERATURE,
        priority=Priority.HIGH,
        description="Article search, PubTator annotations, and PMC full-text handoff",
        ingest_method=IngestMethod.REST_API,
        base_url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/",
        rate_limit_rpm=180,  # 3/sec with API key
        auth_required=True,
        auth_secret_key="ncbi_api_key",
        searchable_fields=["gene_symbol", "disease", "drug", "pmid"],
        tool_description="[LITERATURE] Search PubMed for neuroscience publications. Returns article titles, abstracts, MeSH terms, and PubTator gene/disease annotations.",
        icon="\U0001F4DA",
    ),
    SourceConfig(
        key="clinicaltrials",
        name="ClinicalTrials.gov",
        source_type=SourceType.LITERATURE,
        priority=Priority.HIGH,
        description="Recruiting-study search, eligibility text, and site details",
        ingest_method=IngestMethod.REST_API,
        base_url="https://clinicaltrials.gov/api/v2/",
        rate_limit_rpm=120,
        searchable_fields=["condition", "intervention", "sponsor", "nct_id"],
        tool_description="[LITERATURE] Search clinical trials by condition, drug, or gene target. Returns trial phase, status, eligibility, sites, and endpoints.",
        icon="\U0001F3E5",
    ),
    SourceConfig(
        key="clinvar",
        name="ClinVar",
        source_type=SourceType.LITERATURE,
        priority=Priority.HIGH,
        description="Clinical significance and review-status context for variants",
        ingest_method=IngestMethod.REST_API,
        base_url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/",
        rate_limit_rpm=180,
        auth_required=True,
        auth_secret_key="ncbi_api_key",
        searchable_fields=["gene_symbol", "variant", "condition", "clinical_significance"],
        tool_description="[LITERATURE] Query ClinVar for variant clinical significance, review status, and associated conditions. Links variants to diseases.",
        icon="\U0001F9EC",
    ),
    SourceConfig(
        key="openfda",
        name="OpenFDA",
        source_type=SourceType.LITERATURE,
        priority=Priority.HIGH,
        description="FAERS, recalls, device events, labels, and U.S. approval context",
        ingest_method=IngestMethod.REST_API,
        base_url="https://api.fda.gov/",
        rate_limit_rpm=240,  # 4/sec without key
        searchable_fields=["drug_name", "reaction", "indication", "application_number"],
        tool_description="[LITERATURE] Search FDA adverse events (FAERS), drug labels, approvals, and recalls. Covers orexin antagonists and NDD therapeutics.",
        icon="\U0001F3DB\uFE0F",
    ),
    SourceConfig(
        key="semantic_scholar",
        name="Semantic Scholar",
        source_type=SourceType.LITERATURE,
        priority=Priority.HIGH,
        description="TLDRs, citation graphs, references, and recommendations",
        ingest_method=IngestMethod.REST_API,
        base_url="https://api.semanticscholar.org/graph/v1/",
        rate_limit_rpm=100,
        auth_required=True,
        auth_secret_key="semantic_scholar_api_key",
        searchable_fields=["query", "paper_id", "author", "fields_of_study"],
        tool_description="[LITERATURE] Search Semantic Scholar for papers with TLDRs, citation counts, influential citations, and related papers.",
        icon="\U0001F4D6",
    ),
    SourceConfig(
        key="civic",
        name="CIViC",
        source_type=SourceType.GENETICS,
        priority=Priority.HIGH,
        description="Clinical variant evidence, therapy context, and disease-associated variants",
        ingest_method=IngestMethod.REST_API,
        base_url="https://civicdb.org/api/graphql",
        rate_limit_rpm=60,
        searchable_fields=["gene", "variant", "disease", "therapy", "evidence_level"],
        tool_description="[GENETICS] Query CIViC for clinical evidence on variants — therapy associations, evidence levels, and variant interpretations.",
        icon="\U0001F9EA",
    ),
    SourceConfig(
        key="cbioportal",
        name="cBioPortal",
        source_type=SourceType.GENETICS,
        priority=Priority.HIGH,
        description="Cancer cohort frequencies and local study analytics workflows",
        ingest_method=IngestMethod.REST_API,
        base_url="https://www.cbioportal.org/api/",
        rate_limit_rpm=60,
        searchable_fields=["gene", "study_id", "cancer_type", "mutation_type"],
        tool_description="[GENETICS] Query cBioPortal for cancer mutation frequencies, co-occurrence, and survival data by gene and cancer type.",
        icon="\U0001F52C",
    ),
    SourceConfig(
        key="ema",
        name="EMA",
        source_type=SourceType.LITERATURE,
        priority=Priority.HIGH,
        description="EU regulatory, safety, and shortage context for medicines",
        ingest_method=IngestMethod.REST_API,
        base_url="https://www.ema.europa.eu/en/medicines/",
        rate_limit_rpm=30,
        searchable_fields=["drug_name", "active_substance", "therapeutic_area", "authorisation_status"],
        tool_description="[LITERATURE] Search EMA for EU medicine authorisations, safety signals, shortages, and EPAR summaries.",
        icon="\U0001F1EA\U0001F1FA",
    ),
    SourceConfig(
        key="ncbi_gtr",
        name="NCBI Genetic Testing Registry",
        source_type=SourceType.GENETICS,
        priority=Priority.HIGH,
        description="Gene-centric genetic tests, GTR diagnostic cards, and local bundle lifecycle",
        ingest_method=IngestMethod.REST_API,
        base_url="https://www.ncbi.nlm.nih.gov/gtr/",
        rate_limit_rpm=180,
        auth_required=True,
        auth_secret_key="ncbi_api_key",
        searchable_fields=["gene", "condition", "test_name", "lab_name"],
        tool_description="[GENETICS] Search NCBI GTR for available genetic tests by gene or condition. Returns lab info, test methods, and clinical utility.",
        icon="\U0001F9EA",
    ),
    SourceConfig(
        key="medlineplus",
        name="MedlinePlus",
        source_type=SourceType.LITERATURE,
        priority=Priority.HIGH,
        description="Plain-language disease/symptom context for discover and disease clinical_features",
        ingest_method=IngestMethod.REST_API,
        base_url="https://connect.medlineplus.gov/service",
        rate_limit_rpm=60,
        searchable_fields=["disease", "gene", "symptom"],
        tool_description="[LITERATURE] Query MedlinePlus for plain-language disease descriptions, symptoms, genetic conditions, and patient-facing summaries.",
        icon="\U0001F4CB",
    ),
    SourceConfig(
        key="pharmgkb",
        name="PharmGKB / CPIC",
        source_type=SourceType.LITERATURE,
        priority=Priority.HIGH,
        description="Pharmacogenomic recommendations, frequencies, and clinical annotations",
        ingest_method=IngestMethod.REST_API,
        base_url="https://api.pharmgkb.org/v1/",
        rate_limit_rpm=60,
        searchable_fields=["gene", "drug", "variant", "phenotype"],
        tool_description="[LITERATURE] Query PharmGKB for pharmacogenomic drug-gene interactions, CPIC guidelines, dosing recommendations, and variant annotations.",
        icon="\U0001F48A",
    ),

    # ──────────────────────────────────────────────────────────────────────────
    # EXPRESSION
    # ──────────────────────────────────────────────────────────────────────────
    SourceConfig(
        key="uniprot",
        name="UniProt",
        source_type=SourceType.EXPRESSION,
        priority=Priority.HIGH,
        description="Canonical protein cards and structure-linked context",
        ingest_method=IngestMethod.REST_API,
        base_url="https://rest.uniprot.org/uniprotkb/",
        rate_limit_rpm=600,  # Generous limits
        searchable_fields=["gene", "protein_name", "organism", "function", "subcellular_location"],
        tool_description="[EXPRESSION] Query UniProt for protein function, structure, subcellular localization, domains, and post-translational modifications.",
        icon="\U0001F9F1",
    ),
    SourceConfig(
        key="human_protein_atlas",
        name="Human Protein Atlas",
        source_type=SourceType.EXPRESSION,
        priority=Priority.HIGH,
        description="Tissue expression, localization, and cancer-expression context",
        ingest_method=IngestMethod.REST_API,
        base_url="https://www.proteinatlas.org/api/",
        rate_limit_rpm=60,
        searchable_fields=["gene", "tissue", "cell_type", "cancer_type"],
        tool_description="[EXPRESSION] Query Human Protein Atlas for protein tissue expression, subcellular localization, single-cell type expression, and cancer specificity.",
        icon="\U0001F9CD",
    ),
    SourceConfig(
        key="monarch",
        name="Monarch Initiative",
        source_type=SourceType.EXPRESSION,
        priority=Priority.HIGH,
        description="Phenotype-to-disease matching, disease genes, and model evidence",
        ingest_method=IngestMethod.REST_API,
        base_url="https://api.monarchinitiative.org/v3/api/",
        rate_limit_rpm=120,
        searchable_fields=["gene", "disease", "phenotype", "model_organism"],
        tool_description="[EXPRESSION] Query Monarch Initiative for gene-disease associations, phenotype matching, model organism evidence, and gene-phenotype links.",
        icon="\U0001F451",
    ),

    # ──────────────────────────────────────────────────────────────────────────
    # GENETICS
    # ──────────────────────────────────────────────────────────────────────────
    SourceConfig(
        key="gnomad",
        name="gnomAD",
        source_type=SourceType.GENETICS,
        priority=Priority.HIGH,
        description="Population frequency and gene constraint context",
        ingest_method=IngestMethod.REST_API,
        base_url="https://gnomad.broadinstitute.org/api/",
        rate_limit_rpm=60,
        searchable_fields=["gene", "variant", "region", "population"],
        tool_description="[GENETICS] Query gnomAD for allele frequencies, gene constraint scores (pLI, LOEUF), and population-specific variant data.",
        icon="\U0001F4CA",
    ),
    SourceConfig(
        key="seer",
        name="SEER Explorer",
        source_type=SourceType.GENETICS,
        priority=Priority.HIGH,
        description="Cancer survival statistics and disease survival section output",
        ingest_method=IngestMethod.BULK_DOWNLOAD,
        base_url="https://seer.cancer.gov/data/",
        rate_limit_rpm=10,
        searchable_fields=["cancer_site", "stage", "age_group", "year"],
        tool_description="[GENETICS] Query SEER cancer statistics for incidence, survival rates, and trends by cancer type, stage, age, and demographics.",
        icon="\U0001F4C8",
    ),

    # ──────────────────────────────────────────────────────────────────────────
    # PATHWAYS
    # ──────────────────────────────────────────────────────────────────────────
    SourceConfig(
        key="reactome",
        name="Reactome",
        source_type=SourceType.PATHWAYS,
        priority=Priority.HIGH,
        description="Pathway records, pathway genes, and contained events",
        ingest_method=IngestMethod.REST_API,
        base_url="https://reactome.org/ContentService/",
        rate_limit_rpm=120,
        searchable_fields=["gene", "pathway_name", "pathway_id", "species"],
        tool_description="[PATHWAYS] Query Reactome for biological pathways containing a gene. Returns pathway hierarchy, participating molecules, and reaction steps.",
        icon="\U0001F504",
    ),
    SourceConfig(
        key="kegg",
        name="KEGG",
        source_type=SourceType.PATHWAYS,
        priority=Priority.HIGH,
        description="KEGG pathway IDs, summary cards, and pathway genes",
        ingest_method=IngestMethod.REST_API,
        base_url="https://rest.kegg.jp/",
        rate_limit_rpm=60,
        searchable_fields=["gene", "pathway_id", "compound", "disease"],
        tool_description="[PATHWAYS] Query KEGG for metabolic and signaling pathways, gene-pathway mappings, and compound-pathway links.",
        icon="\U0001F5FA\uFE0F",
    ),

    # ──────────────────────────────────────────────────────────────────────────
    # DRUGGABILITY
    # ──────────────────────────────────────────────────────────────────────────
    SourceConfig(
        key="opentargets",
        name="OpenTargets",
        source_type=SourceType.DRUGGABILITY,
        priority=Priority.HIGH,
        description="Target-disease scores, druggability, and disease-gene evidence",
        ingest_method=IngestMethod.GRAPHQL,
        base_url="https://api.platform.opentargets.org/api/v4/graphql",
        rate_limit_rpm=120,
        searchable_fields=["gene", "disease", "drug", "tractability"],
        tool_description="[DRUGGABILITY] Query OpenTargets for target-disease association scores, druggability tractability, known drugs, and evidence sources.",
        icon="\U0001F3AF",
    ),
    SourceConfig(
        key="ddinter",
        name="DDInter",
        source_type=SourceType.DRUGGABILITY,
        priority=Priority.MEDIUM,
        description="Structured drug-drug interactions, severity levels, and class-oriented partner review",
        ingest_method=IngestMethod.REST_API,
        base_url="https://ddinter.scbdd.com/api/",
        rate_limit_rpm=30,
        searchable_fields=["drug_a", "drug_b", "severity", "mechanism"],
        tool_description="[DRUGGABILITY] Query DDInter for drug-drug interactions, severity classifications, and mechanistic context for polypharmacy review.",
        icon="\u26A0\uFE0F",
    ),

    # ──────────────────────────────────────────────────────────────────────────
    # MEDIUM PRIORITY — LITERATURE
    # ──────────────────────────────────────────────────────────────────────────
    SourceConfig(
        key="who_ivd",
        name="WHO Prequalified IVD",
        source_type=SourceType.LITERATURE,
        priority=Priority.MEDIUM,
        description="Infectious-disease diagnostic products, assay formats, and WHO product-card provenance",
        ingest_method=IngestMethod.REST_API,
        base_url="https://extranet.who.int/pqweb/vitro-diagnostics/",
        rate_limit_rpm=20,
        searchable_fields=["disease", "product_name", "manufacturer"],
        tool_description="[LITERATURE] Search WHO Prequalified IVD list for diagnostic products, assay formats, and manufacturer details.",
        icon="\U0001F9EA",
    ),
]


# ═══════════════════════════════════════════════════════════════════════════════
# Registry Access Helpers
# ═══════════════════════════════════════════════════════════════════════════════

# Keyed lookup
SOURCE_MAP: dict[str, SourceConfig] = {s.key: s for s in SOURCES}

# Grouped by type
SOURCES_BY_TYPE: dict[SourceType, list[SourceConfig]] = {}
for _s in SOURCES:
    SOURCES_BY_TYPE.setdefault(_s.source_type, []).append(_s)

# Grouped by priority
SOURCES_BY_PRIORITY: dict[Priority, list[SourceConfig]] = {}
for _s in SOURCES:
    SOURCES_BY_PRIORITY.setdefault(_s.priority, []).append(_s)


def get_sources(source_type: SourceType | None = None, priority: Priority | None = None) -> list[SourceConfig]:
    """Filter sources by type and/or priority."""
    results = SOURCES
    if source_type:
        results = [s for s in results if s.source_type == source_type]
    if priority:
        results = [s for s in results if s.priority == priority]
    return results


CFG = load_config()

def get_target_tables(catalog: str = CFG.catalog, schema: str = CFG.schema) -> dict[str, str]:
    """Return fully qualified table names for all sources."""
    return {s.key: f"{catalog}.{schema}.{s.target_table}" for s in SOURCES}


# ═══════════════════════════════════════════════════════════════════════════════
# Tool Schema Generation (for dynamic MCP registration)
# ═══════════════════════════════════════════════════════════════════════════════

def generate_tool_schema(source: SourceConfig) -> dict:
    """Generate an OpenAI-compatible tool schema from a SourceConfig."""
    properties = {}
    required = []

    # Every source gets a generic query/search parameter
    properties["query"] = {
        "type": "string",
        "description": f"Search query for {source.name}. Can be a gene symbol, drug name, disease, or free text.",
    }
    required.append("query")

    # Add source-specific searchable fields as optional filters
    for field_name in source.searchable_fields:
        if field_name != "query":
            properties[field_name] = {
                "type": "string",
                "description": f"Filter by {field_name.replace('_', ' ')}",
            }

    # Limit parameter
    properties["limit"] = {
        "type": "integer",
        "description": "Maximum number of results to return (default: 20)",
    }

    return {
        "type": "function",
        "function": {
            "name": source.tool_name,
            "description": source.tool_description,
            "parameters": {
                "type": "object",
                "properties": properties,
                "required": required,
            },
        },
    }


def generate_all_tool_schemas(source_type: SourceType | None = None) -> list[dict]:
    """Generate tool schemas for all (or filtered) sources."""
    sources = get_sources(source_type=source_type)
    return [generate_tool_schema(s) for s in sources]


# ═══════════════════════════════════════════════════════════════════════════════
# Dataset → Tool Mapping (for hard-filtering in app.py)
# ═══════════════════════════════════════════════════════════════════════════════

def get_tool_dataset_map() -> dict[str, str]:
    """Return tool_name → source_type.value mapping for tool filtering.

    Used by app.py to hard-filter which tools are presented to the LLM
    based on the user's sidebar dataset focus selections.
    """
    return {s.tool_name: s.source_type.value for s in SOURCES}


def get_sidebar_config() -> list[tuple[str, str, str]]:
    """Return (type_value, icon + label, i18n_key) tuples for sidebar rendering.

    Groups sources by type for the sidebar checkboxes.
    """
    type_labels = {
        SourceType.LITERATURE: ("\U0001F4DA", "Literature & Regulatory"),
        SourceType.EXPRESSION: ("\U0001F9EC", "Expression & Protein"),
        SourceType.GENETICS: ("\U0001F9EC", "Genetics & Variants"),
        SourceType.PATHWAYS: ("\U0001F504", "Pathways"),
        SourceType.DRUGGABILITY: ("\U0001F3AF", "Druggability & Interactions"),
    }
    result = []
    for st_enum, (icon, label) in type_labels.items():
        result.append((st_enum.value, f"{icon} {label}", f"ds_{st_enum.value}"))
    return result


# ═══════════════════════════════════════════════════════════════════════════════
# Table DDL Generation
# ═══════════════════════════════════════════════════════════════════════════════

TARGET_CATALOG = CFG.catalog
TARGET_SCHEMA = CFG.schema


def generate_base_ddl(source: SourceConfig) -> str:
    """Generate CREATE TABLE DDL for a source's target table.

    All tables share a common base schema + source-specific columns.
    Base columns enable cross-source queries.
    """
    fqn = f"{TARGET_CATALOG}.{TARGET_SCHEMA}.{source.target_table}"
    return f"""CREATE TABLE IF NOT EXISTS {fqn} (
    -- Common base columns (cross-source queryable)
    record_id STRING NOT NULL COMMENT 'Unique record identifier (source-specific ID)',
    source_key STRING NOT NULL COMMENT 'Source registry key: {source.key}',
    gene_symbol STRING COMMENT 'Primary gene symbol (UPPER CASE)',
    disease STRING COMMENT 'Associated disease/condition',
    drug STRING COMMENT 'Associated drug/compound',
    title STRING COMMENT 'Record title or name',
    summary STRING COMMENT 'Brief summary or abstract',

    -- Source-specific payload (VARIANT for flexibility during early ingestion)
    payload VARIANT NOT NULL COMMENT 'Full source record as structured JSON',

    -- Metadata
    ingested_at TIMESTAMP NOT NULL,
    source_updated_at TIMESTAMP COMMENT 'Last update timestamp from source',
    api_version STRING COMMENT 'API version used for ingestion'
)
USING DELTA
COMMENT '{source.name}: {source.description}'
TBLPROPERTIES (
    'delta.enableChangeDataFeed' = 'true',
    'neuroplex.source_type' = '{source.source_type.value}',
    'neuroplex.priority' = '{source.priority.value}',
    'neuroplex.ingest_method' = '{source.ingest_method.value}'
);"""


def generate_all_ddl() -> str:
    """Generate DDL for all source tables."""
    statements = [generate_base_ddl(s) for s in SOURCES]
    return "\n\n".join(statements)


if __name__ == "__main__":
    print(f"NeuroPlex Source Registry: {len(SOURCES)} sources")
    print(f"  High priority: {len(get_sources(priority=Priority.HIGH))}")
    print(f"  Medium priority: {len(get_sources(priority=Priority.MEDIUM))}")
    print("\nBy type:")
    for st, sources in SOURCES_BY_TYPE.items():
        print(f"  {st.value}: {len(sources)} sources — {', '.join(s.name for s in sources)}")
    print(f"\nTarget tables: {len(get_target_tables())}")
    print(f"Tool schemas: {len(generate_all_tool_schemas())}")
