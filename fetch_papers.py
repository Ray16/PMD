"""Fetch highly relevant papers about PMD/MVD enzyme engineering using
the enzyme_agent literature tools (PubMed + OpenAlex with semantic re-ranking).

Pipeline: search -> deduplicate -> relevance filter -> write CSV.
"""

import sys
import csv

sys.path.insert(0, "/nfs/lambda_stor_01/homes/rzhu/Biochemical_agent")
from enzyme_agent.tools.literature import search_pubmed, search_openalex

OUT_CSV = "/nfs/lambda_stor_01/homes/rzhu/PMD/discovered_paper/discovered_papers_curated.csv"
MIN_SEMANTIC_SCORE = 0.70

# DOIs of papers already in paper/ — exclude from output
EXCLUDE_DOIS = {
    "10.1016/j.ymben.2017.03.010",   # Kang 2017 HTP screening
    "10.1016/j.ymben.2015.12.002",   # Kang 2016 IPP-bypass
    "10.1186/s13068-022-02235-6",    # Wang 2022 P. putida isoprenoids
    "10.1128/aem.04033-14",          # Rossoni 2015 P. torridus M3K
    "10.1016/j.mec.2026.e00274",     # Kang 2026 P. putida multi-layered
}

# ---------------------------------------------------------------------------
# Highly specific queries targeting PMD/MDD and closely related enzymes
# ---------------------------------------------------------------------------
QUERIES = [
    {
        "q": "mevalonate diphosphate decarboxylase mutant substrate specificity MVAP",
        "rerank": "engineering phosphomevalonate decarboxylase PMD point mutations to improve activity toward mevalonate 5-phosphate MVAP for isopentenol production",
    },
    {
        "q": "diphosphomevalonate decarboxylase directed evolution isopentenol isoprenol",
        "rerank": "directed evolution or rational design of diphosphomevalonate decarboxylase PMD for isopentenol isoprenol biofuel production via IPP-bypass pathway",
    },
    {
        "q": "Saccharomyces cerevisiae ERG19 MVD1 mevalonate decarboxylase crystal structure",
        "rerank": "crystal structure of Saccharomyces cerevisiae mevalonate diphosphate decarboxylase ERG19 MVD1 active site substrate binding",
    },
    {
        "q": "mevalonate diphosphate decarboxylase crystal structure mechanism catalysis",
        "rerank": "crystal structure and catalytic mechanism of mevalonate diphosphate decarboxylase MDD with substrate analog ATP binding decarboxylation",
    },
    {
        "q": "Staphylococcus epidermidis mevalonate diphosphate decarboxylase structure inhibitor",
        "rerank": "Staphylococcus epidermidis mevalonate diphosphate decarboxylase crystal structure with 6-fluoromevalonate diphosphate inhibitor active site residues",
    },
    {
        "q": "Haloferax volcanii phosphomevalonate decarboxylase isopentenyl phosphate kinase",
        "rerank": "Haloferax volcanii phosphomevalonate decarboxylase that natively decarboxylates mevalonate 5-phosphate to isopentenyl phosphate archaeal mevalonate pathway",
    },
    {
        "q": "Anaerolinea thermophila mevalonate phosphate decarboxylase dual substrate",
        "rerank": "Anaerolinea thermophila enzyme with dual substrate specificity for both mevalonate 5-phosphate and mevalonate diphosphate MPD MDD interconversion",
    },
    {
        "q": "mevalonate 3-kinase Thermoplasma acidophilum structure mechanism",
        "rerank": "mevalonate-3-kinase M3K from Thermoplasma acidophilum crystal structure mechanism distinguishing kinase from decarboxylase GHMP family",
    },
    {
        "q": "mevalonate 3,5-bisphosphate decarboxylase Picrophilus torridus",
        "rerank": "mevalonate 3,5-bisphosphate decarboxylase from Picrophilus torridus crystal structure evolution of decarboxylases in mevalonate pathway",
    },
    {
        "q": "IPP-bypass mevalonate pathway isopentenol isoprenol fed-batch optimization",
        "rerank": "optimization of IPP-bypass mevalonate pathway for isopentenol isoprenol production PMD bottleneck fed-batch fermentation high titer",
    },
    {
        "q": "isopentenyl diphosphate IPP toxicity mevalonate pathway Escherichia coli",
        "rerank": "isopentenyl diphosphate IPP toxicity in E. coli mevalonate pathway growth inhibition stress response mechanism",
    },
    {
        "q": "GHMP kinase superfamily mevalonate decarboxylase phosphorylation mechanism",
        "rerank": "GHMP kinase superfamily mevalonate kinase phosphomevalonate kinase diphosphomevalonate decarboxylase mechanism phosphotransferase",
    },
    {
        "q": "isobutene 3-hydroxyisovalerate mevalonate diphosphate decarboxylase bioproduction",
        "rerank": "production of isobutene from 3-hydroxyisovalerate 3-HIV using mevalonate diphosphate decarboxylase MVD promiscuous decarboxylation activity",
    },
]

# ---------------------------------------------------------------------------
# Title keywords that mark a paper as on-topic for PMD/MDD engineering
# ---------------------------------------------------------------------------
RELEVANT_TITLE_KEYWORDS = [
    "decarboxylase", "phosphomevalonate", "diphosphomevalonate",
    "ipp-bypass", "isopentenol", "isoprenol", "isobutene",
    "mevalonate kinase", "mevalonate-3-kinase", "mevalonate 3-phosphate",
    "mevalonate 3,5-bisphosphate", "ghmp kinase",
    "isoprenoid precursor toxicity", "ipp toxicity",
    "isopentenyl phosphate kinase",
]


def run_searches() -> dict:
    """Run all queries against PubMed + OpenAlex, dedup by DOI."""
    all_papers = {}  # doi_lower -> paper dict

    for i, entry in enumerate(QUERIES):
        q, rr = entry["q"], entry["rerank"]
        print(f"\n[{i+1}/{len(QUERIES)}] {q[:70]}...")

        for source, search_fn, kwargs in [
            ("pubmed", search_pubmed, dict(query=q, max_results=8, rerank_query=rr)),
            ("openalex", search_openalex, dict(query=q, max_results=8,
                                                sort_by_citations=True, rerank_query=rr)),
        ]:
            try:
                result = search_fn(**kwargs)
                for p in result.get("data", []):
                    doi = (p.get("doi") or "").strip()
                    doi_key = doi.lower().rstrip("/")
                    if doi and doi_key not in all_papers:
                        p["_source"] = source
                        all_papers[doi_key] = p
            except Exception as e:
                print(f"  {source} error: {e}")

    print(f"\nTotal unique papers (by DOI): {len(all_papers)}")
    return all_papers


def filter_papers(papers: dict) -> list[dict]:
    """Keep only papers with high semantic score AND on-topic title keywords.
    Excludes papers already in paper/ directory (by DOI)."""
    exclude_norm = {d.lower().rstrip("/") for d in EXCLUDE_DOIS}
    filtered = []
    for doi_key, p in papers.items():
        if doi_key in exclude_norm:
            continue
        score = p.get("semantic_score", 0)
        if score < MIN_SEMANTIC_SCORE:
            continue
        title_lower = (p.get("title") or "").lower()
        if not any(kw in title_lower for kw in RELEVANT_TITLE_KEYWORDS):
            continue
        filtered.append(p)

    filtered.sort(key=lambda x: x.get("semantic_score", 0), reverse=True)
    print(f"After relevance filter (excl. {len(exclude_norm)} existing): {len(filtered)} papers")
    return filtered


def write_csv(papers: list[dict]) -> None:
    """Write final curated CSV."""
    with open(OUT_CSV, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "Title", "Authors", "Journal", "Year", "DOI",
            "Semantic_Score", "Citations", "Abstract_Snippet", "Source",
        ])
        for p in papers:
            authors = p.get("authors", [])
            if isinstance(authors, list):
                authors = "; ".join(str(a) for a in authors[:5])
            abstract = (p.get("abstract") or "")[:300].replace("\n", " ")
            w.writerow([
                p.get("title", ""),
                authors,
                p.get("journal", ""),
                p.get("year", ""),
                p.get("doi", ""),
                p.get("semantic_score", ""),
                p.get("cited_by_count", ""),
                abstract,
                p.get("_source", ""),
            ])
    print(f"Wrote {len(papers)} papers to {OUT_CSV}")


if __name__ == "__main__":
    raw = run_searches()
    curated = filter_papers(raw)
    write_csv(curated)
