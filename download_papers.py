"""Download discovered papers for PMD engineering project.

Attempts multiple open-access sources: PMC, Unpaywall, direct publisher links.
Validates that downloaded files are real PDFs (>50KB and start with %PDF).
"""

import csv
import os
import subprocess
import time

CURATED_CSV = "/nfs/lambda_stor_01/homes/rzhu/PMD/discovered_paper/discovered_papers_curated.csv"
OUT_DIR = "/nfs/lambda_stor_01/homes/rzhu/PMD/discovered_paper"
MIN_PDF_SIZE = 50_000  # 50KB minimum to be considered a real PDF

# DOIs already present in paper/ (the 5 original papers) — skip these
ALREADY_HAVE = {
    "10.1016/j.ymben.2017.03.010",   # Kang 2017 HTP screening
    "10.1016/j.ymben.2015.12.002",   # Kang 2016 IPP-bypass
    "10.1186/s13068-022-02235-6",    # Wang 2022 P. putida isoprenoids
    "10.1128/aem.04033-14",          # Rossoni 2015 P. torridus M3K
    "10.1016/j.mec.2026.e00274",     # Kang 2026 P. putida multi-layered
}

# Map DOI -> desired filename
DOI_TO_FILENAME = {
    "10.1128/JB.01230-13": "Vannice_2014_JBacteriol_Hv_PMD_IPK.pdf",
    "10.1021/bi300591x": "Barta_2012_Biochemistry_SeMDD_nucleotide_binding.pdf",
    "10.1016/j.jbc.2022.102111": "Vinokur_2022_JBC_MVA35BD_structure.pdf",
    "10.1038/s41467-020-17733-0": "Chen_2020_NatCommun_visualizing_MDD_mechanism.pdf",
    "10.1007/s00894-015-2873-0": "Addo_2016_JMolModel_SaMDD_inhibitor_docking.pdf",
    "10.1074/jbc.M110787200": "Fu_2002_JBC_MjMVK_GHMP_structure.pdf",
    "10.1074/jbc.m114.562686": "Azami_2014_JBC_MVA3P_Thermoplasma.pdf",
    "10.1074/jbc.m111.242016": "Barta_2011_JBC_SeMDD_inhibitor_structures.pdf",
    "10.2210/pdb6n0x/pdb": "Thomas_2018_PDB_AtMPD_MVAP_structure.pdf",
    "10.1021/acschembio.9b00322": "Lubin_2019_ACSChemBiol_MPD_engineering.pdf",
    "10.1016/j.ymben.2019.09.003": "Kang_2019_MetEng_IPPbypass_fedbatch.pdf",
    "10.1016/j.abb.2014.12.002": "Skaff_2015_ABB_MDD_eriochrome_inhibition.pdf",
    "10.1074/jbc.m117.802223": "Chen_2017_JBC_EfMDD_ATP_binding.pdf",
    "10.1039/c9ob02254f": "McClory_2019_OrgBiomolChem_MDD_QMMM.pdf",
    "10.1074/jbc.m116.752535": "Motoyama_2017_JBC_single_AA_decarboxylase_to_kinase.pdf",
    "10.1002/pro.2607": "Vinokur_2015_ProteinSci_M3K_structure.pdf",
    "10.1016/j.ymben.2024.02.004": "Banerjee_2024_MetEng_Pputida_isoprenol.pdf",
}


def normalize_doi(doi):
    return doi.lower().strip().rstrip("/")


def is_valid_pdf(filepath):
    """Check file exists, is large enough, and starts with %PDF."""
    if not os.path.exists(filepath):
        return False
    size = os.path.getsize(filepath)
    if size < MIN_PDF_SIZE:
        return False
    try:
        with open(filepath, "rb") as f:
            header = f.read(5)
        return header == b"%PDF-"
    except Exception:
        return False


def try_download(url, filepath, timeout=30):
    """Download URL to filepath using curl. Returns True if valid PDF."""
    try:
        subprocess.run(
            ["curl", "-L", "-f", "-s", "-o", filepath,
             "-H", "User-Agent: Mozilla/5.0 (research; mailto:enzyme-agent@example.edu)",
             "--max-time", str(timeout), url],
            check=True, capture_output=True
        )
        return is_valid_pdf(filepath)
    except subprocess.CalledProcessError:
        if os.path.exists(filepath):
            os.remove(filepath)
        return False


def download_paper(doi, filename):
    """Try multiple sources to download a paper by DOI."""
    filepath = os.path.join(OUT_DIR, filename)

    # Skip if already have a valid PDF
    if is_valid_pdf(filepath):
        print(f"  SKIP (already valid): {filename}")
        return True

    # Remove any broken file
    if os.path.exists(filepath):
        os.remove(filepath)

    doi_clean = doi.strip()

    # --- Source 1: Unpaywall (finds legal open-access PDFs) ---
    try:
        import json, urllib.request
        unpaywall_url = f"https://api.unpaywall.org/v2/{doi_clean}?email=enzyme-agent@example.edu"
        req = urllib.request.Request(unpaywall_url)
        with urllib.request.urlopen(req, timeout=10) as resp:
            data = json.loads(resp.read())
        # Try best OA location
        best = data.get("best_oa_location") or {}
        pdf_url = best.get("url_for_pdf") or best.get("url")
        if pdf_url and try_download(pdf_url, filepath):
            print(f"  OK (unpaywall): {filename}")
            return True
        # Try all OA locations
        for loc in data.get("oa_locations", []):
            pdf_url = loc.get("url_for_pdf") or loc.get("url")
            if pdf_url and try_download(pdf_url, filepath):
                print(f"  OK (unpaywall-alt): {filename}")
                return True
    except Exception:
        pass

    # --- Source 2: PMC (for NIH-funded papers) ---
    try:
        import json, urllib.request
        pmc_url = f"https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?ids={doi_clean}&format=json"
        req = urllib.request.Request(pmc_url)
        with urllib.request.urlopen(req, timeout=10) as resp:
            data = json.loads(resp.read())
        for rec in data.get("records", []):
            pmcid = rec.get("pmcid")
            if pmcid:
                pdf_link = f"https://pmc.ncbi.nlm.nih.gov/articles/{pmcid}/pdf/"
                if try_download(pdf_link, filepath):
                    print(f"  OK (PMC {pmcid}): {filename}")
                    return True
    except Exception:
        pass

    # --- Source 3: Direct publisher PDF (works for some open-access journals) ---
    direct_urls = []
    if "nature.com" in doi_clean or "s41467" in doi_clean:
        direct_urls.append(f"https://www.nature.com/articles/{doi_clean.split('/')[-1]}.pdf")
    if "10.1074/jbc" in doi_clean:
        direct_urls.append(f"https://www.jbc.org/article/{doi_clean}/fulltext")
    if "10.1002/pro" in doi_clean:
        direct_urls.append(f"https://onlinelibrary.wiley.com/doi/pdfdirect/{doi_clean}")

    for url in direct_urls:
        if try_download(url, filepath):
            print(f"  OK (direct): {filename}")
            return True

    # --- Source 4: EuropePMC ---
    try:
        import json, urllib.request
        epmc_url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=DOI:{doi_clean}&format=json"
        req = urllib.request.Request(epmc_url)
        with urllib.request.urlopen(req, timeout=10) as resp:
            data = json.loads(resp.read())
        for result in data.get("resultList", {}).get("result", []):
            pmcid = result.get("pmcid")
            if pmcid:
                pdf_link = f"https://europepmc.org/backend/ptpmcrender.fcgi?accid={pmcid}&blobtype=pdf"
                if try_download(pdf_link, filepath):
                    print(f"  OK (EuropePMC {pmcid}): {filename}")
                    return True
    except Exception:
        pass

    print(f"  FAIL: {filename} (DOI: {doi_clean})")
    return False


def main():
    # Read curated CSV
    with open(CURATED_CSV) as f:
        papers = list(csv.DictReader(f))

    print(f"Curated papers: {len(papers)}")
    print(f"Already have (in paper/): {len(ALREADY_HAVE)}\n")

    success, fail = 0, 0
    for p in papers:
        doi = p["DOI"].strip()
        doi_norm = normalize_doi(doi)

        # Skip papers already in paper/ folder
        if any(normalize_doi(d) == doi_norm for d in ALREADY_HAVE):
            print(f"SKIP (in paper/): {p['Title'][:70]}")
            continue

        # Get filename
        filename = None
        for map_doi, map_name in DOI_TO_FILENAME.items():
            if normalize_doi(map_doi) == doi_norm:
                filename = map_name
                break
        if not filename:
            # Auto-generate filename from first author and year
            author = (p.get("Authors") or "Unknown").split(";")[0].split()[-1].strip()
            year = p.get("Year", "")
            filename = f"{author}_{year}.pdf"

        print(f"\nDownloading: {p['Title'][:80]}")
        if download_paper(doi, filename):
            success += 1
        else:
            fail += 1
        time.sleep(1)  # rate limit

    print(f"\n{'='*60}")
    print(f"Downloaded: {success}, Failed: {fail}")

    # Final validation: check all PDFs in directory
    print(f"\nValidating all PDFs in {OUT_DIR}:")
    for f in sorted(os.listdir(OUT_DIR)):
        if f.endswith(".pdf"):
            fp = os.path.join(OUT_DIR, f)
            size = os.path.getsize(fp)
            valid = is_valid_pdf(fp)
            status = "OK" if valid else "BROKEN"
            print(f"  [{status}] {size:>10,} bytes  {f}")


if __name__ == "__main__":
    main()
