#!/usr/bin/env python3
"""Convert comprehensive_guide.md to PDF with embedded images."""
import os
import markdown
from weasyprint import HTML

ANALYSIS_DIR = os.path.dirname(os.path.abspath(__file__))
MD_FILE = os.path.join(ANALYSIS_DIR, "comprehensive_guide.md")
PDF_FILE = os.path.join(ANALYSIS_DIR, "comprehensive_guide.pdf")

with open(MD_FILE) as f:
    md_text = f.read()

html_body = markdown.markdown(md_text, extensions=["tables", "fenced_code"])

# Resolve relative image paths to absolute file:// URIs
for img_file in os.listdir(ANALYSIS_DIR):
    if img_file.endswith(".png"):
        abs_path = os.path.join(ANALYSIS_DIR, img_file)
        html_body = html_body.replace(f'src="{img_file}"', f'src="file://{abs_path}"')

html_doc = f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<style>
  body {{
    font-family: "Helvetica Neue", Helvetica, Arial, sans-serif;
    font-size: 11pt;
    line-height: 1.5;
    max-width: 7.5in;
    margin: 0.75in auto;
    color: #222;
  }}
  h1 {{ font-size: 20pt; margin-top: 0; }}
  h2 {{ font-size: 15pt; margin-top: 1.5em; border-bottom: 1px solid #ccc; padding-bottom: 4px; }}
  h3 {{ font-size: 12pt; margin-top: 1.2em; }}
  table {{
    border-collapse: collapse;
    width: 100%;
    margin: 1em 0;
    font-size: 10pt;
  }}
  th, td {{
    border: 1px solid #bbb;
    padding: 5px 8px;
    text-align: left;
  }}
  th {{ background-color: #f0f0f0; font-weight: bold; }}
  code {{
    background-color: #f5f5f5;
    padding: 1px 4px;
    border-radius: 3px;
    font-size: 10pt;
  }}
  pre {{
    background-color: #f5f5f5;
    padding: 10px;
    border-radius: 4px;
    overflow-x: auto;
    font-size: 9pt;
    line-height: 1.4;
  }}
  pre code {{ background: none; padding: 0; }}
  img {{
    max-width: 100%;
    height: auto;
    display: block;
    margin: 1em auto;
  }}
  hr {{ border: none; border-top: 1px solid #ddd; margin: 1.5em 0; }}
  strong {{ color: #111; }}
  @page {{
    size: letter;
    margin: 0.75in;
  }}
</style>
</head>
<body>
{html_body}
</body>
</html>"""

HTML(string=html_doc).write_pdf(PDF_FILE)
print(f"PDF written to {PDF_FILE}")
