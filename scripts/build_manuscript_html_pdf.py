import os
import re
import subprocess
import markdown

def math_protect(text):
    math_blocks = []
    def save_block(match):
        idx = len(math_blocks)
        math_blocks.append(match.group(0))
        return f"MATHBLOCKPROTECTED{idx}ENDMATH"
    
    # Protect $$ ... $$
    text = re.sub(r'\$\$(.*?)\$\$', save_block, text, flags=re.DOTALL)
    # Protect $ ... $
    text = re.sub(r'(?<!\\)\$(.*?)(?<!\\)\$', save_block, text, flags=re.DOTALL)
    return text, math_blocks

def math_restore(html, math_blocks):
    for idx, block in enumerate(math_blocks):
        html = html.replace(f"MATHBLOCKPROTECTED{idx}ENDMATH", block)
    return html

def build_manuscript():
    src_md = r"c:\Users\matth\Documents\GitHub\Palaeoprot-Publications\projects\Ea Collagen Hydrolysis\manuscript_thermal_age_collagen_hydrolysis.md"
    dst_dir_proj = r"c:\Users\matth\Documents\GitHub\Palaeoprot-Publications\projects\Ea Collagen Hydrolysis"
    dst_dir_manu = r"C:\Users\matth\Documents\GitHub\manuscripts\07_Thermal_Age_and_Collagen_Hydrolysis_Ea"
    
    with open(src_md, "r", encoding="utf-8") as f:
        md_content = f.read()

    protected_md, math_blocks = math_protect(md_content)
    
    html_body = markdown.markdown(
        protected_md,
        extensions=['tables', 'fenced_code', 'toc', 'sane_lists', 'attr_list']
    )
    
    html_body = math_restore(html_body, math_blocks)
    
    # Process images and captions
    # Convert <p><img ... alt="Caption" /></p><p><em><strong>Figure X:</strong> ...</em></p> into structured figure cards
    html_body = re.sub(
        r'<p><img\s+alt="([^"]*)"\s+src="([^"]*)"\s*/?></p>\s*<p><em><strong>(Figure\s+\d+:?)</strong>(.*?)</em></p>',
        r'<div class="figure-container"><img src="\2" alt="\1" class="figure-img" /><p class="figure-caption"><strong>\3</strong>\4</p></div>',
        html_body,
        flags=re.DOTALL
    )

    template = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Deconvolving the Radiocarbon Wall: Empirical Activation Energy of Bone Collagen Hydrolysis and Deep-Time Thermal Age</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link href="https://fonts.googleapis.com/css2?family=Cinzel:wght@600;700;800&family=JetBrains+Mono:wght@400;500;600&family=Latin+Modern+Roman:ital,wght@0,400;0,700;1,400&family=Newsreader:ital,opsz,wght@0,6..72,400;0,6..72,600;0,6..72,700;1,6..72,400;1,6..72,600&display=swap" rel="stylesheet">
<script>
window.MathJax = {{
  tex: {{
    inlineMath: [['$', '$'], ['\\\\(', '\\\\)']],
    displayMath: [['$$', '$$'], ['\\\\[', '\\\\]']],
    processEscapes: true
  }},
  options: {{
    skipHtmlTags: ['script', 'noscript', 'style', 'textarea', 'pre', 'code']
  }}
}};
</script>
<script id="MathJax-script" async src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script>
<style>
@page {{
  size: A4;
  margin: 20mm 18mm 22mm 18mm;
  @bottom-right {{
    content: counter(page);
    font-family: 'Latin Modern Roman', 'Computer Modern', 'Newsreader', Georgia, serif;
    font-size: 9pt;
    color: #64748B;
  }}
}}

*, *::before, *::after {{
  box-sizing: border-box;
}}

body {{
  font-family: 'Latin Modern Roman', 'Newsreader', 'Computer Modern', 'Times New Roman', Georgia, serif;
  font-size: 10pt;
  line-height: 1.55;
  color: #111827;
  background-color: #FFFFFF;
  margin: 0 auto;
  max-width: 190mm;
  padding: 10mm 15mm;
  text-rendering: optimizeLegibility;
}}

h1 {{
  font-family: 'Latin Modern Roman', 'Newsreader', 'Times New Roman', serif;
  font-size: 19pt;
  font-weight: 700;
  line-height: 1.25;
  text-align: center;
  color: #0F172A;
  margin-top: 15pt;
  margin-bottom: 12pt;
  letter-spacing: -0.01em;
}}

.timestamp {{
  text-align: center;
  font-size: 8.5pt;
  color: #64748B;
  margin-bottom: 14pt;
}}

.authors {{
  text-align: center;
  font-size: 10.5pt;
  font-weight: 600;
  color: #1E293B;
  margin-bottom: 4pt;
}}

.affiliations {{
  text-align: center;
  font-size: 8.5pt;
  font-style: italic;
  color: #475569;
  margin-bottom: 16pt;
  line-height: 1.4;
}}

hr {{
  border: none;
  border-top: 1pt solid #E2E8F0;
  margin: 18pt 0;
}}

h2 {{
  font-family: 'Latin Modern Roman', 'Newsreader', 'Times New Roman', serif;
  font-size: 13pt;
  font-weight: 700;
  color: #0F172A;
  border-bottom: 1.2pt solid #CBD5E1;
  padding-bottom: 3.5pt;
  margin-top: 22pt;
  margin-bottom: 10pt;
  page-break-after: avoid;
}}

h3 {{
  font-family: 'Latin Modern Roman', 'Newsreader', 'Times New Roman', serif;
  font-size: 11pt;
  font-weight: 700;
  color: #1E293B;
  margin-top: 14pt;
  margin-bottom: 6pt;
  page-break-after: avoid;
}}

h4 {{
  font-size: 10pt;
  font-weight: 700;
  color: #334155;
  margin-top: 10pt;
  margin-bottom: 4pt;
  page-break-after: avoid;
}}

p {{
  margin-top: 0;
  margin-bottom: 8pt;
  text-align: justify;
  hyphens: auto;
}}

blockquote {{
  margin: 10pt 0;
  padding: 8pt 14pt;
  background-color: #F8FAFC;
  border-left: 3pt solid #3B82F6;
  color: #334155;
  font-style: italic;
}}

table {{
  width: 100%;
  border-collapse: collapse;
  margin: 14pt 0;
  font-size: 8.5pt;
  line-height: 1.35;
  page-break-inside: avoid;
}}

th, td {{
  border: 0.5pt solid #CBD5E1;
  padding: 5pt 7pt;
  text-align: left;
}}

th {{
  background-color: #F1F5F9;
  font-weight: 700;
  color: #0F172A;
}}

tr:nth-child(even) {{
  background-color: #F8FAFC;
}}

pre, code {{
  font-family: 'JetBrains Mono', 'Courier New', monospace;
  font-size: 8pt;
}}

pre {{
  background-color: #0F172A;
  color: #F8FAFC;
  padding: 10pt 12pt;
  border-radius: 4pt;
  overflow-x: auto;
  margin: 12pt 0;
  line-height: 1.35;
  page-break-inside: avoid;
}}

code {{
  background-color: #F1F5F9;
  color: #0F172A;
  padding: 1.5pt 3pt;
  border-radius: 2pt;
}}

pre code {{
  background-color: transparent;
  color: inherit;
  padding: 0;
}}

.figure-container {{
  margin: 16pt 0;
  text-align: center;
  page-break-inside: avoid;
}}

.figure-img {{
  max-width: 100%;
  height: auto;
  border-radius: 3pt;
  box-shadow: 0 1px 3px rgba(0,0,0,0.1);
  margin-bottom: 6pt;
}}

.figure-caption {{
  font-size: 8.5pt;
  line-height: 1.4;
  color: #334155;
  text-align: justify;
  margin: 4pt 10pt 12pt 10pt;
  padding: 4pt 8pt;
  background-color: #F8FAFC;
  border-left: 2pt solid #64748B;
}}

ul, ol {{
  margin-top: 0;
  margin-bottom: 8pt;
  padding-left: 18pt;
}}

li {{
  margin-bottom: 3pt;
}}

.abstract-box {{
  background-color: #F8FAFC;
  border: 1pt solid #E2E8F0;
  border-radius: 4pt;
  padding: 12pt 16pt;
  margin: 14pt 0;
}}

.abstract-title {{
  font-weight: 700;
  font-size: 10.5pt;
  color: #0F172A;
  margin-bottom: 6pt;
  text-transform: uppercase;
  letter-spacing: 0.05em;
}}

@media print {{
  body {{
    padding: 0;
    max-width: 100%;
  }}
  .figure-img {{
    box-shadow: none;
  }}
}}
</style>
</head>
<body>
{html_body}
</body>
</html>
"""

    # Format Abstract box nicely in HTML
    template = template.replace(
        "<h3>Abstract</h3>",
        "<div class=\"abstract-box\"><div class=\"abstract-title\">Abstract</div>"
    )
    # Close abstract box before keywords or hr
    template = template.replace(
        "<p><strong>Keywords:</strong>",
        "</div>\n<p><strong>Keywords:</strong>"
    )

    out_html_proj = os.path.join(dst_dir_proj, "manuscript_thermal_age_collagen_hydrolysis.html")
    out_html_manu = os.path.join(dst_dir_manu, "manuscript_thermal_age_collagen_hydrolysis.html")
    
    with open(out_html_proj, "w", encoding="utf-8") as f:
        f.write(template)
    with open(out_html_manu, "w", encoding="utf-8") as f:
        f.write(template)
        
    print(f"HTML generated successfully: {out_html_proj}")
    
    # Generate PDF using Microsoft Edge / Chrome headless
    out_pdf_proj = os.path.join(dst_dir_proj, "manuscript_thermal_age_collagen_hydrolysis.pdf")
    out_pdf_manu = os.path.join(dst_dir_manu, "manuscript_thermal_age_collagen_hydrolysis.pdf")
    
    edge_exe = r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe"
    chrome_exe = r"C:\Program Files\Google\Chrome\Application\chrome.exe"
    browser_exe = edge_exe if os.path.exists(edge_exe) else chrome_exe
    
    print(f"Compiling PDF via {browser_exe}...")

    # The legacy --headless flag silently no-ops (exit 0, no file written) when
    # another browser instance is already running, so use --headless=new with an
    # isolated profile, then VERIFY the file actually changed on disk.
    import shutil, tempfile, time
    before = os.path.getmtime(out_pdf_proj) if os.path.exists(out_pdf_proj) else 0.0
    if os.path.exists(out_pdf_proj):
        try:
            os.remove(out_pdf_proj)
        except OSError:
            pass

    profile = tempfile.mkdtemp(prefix="pdfbuild_")
    file_url = "file:///" + out_html_proj.replace("\\", "/").lstrip("/")
    cmd = [
        browser_exe,
        "--headless=new",
        "--disable-gpu",
        "--no-sandbox",
        "--no-first-run",
        "--no-pdf-header-footer",
        "--run-all-compositor-stages-before-draw",
        "--virtual-time-budget=20000",
        f"--user-data-dir={profile}",
        f"--print-to-pdf={out_pdf_proj}",
        file_url,
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    for _ in range(20):                      # the write can lag process exit
        if os.path.exists(out_pdf_proj) and os.path.getsize(out_pdf_proj) > 1000:
            break
        time.sleep(0.5)
    shutil.rmtree(profile, ignore_errors=True)

    if not os.path.exists(out_pdf_proj) or os.path.getsize(out_pdf_proj) < 1000:
        raise SystemExit(
            "PDF generation FAILED - no output written. "
            f"exit={proc.returncode} stderr={proc.stderr[-600:]}"
        )
    after = os.path.getmtime(out_pdf_proj)
    if after <= before:
        raise SystemExit("PDF generation FAILED - existing file was not replaced.")

    shutil.copy2(out_pdf_proj, out_pdf_manu)
    size_mb = os.path.getsize(out_pdf_proj) / 1e6
    print(f"PDF compiled successfully ({size_mb:.1f} MB): {out_pdf_proj}")
    print(f"PDF mirrored successfully: {out_pdf_manu}")

if __name__ == "__main__":
    build_manuscript()
