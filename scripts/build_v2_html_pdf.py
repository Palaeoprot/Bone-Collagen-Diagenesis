import os, subprocess, markdown

src_md = r'c:\Users\matth\Documents\GitHub\Palaeoprot-Publications\projects\Ea Collagen Hydrolysis\manuscript_thermal_age_collagen_hydrolysis_v2.0_revision_proposal.md'
out_html = r'c:\Users\matth\Documents\GitHub\Palaeoprot-Publications\projects\Ea Collagen Hydrolysis\manuscript_thermal_age_collagen_hydrolysis_v2.0_revision_proposal.html'
out_pdf = r'c:\Users\matth\Documents\GitHub\Palaeoprot-Publications\projects\Ea Collagen Hydrolysis\manuscript_thermal_age_collagen_hydrolysis_v2.0_revision_proposal.pdf'

with open(src_md, 'r', encoding='utf-8') as f:
    text = f.read()

html_body = markdown.markdown(text, extensions=['tables', 'fenced_code'])

header = """<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>Deconvolving the Radiocarbon Wall (Version 2.0 Revision Proposal)</title>
<style>
body { font-family: 'Segoe UI', Arial, sans-serif; line-height: 1.6; max-width: 900px; margin: 40px auto; padding: 0 20px; color: #222; }
h1, h2, h3 { color: #112233; }
s { color: #c0392b; text-decoration: line-through; opacity: 0.75; }
table { border-collapse: collapse; width: 100%; margin: 20px 0; }
th, td { border: 1px solid #ccc; padding: 8px 12px; text-align: left; }
th { background-color: #f2f4f8; }
pre, code { background-color: #f8f9fa; font-family: Consolas, monospace; }
pre { padding: 12px; border-radius: 4px; overflow-x: auto; }
</style>
</head>
<body>
"""

footer = """
</body>
</html>
"""

with open(out_html, 'w', encoding='utf-8') as f:
    f.write(header + html_body + footer)

print("HTML built successfully:", out_html)

edge_path = r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe"
if os.path.exists(edge_path):
    subprocess.run([edge_path, "--headless", "--disable-gpu", f"--print-to-pdf={out_pdf}", out_html], check=True)
    print("PDF built successfully:", out_pdf)
