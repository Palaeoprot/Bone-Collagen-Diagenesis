import os, sys, re, email
from email import policy
from bs4 import BeautifulSoup
from markdownify import markdownify as md
import datetime

def clean_one_mhtml(mhtml_path, out_md_path, title_override=None, author_override=None, year_override=None):
    with open(mhtml_path, "rb") as fp:
        msg = email.message_from_binary_file(fp, policy=policy.default)

    html_content = ""
    for part in msg.walk():
        if "html" in part.get_content_type():
            payload = part.get_payload(decode=True)
            charset = part.get_content_charset() or "utf-8"
            html_content = payload.decode(charset, errors="replace")
            break

    if not html_content:
        print(f"Warning: No HTML in {mhtml_path}")
        return False

    soup = BeautifulSoup(html_content, "html.parser")
    
    # Title resolution
    doc_title = title_override
    if not doc_title:
        title_meta = soup.find('meta', {'name': re.compile(r'citation_title|dc\.title', re.I)})
        if title_meta and 'content' in title_meta.attrs:
            doc_title = title_meta['content'].strip()
        elif soup.title and soup.title.string:
            doc_title = soup.title.string.strip()
        else:
            doc_title = os.path.splitext(os.path.basename(mhtml_path))[0]
            
    doc_title = re.sub(r"[\r\n\t]+", " ", doc_title)
    doc_title = re.sub(r"\s*\|\s*(PLOS ONE|Nature.*|PLOS One|PMC|ScienceDirect|Springer.*|Cambridge Core|Wiley.*|ACS Publications|Antiquity|Radiocarbon|Communications.*).*$", "", doc_title, flags=re.I)
    doc_title = re.sub(r"\s*-\s*(PMC|Wiley Online Library|ScienceDirect|ACS Publications|Springer Nature Link|Antiquity|Radiocarbon).*$", "", doc_title, flags=re.I)
    doc_title = re.sub(r"\s*_\s*(Nature|PLOS One|Journal.*|ACS.*|Cambridge Core|Springer Nature Link|Physical Chemistry.*|Soft Matter.*|Communications.*).*$", "", doc_title, flags=re.I)
    doc_title = re.sub(r"^Full article:\s*", "", doc_title, flags=re.I).strip()

    # Target content area
    nature_main = soup.find("main", class_=re.compile(r"c-article-main-column"))
    pmc_main = soup.find("main", id="main-content")
    plos_main = soup.find("main", id="main-content")
    wiley_main = soup.find("section", class_=re.compile(r"article-section__content"))
    cambridge_main = soup.find("div", class_=re.compile(r"article-content|core-content"))

    if nature_main:
        main_elem = nature_main
        for bad in main_elem.find_all(class_=re.compile(r"(c-article-recommendations|c-context-bar|c-article-share-box|c-site-messages)", re.I)):
            bad.decompose()
    elif pmc_main:
        main_elem = pmc_main
        for bad in main_elem.find_all(class_=re.compile(r"(usa-banner|pmc-footer|share-tool|social-share)", re.I)):
            bad.decompose()
    elif plos_main:
        main_elem = plos_main
        for tab in main_elem.find_all(class_=re.compile(r"(tab-list|article-tabs|nav-tabs)", re.I)):
            tab.decompose()
    elif wiley_main:
        main_elem = wiley_main
    elif cambridge_main:
        main_elem = cambridge_main
    else:
        main_elem = (
            soup.find("article", class_=re.compile(r"article|content|paper", re.I))
            or soup.find("main")
            or soup.find("article")
            or soup.body
        )

    if not main_elem:
        print(f"Warning: No main element found in {mhtml_path}")
        return False

    # Strip navigation, scripts, ads
    for tag in main_elem(["script", "style", "nav", "footer", "noscript", "svg", "form", "iframe", "header"]):
        tag.decompose()

    for tag in main_elem.find_all(class_=re.compile(r"(sidebar|navigation|menu|advert|banner|cookie|share|social|metric)", re.I)):
        tag.decompose()

    # Convert to Markdown
    md_content = md(str(main_elem), heading_style="ATX", strip=["img"])
    md_content = re.sub(r"\n{3,}", "\n\n", md_content).strip()

    # Apply Document Timestamping Directive
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S (+02:00)")
    author_str = f"**Authors:** {author_override}\n\n" if author_override else ""
    year_str = f"**Year:** {year_override}\n\n" if year_override else ""
    header = f"# {doc_title}\n**Date & Time:** {now}\n\n{author_str}{year_str}"

    final_md = header + md_content

    os.makedirs(os.path.dirname(os.path.abspath(out_md_path)), exist_ok=True)
    with open(out_md_path, "w", encoding="utf-8") as out_fp:
        out_fp.write(final_md)

    return True

if __name__ == "__main__":
    import json
    proj_dir = r"c:\Users\matth\Documents\GitHub\Palaeoprot-Publications\projects\Ea Collagen Hydrolysis"
    with open(os.path.join(proj_dir, "papers_metadata.json"), "r", encoding="utf-8") as f:
        papers = json.load(f)

    success_count = 0
    for p in papers:
        if p["ext"] == ".mhtml":
            fn = p["filename"]
            folder = p.get("folder", "")
            src_mhtml = os.path.join(proj_dir, folder, fn)
            # Create cleaned markdown filename
            auth = (p.get("author") or "Unknown").split(",")[0].split(" ")[-1]
            yr = p.get("year") or "ND"
            safe_title = re.sub(r'[^a-zA-Z0-9_\- ]', '', p.get("title") or fn)[:40].strip()
            out_name = f"{auth}_{yr}_{safe_title}.md"
            out_path = os.path.join(proj_dir, folder, "cleaned_markdown", out_name)
            
            res = clean_one_mhtml(
                src_mhtml,
                out_path,
                title_override=p.get("title"),
                author_override=p.get("author"),
                year_override=p.get("year")
            )
            if res:
                success_count += 1
                p["cleaned_md"] = os.path.relpath(out_path, proj_dir)

    print(f"Successfully cleaned and converted {success_count} mhtml files to markdown.")
    with open(os.path.join(proj_dir, "papers_metadata.json"), "w", encoding="utf-8") as f:
        json.dump(papers, f, indent=2, ensure_ascii=False)
