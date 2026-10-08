"""Compare manuscript publication fields with saved Crossref responses.

Read-only and offline. This checks bibliographic transcription, not the cited
mathematics or priority. Raw responses were obtained directly
from https://api.crossref.org/works/<DOI> on 2026-10-02 and 2026-10-03.
"""
from pathlib import Path
import hashlib
import json
import re
import unicodedata

HERE = Path(__file__).resolve().parent


def normalize(value):
    value = value.replace(r"\infty", "infinity").replace("∞", "infinity")
    value = unicodedata.normalize("NFKD", value).casefold()
    return "".join(c for c in value if c.isascii() and c.isalnum())


def require(condition, message):
    if not condition:
        raise ValueError(message)


def check(source, metadata=None):
    if metadata is None:
        metadata = {p.stem: json.loads(p.read_text(encoding="utf-8"))["message"]
                    for p in (HERE/"reference_metadata").glob("*.json")}
    entries = re.findall(
        r"\\bibitem\{([^}]+)\}(.*?)(?=\\bibitem|\\end\{thebibliography\})",
        source, re.S)
    require(len(entries) == len(metadata) == 11, "reference inventory")
    require({key for key, _ in entries} == set(metadata), "reference keys")
    rows = []
    for key, entry in entries:
        record = metadata[key]
        doi = re.search(r"\\url\{https://doi\.org/([^}]+)\}", entry)
        require(doi is not None and doi[1].casefold() == record["DOI"].casefold(),
                key+": DOI mismatch")
        author_text, content = entry.split(":", 1)
        for author in record["author"]:
            require(normalize(author["family"]) in normalize(author_text),
                    key+": author mismatch")
        title, journal_part = content.split(r"\emph{", 1)
        # Crossref's Gosse title includes a publisher-exported funding footnote.
        registered_title = record["title"][0].split("⁎", 1)[0]
        require(normalize(title) == normalize(registered_title), key+": title mismatch")
        journal = journal_part.split("}", 1)[0]
        require(normalize(journal) == normalize(record["container-title"][0]),
                key+": journal mismatch")
        volume = re.search(r"\\textbf\{([^}]+)\}", entry)
        require(volume is not None and volume[1] == record["volume"], key+": volume mismatch")
        issue = re.search(r"\\textbf\{[^}]+\}\(([^)]+)\)", entry)
        # The journal issue may be omitted, but an issue printed in the draft must match.
        if issue:
            require(issue[1].replace("--", "-") == record.get("issue"), key+": issue mismatch")
        pages = re.search(r"(\d+)--(\d+)\s*\((\d{4})\)", entry)
        require(pages is not None, key+": page/year field missing")
        require(pages[1]+"-"+pages[2] == record["page"], key+": pages mismatch")
        year_kind = "published-print" if "published-print" in record else "published"
        year = record[year_kind]["date-parts"][0][0]
        require(int(pages[3]) == year, key+": volume year mismatch")
        rows.append(dict(key=key, doi=doi[1], volume=volume[1], pages=record["page"],
                         cited_year=year, year_source=year_kind))
    return dict(references_checked=len(rows), references=rows,
                metadata_source="Crossref DOI registration responses, retrieved 2026-10-02 and 2026-10-03",
                citation_relevance_automatically_verified=False)


def main():
    source_path = HERE/"representability.tex"
    report = check(source_path.read_text(encoding="utf-8"))
    report["manuscript_sha256"] = hashlib.sha256(source_path.read_bytes()).hexdigest()
    report["metadata_sha256"] = {
        p.name: hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted((HERE/"reference_metadata").glob("*.json"))
    }
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
