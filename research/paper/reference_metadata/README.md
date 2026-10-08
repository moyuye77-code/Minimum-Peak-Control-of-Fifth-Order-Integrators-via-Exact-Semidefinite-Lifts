# Bibliographic registration snapshots

The ten JSON files are unedited Crossref API responses retrieved on 2026-10-02.
Each message contains its DOI and registration metadata; the retrieval route was
`https://api.crossref.org/works/<URL-encoded DOI>`.

The offline command `python -B -S research/paper/check_references.py` compares
the current manuscript's author surnames, title, journal, volume, optional issue,
pages and volume year against these responses. It reports file hashes.

Crossref's earliest publication date is not necessarily the volume year:
Nie's rational paper uses 2012 and Fawzi's paper uses 2019 from published-print.
The Gosse title includes a funding footnote after a star in the raw response;
only that footnote is excluded from normalized title comparison.

These are metadata records, not proof checks or a guarantee that the sources
support every statement. The separate scope review is in
[the reference audit](../../../docs/REFERENCE_VERIFICATION.md).
