"""Limited AMO/Springer source checks, not compilation or proof validation."""
from pathlib import Path
import re
import json

def braced_argument(source, command):
    """Read a balanced macro argument without executing TeX."""
    match = re.search(re.escape(command) + r'\s*\{', source)
    assert match, ('missing macro', command)
    start = match.end()
    depth = 1
    for index in range(start, len(source)):
        char = source[index]
        if char in '{}' and source[index - 1] != '\\':
            depth += 1 if char == '{' else -1
            if depth == 0:
                return source[start:index]
    raise AssertionError(('unbalanced macro', command))


def check_source(source):
    source = re.sub(r'(?<!\\)%[^\n]*', '', source)
    assert r'\documentclass[pdflatex,sn-mathphys-num]{sn-jnl}' in source
    labels = re.findall(r'\\label\{([^}]+)\}', source)
    assert len(labels) == len(set(labels)), 'duplicate labels'
    refs = re.findall(r'\\(?:eqref|ref)\{([^}]+)\}', source)
    assert set(refs) <= set(labels), sorted(set(refs) - set(labels))
    items = re.findall(
        r'\\bibitem(?:\[[^]]*\])?\{([^}]+)\}(.*?)(?=\\bibitem|\\end\{thebibliography\})',
        source, re.S)
    bibs = [key for key, _ in items]
    cites = [key.strip() for group in re.findall(
        r'\\cite(?:\[[^]]*\])?\{([^}]+)\}', source) for key in group.split(',')]
    assert len(bibs) == len(set(bibs)), 'duplicate bibliography keys'
    assert set(cites) == set(bibs), 'uncited or missing bibliography entry'
    sort_keys = []
    for key, entry in items:
        author_text = entry.split(':', 1)[0]
        assert ':' in entry and author_text.strip(), ('missing author field', key)
        author_text = re.sub(r'\s+', ' ', author_text).strip().casefold()
        year = re.search(r'\((\d{4})\)', entry)
        assert year, ('missing publication year', key)
        assert r'\url{https://doi.org/' in entry, ('missing DOI link', key)
        sort_keys.append((author_text, int(year.group(1))))
    assert sort_keys == sorted(sort_keys), 'bibliography is not alphabetic/chronological'
    abstract = braced_argument(source, r'\abstract')
    assert '$' not in abstract and r'\cite' not in abstract, 'abstract must stand alone'
    abstract_words = len(abstract.split())
    assert 150 <= abstract_words <= 250, ('abstract word count', abstract_words)
    keywords = [part.strip() for part in braced_argument(source, r'\keywords').split(',')]
    assert 4 <= len(keywords) <= 6 and all(keywords), ('keywords', keywords)
    assert r'\pacs[Mathematics Subject Classification (2020)]' in source
    preprints = ('math/0607463v2', '2606.14063v1')
    body = source.split(r'\begin{thebibliography}', 1)[0]
    for identifier in preprints:
        assert 'https://arxiv.org/abs/' + identifier in body, ('missing preprint attribution', identifier)
    sections = re.findall(r'^\\section\{([^}]+)\}', source, re.M)
    assert sections == [
        'Introduction', 'System modeling', 'Theoretical derivation',
        'Numerical validation of the theoretical results',
    ], 'requested section order'
    assert r'\appendix' not in source, 'appendices must be integrated'
    theory = source.index(r'\label{sec:theory}')
    validation = source.index(r'\label{sec:validation}')
    validation_end = source.index(r'\backmatter')
    bibliography = source.index(r'\begin{thebibliography}')
    heading = r'\section*{Statements and Declarations}'
    assert source.count(heading) == 1, 'missing or duplicate declarations heading'
    declarations = source.index(heading)
    bibliography_end = source.index(r'\end{thebibliography}')
    assert theory < validation < validation_end < bibliography < bibliography_end < declarations, \
        'declarations must follow references'
    assert source.split(r'\end{thebibliography}', 1)[1].lstrip().startswith(heading), \
        'unexpected prose between references and declarations'
    for match in re.finditer(r'\\begin\{(?:theorem|lemma|corollary|proof)\}', source):
        assert theory < match.start() < validation, 'proof outside theoretical derivation'
    for label in ('sec:support-validation', 'sec:five-pilot', 'sec:matched',
                  'app:four-pilot', 'sec:five-feedback', 'sec:four-feedback',
                  'app:reproducibility'):
        assert validation < source.index(r'\label{' + label + '}') < validation_end, \
            ('validation outside numerical section', label)
    for heading in ('Funding', 'Competing interests', 'Data and code availability'):
        assert '\\bmhead{' + heading + '}' in source, ('missing declaration', heading)
    disclosure_present = r'\bmhead{Use of artificial intelligence}' in source
    assert not re.search(r'\\(?:input|include|includegraphics)\b', source), \
        'unexpected external manuscript input'
    assert r'\path{research/' not in source, 'internal research path in manuscript'
    stack = []
    count = 0
    for kind, name in re.findall(r'\\(begin|end)\{([^}]+)\}', source):
        if kind == 'begin':
            stack.append(name)
            count += 1
        else:
            assert stack and stack.pop() == name, (kind, name, stack)
    assert not stack
    return dict(labels=len(labels), bibliography_keys=len(bibs),
                attributed_preprints_in_text=len(preprints), environment_pairs=count,
                abstract_words=abstract_words, keywords=len(keywords),
                bibliography_alphabetic=True, missing_references=[],
                numbered_sections=sections, references_last=False,
                declarations_after_references=True,
                source_layout_valid=True,
                ai_use_disclosure_present=disclosure_present,
                publication_compliance_verified=False,
                check_scope="Source structure and bibliography layout only",
                compilation_verified=False)


if __name__ == '__main__':
    path = Path(__file__).with_name('representability.tex')
    print(json.dumps(check_source(path.read_text(encoding='utf-8'))))
