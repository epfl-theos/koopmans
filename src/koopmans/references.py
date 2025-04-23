"""References module for koopmans."""

from pathlib import Path

from pybtex.database.input import bibtex

parser = bibtex.Parser()
bib_file = Path(__file__).parent / 'references.bib'
bib_data = parser.parse_file(bib_file)
