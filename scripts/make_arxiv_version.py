"""The arXiv version of the paper (docs/arxiv), generated from the JCAP version (docs/paper).

Edit only docs/paper/ngpbh_paper.tex; this script rewrites it into the PRIMEarxiv template:
the template's preamble, author block, abstract and keywords, capitalised cross-references, table
captions above the tables, figures copied from docs/paper/figures to docs/arxiv/media, and
references.bib adapted to unsrt (arXiv numbers in the note field, arXiv-only papers as misc).
Then, in docs/arxiv: pdflatex ngpbh_arxiv && bibtex ngpbh_arxiv && pdflatex ngpbh_arxiv (twice).

Usage: python scripts/make_arxiv_version.py
"""

import re
import shutil
from pathlib import Path

for fig in Path("docs/paper/figures").glob("*.pdf"):
    shutil.copy2(fig, Path("docs/arxiv/media") / fig.name)

SRC = Path("docs/paper/ngpbh_paper.tex").read_text()


def braced(s, start):
    """Content of the brace group that opens at s[start] == '{'; returns (content, end_index)."""
    assert s[start] == "{"
    depth = 0
    for i in range(start, len(s)):
        if s[i] == "{":
            depth += 1
        elif s[i] == "}":
            depth -= 1
            if depth == 0:
                return s[start + 1:i], i + 1
    raise ValueError("unbalanced braces")


abstract, _ = braced(SRC, SRC.index("\\abstract{") + len("\\abstract"))
macros = re.search(r"(\\newcommand\{\\Cg\}.*?\\newcommand\{\\Pz\}\{[^\n]*\})", SRC, re.DOTALL).group(1)
body = SRC.split("\\maketitle\n\\flushbottom\n", 1)[1].split("\\end{document}")[0]

# JCAP commands -> standard article / template
body = body.replace("\\acknowledgments\n", "\\section*{Acknowledgments}\n")
body = body.replace("\\bibliographystyle{JHEP}\n\\bibliography{refs}",
                    "%Bibliography\n\\bibliographystyle{unsrt}\n\\bibliography{references}")

# capitalised cross-references, as in the template ("See Section \ref", "Figure \ref", "Table~\ref")
caps = [(r"\bsections~(?=\\ref)", "Sections~"), (r"\bsection~(?=\\ref)", "Section~"),
        (r"\bfigures~(?=\\ref|[0-9S])", "Figures~"), (r"\bfigure~(?=\\ref|[0-9S])", "Figure~"),
        (r"\btable~(?=\\ref)", "Table~"), (r"\bappendix~(?=\\ref)", "Appendix~"),
        (r"\beqs\.~(?=\\eqref|[0-9S(IVX])", "Eqs.~"), (r"\beq\.~(?=\\eqref|[0-9S(IVX])", "Eq.~"),
        (r"\brefs\.~(?=\\cite)", "Refs.~"), (r"\bref\.~(?=\\cite)", "Ref.~")]
for pat, rep in caps:
    body = re.sub(pat, rep, body)

# table captions above the tabular, as in the template
def caption_on_top(block):
    i = block.index("\\caption{")
    cap, end = braced(block, i + len("\\caption"))
    lab = re.match(r"\s*\\label\{[^}]*\}", block[end:])
    label = lab.group(0).strip()
    block = block[:i] + block[end + lab.end():]
    block = block.replace("\\centering\\small\n", "\\caption{" + cap + "}\n" + label + "\n\\centering\\small\n", 1)
    return re.sub(r"\n\n+\\end\{table\}", lambda m: "\n\\end{table}", block)

body = re.sub(r"\\begin\{table\}\[t\].*?\\end\{table\}", lambda m: caption_on_top(m.group(0)), body, flags=re.DOTALL)

preamble = r"""\documentclass{article}


\usepackage{PRIMEarxiv}

\usepackage[utf8]{inputenc} % allow utf-8 input
\usepackage[T1]{fontenc}    % use 8-bit T1 fonts
\usepackage{amsmath,amssymb,bm} % mathematics
\usepackage{hyperref}       % hyperlinks
\usepackage{url}            % simple URL typesetting
\usepackage{booktabs}       % professional-quality tables
\usepackage{array}          % column types for the conventions table
\usepackage{amsfonts}       % blackboard math symbols
\usepackage{nicefrac}       % compact symbols for 1/2, etc.
\usepackage{microtype}      % microtypography
\usepackage{fancyhdr}       % header
\usepackage{graphicx}       % graphics
\usepackage{orcidlink}      % ORCID logo
\graphicspath{{media/}}     % the figures are under media/
\urlstyle{same}             % DOIs in the text font, breakable at / and .
\hypersetup{colorlinks=true,linkcolor=blue,citecolor=blue,urlcolor=blue,
  pdftitle={Peaks or thresholds? Non-Gaussian primordial black hole abundances for broad spectra and
  the pulsar-timing-array tension},pdfauthor={Duong Ngoc Khoa}}
\setlength{\emergencystretch}{1.5em}

%Header
\pagestyle{fancy}
\thispagestyle{empty}
\rhead{ \textit{ }}

% Update your Headers here
\fancyhead[LO]{Peaks or thresholds? Non-Gaussian PBH abundances and the PTA tension}

""" + macros + r"""

%% Title
\title{Peaks or thresholds? Non-Gaussian primordial black hole abundances for broad spectra and the
pulsar-timing-array tension
%%%% Cite as
%%%% Update your official citation here when published
%\thanks{\textit{\underline{Citation}}: \textbf{D.~N.~Khoa. Peaks or thresholds? ... DOI:...}}
}

\author{
  Duong Ngoc Khoa\,\orcidlink{0009-0009-9409-3298} \\
  Phenikaa Institute for Advanced Study \\
  Phenikaa University \\
  Hanoi, Vietnam \\
  \texttt{yuutoduong@gmail.com} \\
}


\begin{document}
\maketitle


\begin{abstract}
""" + abstract.strip() + r"""
\end{abstract}


% keywords can be removed
\keywords{Primordial black holes \and Non-Gaussianity \and Peaks theory \and Pulsar timing arrays \and
Induced gravitational waves}


"""
Path("docs/arxiv/ngpbh_arxiv.tex").write_text(preamble + body.lstrip("\n") + "\\end{document}\n")

# bibliography for unsrt: arXiv numbers (linked) in the note, arXiv-only entries as @misc
bib = Path("docs/paper/refs.bib").read_text()
entries = re.split(r"\n(?=@)", bib)
out = ["% References of \"Peaks or thresholds?\" (D. N. Khoa), from INSPIRE-HEP, adapted to unsrt:",
       "% arXiv numbers in the note field, arXiv-only papers as misc entries.", ""]
for e in entries:
    if not e.lstrip().startswith("@"):
        continue
    e = re.sub(r"\n\s*month = \"[^\"]*\",?", "", e.strip())
    assert e.endswith("\n}")
    ep = re.search(r'eprint = "([^"]+)"', e)
    has_journal = re.search(r'\bjournal = "', e) is not None
    if ep:
        link = "\\href{https://arxiv.org/abs/" + ep.group(1) + "}{arXiv:" + ep.group(1) + "}"
        if has_journal:
            m = re.search(r'note = "([^"]*)"', e)
            if m:
                e = e.replace(m.group(0), 'note = "' + link + " " + m.group(1) + '"')
            else:
                e = e[:-2] + ',\n    note = "' + link + '"\n}'  # e ends with "\n}"
        else:
            e = re.sub(r"^@article\{", "@misc{", e)
            e = e[:-2] + ',\n    howpublished = "' + link + '"\n}'
    e = re.sub(r",\s*,", ",", e)
    out += [e, ""]
Path("docs/arxiv/references.bib").write_text("\n".join(out))
print("wrote docs/arxiv/ngpbh_arxiv.tex and docs/arxiv/references.bib")
