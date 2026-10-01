"""
Export the dissertation to a Word document for review.

Converts the LaTeX sources under thesis/ into a .docx with resolved cross-references,
author-year citations and a reference list built from dissertation.bib, so that a reader
without a LaTeX toolchain can read and comment on the text. The PDF built from the LaTeX
sources remains the authoritative document: this export carries the text, tables and
figures, not the template's layout.

Usage:
    python scripts/thesis_to_docx.py [out.docx]
"""
import re
import sys
import unicodedata
from pathlib import Path

from docx import Document
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.shared import Inches, Pt, RGBColor

BASE_DIR = Path(__file__).resolve().parent.parent
THESIS = BASE_DIR / "thesis"
ORDER = ["preamble/Abstract.tex", "chapters/Introduction.tex", "chapters/StateOfTheArt.tex",
         "chapters/Methods.tex", "chapters/ResultsAndDiscussion.tex", "chapters/Conclusion.tex"]

# ---------------------------------------------------------------- LaTeX -> text
SYMBOLS = {
    r"\AA": "Å", r"\times": "×", r"\lambda": "λ", r"\Delta": "Δ", r"\mu": "µ",
    r"\pi": "π", r"\varepsilon": "ε", r"\sqrt": "√", r"\geq": "≥", r"\leq": "≤",
    r"\pm": "±", r"\sim": "~", r"\rightarrow": "→", r"\cdot": "·", r"\ldots": "…",
    r"\circ": "°", r"\%": "%", r"\&": "&", r"\_": "_", r"\#": "#", r"\$": "$",
    r"\,": "\u202f", r"\;": " ", r"\!": "", r"\:": " ", r"\ ": " ", r"\\": " ",
}
SUBSCRIPT = {"0": "₀", "1": "₁", "2": "₂", "3": "₃", "4": "₄", "5": "₅", "6": "₆",
             "7": "₇", "8": "₈", "9": "₉", "D": "_D", "m": "_m"}


SUPERSCRIPT = str.maketrans("0123456789+-=()n", "⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻⁼⁽⁾ⁿ")


def superscript(t):
    """Unicode superscript where every character has one; otherwise a caret-free fallback."""
    t = t.replace("\\", "").strip()
    out = t.translate(SUPERSCRIPT)
    return out if all(c not in "0123456789+-=()n" for c in out) else out


def strip_math(s):
    """Render inline math as readable plain text."""
    def one(m):
        t = m.group(1)
        t = re.sub(r"\\text(?:rm|it|bf)?\{([^{}]*)\}", r"\1", t)
        t = re.sub(r"\\mathrm\{([^{}]*)\}", r"\1", t)
        for k, v in SYMBOLS.items():
            t = t.replace(k, v)
        t = re.sub(r"_\{?([0-9A-Za-z]+)\}?", lambda g: "".join(SUBSCRIPT.get(c, c) for c in g.group(1)), t)
        t = re.sub(r"\^\{([^{}]*)\}|\^(\S)",
                   lambda g: superscript(g.group(1) if g.group(1) is not None else g.group(2)), t)
        t = t.replace("{", "").replace("}", "").replace("\\", "")
        return t.strip()
    return re.sub(r"\$([^$]*)\$", one, s)


def clean(s, refs=None, cites=None, acronyms=None):
    """LaTeX fragment -> plain text, resolving refs, citations and acronyms."""
    s = re.sub(r"(?<!\\)%.*", "", s)
    s = s.replace(r"\,^{\circ}", "°").replace(r"^{\circ}", "°")
    s = strip_math(s)
    if acronyms:
        s = re.sub(r"\\acrfull\{([^}]+)\}",
                   lambda m: f"{acronyms.get(m.group(1), ('', m.group(1)))[1]} "
                             f"({acronyms.get(m.group(1), (m.group(1), ''))[0]})", s)
        s = re.sub(r"\\acr(?:short|long)\{([^}]+)\}",
                   lambda m: acronyms.get(m.group(1), (m.group(1), ""))[0], s)
    if refs is not None:
        s = re.sub(r"\\(?:ref|autoref)\{([^}]+)\}", lambda m: refs.get(m.group(1), "?"), s)
    if cites is not None:
        def citep(m):
            keys = [k.strip() for k in m.group(1).split(",")]
            return "(" + "; ".join(cites.get(k, k) for k in keys) + ")"

        def citet(m):
            keys = [k.strip() for k in m.group(1).split(",")]
            return ", ".join(cites.get(k, k).replace(", ", " (") + ")" for k in keys)
        s = re.sub(r"\\citep\{([^}]+)\}", citep, s)
        s = re.sub(r"\\citet\{([^}]+)\}", citet, s)
        s = re.sub(r"\\cite\{([^}]+)\}", citep, s)
    s = re.sub(r"\\(?:label|index|vspace|hspace|noindent|centering|cleardoublepage|newpage)\*?\{?[^}\n]*\}?", "", s)
    s = re.sub(r"\\(?:emph|textit|textbf|textsc|texttt|text|mbox|underline)\{([^{}]*)\}", r"\1", s)
    s = re.sub(r"\\(?:emph|textit|textbf|textsc|texttt)\{(.*?)\}", r"\1", s)
    s = re.sub(r"\\url\{([^}]*)\}", r"\1", s)
    s = re.sub(r"\\footnote\{([^{}]*)\}", r" [\1]", s)
    for k, v in SYMBOLS.items():
        s = s.replace(k, v)
    s = s.replace("---", "—").replace("--", "–")
    s = s.replace("``", "\u201c").replace("''", "\u201d").replace("~", " ")
    s = re.sub(r"\\[A-Za-z]+\*?", "", s)
    s = s.replace("{", "").replace("}", "")
    return re.sub(r"[ \t]+", " ", s).strip()


# ---------------------------------------------------------------- bibliography
def parse_bib(path):
    txt = path.read_text(encoding="utf-8")
    entries = {}
    for m in re.finditer(r"@(\w+)\s*\{\s*([^,]+),", txt):
        start = m.end()
        depth, i = 1, start
        while depth and i < len(txt):
            depth += (txt[i] == "{") - (txt[i] == "}")
            i += 1
        body = txt[start:i - 1]
        fields = {}
        for fm in re.finditer(r"(\w+)\s*=\s*", body):
            j = fm.end()
            if j >= len(body):
                continue
            if body[j] == "{":
                depth, k = 1, j + 1
                while depth and k < len(body):
                    depth += (body[k] == "{") - (body[k] == "}")
                    k += 1
                val = body[j + 1:k - 1]
            else:
                k = body.find(",", j)
                val = body[j:k if k > 0 else len(body)]
            fields[fm.group(1).lower()] = strip_math(re.sub(r"\s+", " ", val).strip())
        entries[m.group(2).strip()] = fields
    return entries


def delatex(s):
    rep = {r"{\'e}": "é", r"{\^o}": "ô", r"{\'a}": "á", r"{\'i}": "í", r"{\'o}": "ó",
           r"{\'u}": "ú", r"{\~a}": "ã", r"{\~n}": "ñ", r"{\`e}": "è", r'{\"o}': "ö",
           r'{\"u}': "ü", r'{\"a}': "ä", r"\c{c}": "ç", r"\c{t}": "ţ", r"\u{a}": "ă"}
    for k, v in rep.items():
        s = s.replace(k, v)
    return re.sub(r"[{}\\]", "", s).replace("~", " ")


def surname(author_field):
    first = delatex(author_field).split(" and ")[0].strip()
    return first.split(",")[0].strip() if "," in first else first.split()[-1]


def cite_label(fields):
    authors = delatex(fields.get("author", "")).split(" and ")
    year = fields.get("year", "n.d.")
    if not authors or not authors[0]:
        return year
    s = surname(fields["author"])
    if len(authors) == 1:
        return f"{s}, {year}"
    if len(authors) == 2:
        second = authors[1].strip()
        s2 = second.split(",")[0].strip() if "," in second else second.split()[-1]
        return f"{s} & {s2}, {year}"
    return f"{s} et al., {year}"


def reference_entry(fields):
    authors = [a.strip() for a in delatex(fields.get("author", "")).split(" and ") if a.strip()]
    def fmt(a):
        if "," in a:
            fam, giv = [x.strip() for x in a.split(",", 1)]
            ini = " ".join(f"{p[0]}." for p in giv.replace(".", " ").split() if p)
            return f"{fam}, {ini}"
        parts = a.split()
        return f"{parts[-1]}, " + " ".join(f"{p[0]}." for p in parts[:-1]) if len(parts) > 1 else a
    names = "; ".join(fmt(a) for a in authors) if authors else ""
    out = f"{names} ({fields.get('year', 'n.d.')}). {delatex(fields.get('title', ''))}."
    jour = delatex(fields.get("journal", "") or fields.get("booktitle", ""))
    if jour:
        out += f" {jour}"
        if fields.get("volume"):
            out += f", {fields['volume']}"
            if fields.get("number"):
                out += f"({fields['number']})"
        if fields.get("pages"):
            out += f", {delatex(fields['pages']).replace('--', '–')}"
        out += "."
    if fields.get("howpublished"):
        out += f" {delatex(fields['howpublished'])}."
    if fields.get("doi"):
        out += f" https://doi.org/{fields['doi']}"
    return re.sub(r"\s+", " ", out).strip()


# ---------------------------------------------------------------- pass 1: numbering
def collect(sources):
    """Resolve \\label -> '3.2', 'Table 4.1', 'Figure 4.2' and read the acronyms."""
    refs, acronyms = {}, {}
    chap = 0
    sec = sub = 0
    tab = fig = 0
    for rel in sources:
        text = (THESIS / rel).read_text(encoding="utf-8")
        for m in re.finditer(r"\\newacronym\{([^}]+)\}\{([^}]+)\}\{([^}]+)\}", text):
            acronyms[m.group(1)] = (strip_math(m.group(2)), m.group(3))
        pending = None                      # what the next \label refers to
        for line in text.splitlines():
            if re.match(r"\s*\\chapter\*", line):
                sec = sub = tab = fig = 0; pending = None
            elif re.match(r"\s*\\chapter\b", line):
                chap += 1; sec = sub = tab = fig = 0; pending = ("sec", f"{chap}")
            elif re.match(r"\s*\\section\b", line):
                sec += 1; sub = 0; pending = ("sec", f"{chap}.{sec}")
            elif re.match(r"\s*\\subsection\b", line):
                sub += 1; pending = ("sec", f"{chap}.{sec}.{sub}")
            elif r"\begin{table}" in line:
                tab += 1; pending = ("tab", f"{chap}.{tab}")
            elif r"\begin{figure}" in line:
                fig += 1; pending = ("fig", f"{chap}.{fig}")
            for m in re.finditer(r"\\label\{([^}]+)\}", line):
                if pending:
                    kind, val = pending
                    refs[m.group(1)] = val if kind == "sec" else val
    return refs, acronyms


# ---------------------------------------------------------------- pass 2: build
def add_caption(doc, text, bold_prefix):
    p = doc.add_paragraph()
    r = p.add_run(bold_prefix + ". ")
    r.bold = True
    r.font.size = Pt(9)
    r2 = p.add_run(text)
    r2.font.size = Pt(9)
    p.paragraph_format.space_after = Pt(12)
    return p


def emit_table(doc, body, refs, cites, acronyms, number):
    rows = []
    caption = ""
    for m in re.finditer(r"\\caption\{(.*?)\}\s*(?:\\label|\\end\{table\})", body, re.S):
        caption = clean(m.group(1), refs, cites, acronyms)
        break
    tab = re.search(r"\\begin\{tabular\}\{[^}]*\}(.*?)\\end\{tabular\}", body, re.S)
    if not tab:
        return
    for raw in tab.group(1).split(r"\\"):
        raw = re.sub(r"\\(?:hline|toprule|midrule|bottomrule)", "", raw).strip()
        if not raw:
            continue
        cells = []
        for c in raw.split("&"):
            c = re.sub(r"\\multicolumn\{(\d+)\}\{[^}]*\}\{(.*)\}", r"\2", c.strip(), flags=re.S)
            cells.append(clean(c, refs, cites, acronyms))
        if any(cells):
            rows.append(cells)
    if not rows:
        return
    ncol = max(len(r) for r in rows)
    t = doc.add_table(rows=0, cols=ncol)
    t.style = "Table Grid"
    t.alignment = WD_TABLE_ALIGNMENT.CENTER
    for i, r in enumerate(rows):
        cells = t.add_row().cells
        for j in range(ncol):
            txt = r[j] if j < len(r) else ""
            cells[j].text = txt
            for par in cells[j].paragraphs:
                for run in par.runs:
                    run.font.size = Pt(9)
                    if i == 0:
                        run.bold = True
    if caption:
        add_caption(doc, caption, f"Table {number}")


def emit_figure(doc, body, refs, cites, acronyms, number):
    caption = ""
    m = re.search(r"\\caption\{(.*?)\}\s*(?:\\label|\\end\{figure\})", body, re.S)
    if m:
        caption = clean(m.group(1), refs, cites, acronyms)
    for img in re.findall(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}", body):
        path = THESIS / img
        if not path.exists():
            continue
        p = doc.add_paragraph()
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        p.add_run().add_picture(str(path), width=Inches(6.0))
    if caption:
        add_caption(doc, caption, f"Figure {number}")


def main(out_path):
    refs, acronyms = collect(ORDER)
    bib = parse_bib(THESIS / "dissertation.bib")
    cites = {k: cite_label(v) for k, v in bib.items()}

    doc = Document()
    st = doc.styles["Normal"]
    st.font.name = "Calibri"
    st.font.size = Pt(11)
    st.paragraph_format.space_after = Pt(8)

    title = re.search(r"\\titleA\{([^}]*)\}\s*\n\s*\\titleB\{([^}]*)\}",
                      (THESIS / "dissertation.tex").read_text(encoding="utf-8"))
    h = doc.add_heading(f"{title.group(1)} {title.group(2)}" if title else "Dissertation", level=0)
    sub = doc.add_paragraph("Romeu Alexandre Ribeiro Fernandes — Master's Dissertation in Bioinformatics, "
                            "University of Minho")
    sub.runs[0].italic = True
    note = doc.add_paragraph("Working copy exported from the LaTeX sources for review; the compiled PDF "
                             "remains the authoritative document.")
    note.runs[0].font.size = Pt(9)
    note.runs[0].font.color.rgb = RGBColor(0x66, 0x66, 0x66)

    used = set()
    chap = tab_n = fig_n = 0
    for rel in ORDER:
        text = (THESIS / rel).read_text(encoding="utf-8")
        text = re.sub(r"(?m)^\s*%.*$", "", text)
        text = re.sub(r"\\newacronym\{[^}]*\}\{[^}]*\}\{[^}]*\}", "", text)
        for key in re.findall(r"\\cite[pt]?\{([^}]+)\}", text):
            used.update(k.strip() for k in key.split(","))

        i = 0
        buf = []

        def flush():
            if not buf:
                return
            txt = clean(" ".join(buf), refs, cites, acronyms)
            buf.clear()
            if txt:
                doc.add_paragraph(txt)

        lines = text.split("\n")
        while i < len(lines):
            line = lines[i]
            mchap = re.match(r"\s*\\chapter\*?\{(.*)\}", line)
            msec = re.match(r"\s*\\section\*?\{(.*)\}", line)
            msub = re.match(r"\s*\\subsection\*?\{(.*)\}", line)
            mpar = re.match(r"\s*\\paragraph\{(.*?)\}\s*(.*)", line)
            if mchap:
                flush()
                if "\\chapter*" not in line:
                    chap += 1
                tab_n = fig_n = 0
                doc.add_page_break()
                doc.add_heading(clean(mchap.group(1), refs, cites, acronyms), level=1)
            elif msec:
                flush(); doc.add_heading(clean(msec.group(1), refs, cites, acronyms), level=2)
            elif msub:
                flush(); doc.add_heading(clean(msub.group(1), refs, cites, acronyms), level=3)
            elif mpar:
                flush()
                p = doc.add_paragraph()
                r = p.add_run(clean(mpar.group(1), refs, cites, acronyms) + " ")
                r.bold = True
                rest = clean(mpar.group(2), refs, cites, acronyms)
                j = i + 1
                while j < len(lines) and lines[j].strip() and not lines[j].lstrip().startswith("\\"):
                    rest += " " + clean(lines[j], refs, cites, acronyms)
                    j += 1
                i = j - 1
                p.add_run(rest)
            elif r"\begin{table}" in line:
                flush()
                j = i
                while r"\end{table}" not in lines[j]:
                    j += 1
                tab_n += 1
                emit_table(doc, "\n".join(lines[i:j + 1]), refs, cites, acronyms, f"{chap}.{tab_n}")
                i = j
            elif r"\begin{figure}" in line:
                flush()
                j = i
                while r"\end{figure}" not in lines[j]:
                    j += 1
                fig_n += 1
                emit_figure(doc, "\n".join(lines[i:j + 1]), refs, cites, acronyms, f"{chap}.{fig_n}")
                i = j
            elif r"\begin{itemize}" in line or r"\begin{enumerate}" in line:
                flush()
                style = "List Bullet" if "itemize" in line else "List Number"
                j = i + 1
                item = []
                while r"\end{itemize}" not in lines[j] and r"\end{enumerate}" not in lines[j]:
                    if lines[j].lstrip().startswith(r"\item"):
                        if item:
                            doc.add_paragraph(clean(" ".join(item), refs, cites, acronyms), style=style)
                        item = [lines[j].lstrip()[5:]]
                    elif item:
                        item.append(lines[j])
                    j += 1
                if item:
                    doc.add_paragraph(clean(" ".join(item), refs, cites, acronyms), style=style)
                i = j
            elif re.match(r"\s*\\\[", line):                       # display math
                flush()
                j = i
                chunk = []
                while r"\]" not in lines[j]:
                    chunk.append(lines[j]); j += 1
                chunk.append(lines[j])
                body = " ".join(chunk).replace(r"\[", "").replace(r"\]", "")
                p = doc.add_paragraph(clean("$" + body + "$", refs, cites, acronyms))
                p.alignment = WD_ALIGN_PARAGRAPH.CENTER
                i = j
            elif not line.strip():
                flush()
            elif re.match(r"\s*\\(begin|end|include|input|clearpage|cleardoublepage)", line):
                flush()
            else:
                buf.append(line)
            i += 1
        flush()

    doc.add_page_break()
    doc.add_heading("References", level=1)
    for key in sorted(used, key=lambda k: (surname(bib[k].get("author", "")) if k in bib else k,
                                           bib.get(k, {}).get("year", ""))):
        if key in bib:
            p = doc.add_paragraph(reference_entry(bib[key]))
            p.paragraph_format.left_indent = Inches(0.4)
            p.paragraph_format.first_line_indent = Inches(-0.4)
            p.paragraph_format.space_after = Pt(6)
    doc.save(out_path)
    print(f"wrote {out_path}")
    print(f"  {chap} chapters, {len(used)} references cited, {len(refs)} cross-references resolved")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else str(BASE_DIR / "Thesis.docx"))
