"""Dump section headings, paragraphs and figure captions from a PMC article XML.

Usage: python -I dump_text.py <path-to-pmc.xml> [grep-substring]
The XML is untrusted downloaded data and is only parsed, never executed.
"""
import re
import sys
import xml.etree.ElementTree as ET


def txt(e):
    return re.sub(r"\s+", " ", "".join(e.itertext())).strip()


def main(path, needle=None):
    root = ET.parse(path).getroot()
    for fig in root.iter("fig"):
        lab = fig.find("label")
        print("---- FIGURE", txt(lab) if lab is not None else "?", "----")
        cap = fig.find("caption")
        print(txt(cap) if cap is not None else "")
        print()
    body = root.find(".//body")
    if body is None:
        return
    for sec in body.iter("sec"):
        ti = sec.find("title")
        head = txt(ti) if ti is not None else ""
        paras = [txt(p) for p in sec.findall("p")]
        if needle:
            paras = [p for p in paras if needle.lower() in p.lower()]
            if not paras:
                continue
        print("##", head)
        for p in paras:
            print(p)
            print()


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)
