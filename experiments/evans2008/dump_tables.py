"""Dump every table-wrap in a PMC article XML as pipe-delimited rows.

Usage: python -I dump_tables.py <path-to-pmc.xml>
The XML is untrusted downloaded data and is only parsed, never executed.
"""
import re
import sys
import xml.etree.ElementTree as ET


def txt(e):
    return re.sub(r"\s+", " ", "".join(e.itertext())).strip()


def main(path):
    root = ET.parse(path).getroot()
    for tw in root.iter("table-wrap"):
        lab = tw.find("label")
        cap = tw.find("caption")
        print("=" * 100)
        print("LABEL:", txt(lab) if lab is not None else "?")
        print("CAPTION:", txt(cap) if cap is not None else "?")
        for tbl in tw.iter("table"):
            for row in tbl.iter("tr"):
                cells = []
                for c in row:
                    if c.tag in ("td", "th"):
                        s = txt(c)
                        span = c.get("colspan")
                        if span and span != "1":
                            s += "[cs=%s]" % span
                        cells.append(s)
                print(" | ".join(cells))
        for f in tw.findall(".//table-wrap-foot"):
            print("FOOT:", txt(f))


if __name__ == "__main__":
    main(sys.argv[1])
