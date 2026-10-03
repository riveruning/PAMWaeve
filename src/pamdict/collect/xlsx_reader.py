"""Minimal dependency-free XLSX reader (XLSX = zip of XML).

Reads shared strings and converts a worksheet to a list of rows (list of cell
values). Avoids openpyxl so the collector stays install-free in the sandbox
where pip is unreliable.
"""

from __future__ import annotations

import re
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

_NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
_REL_NS = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"
_REL_ID = f"{_REL_NS}id"


def _col_to_idx(ref: str) -> int:
    """Convert a column letter sequence ('A'..'ZZ') to a 0-based index."""
    n = 0
    for ch in ref:
        n = n * 26 + (ord(ch.upper()) - ord("A") + 1)
    return n - 1


def _cell_ref(coord: str) -> tuple[int, int]:
    m = re.match(r"([A-Z]+)(\d+)", coord)
    return int(m.group(2)) - 1, _col_to_idx(m.group(1))  # (row, col)


def read_sheet(path: str | Path, sheet_name: str | None = None) -> list[list[str]]:
    """Return a worksheet as rows of string cell values."""
    z = zipfile.ZipFile(str(path))

    # shared strings
    shared: list[str] = []
    if "xl/sharedStrings.xml" in z.namelist():
        root = ET.fromstring(z.read("xl/sharedStrings.xml"))
        for si in root.iter(f"{_NS}si"):
            shared.append("".join(t.text or "" for t in si.iter(f"{_NS}t")))

    # pick sheet path
    wb = ET.fromstring(z.read("xl/workbook.xml"))
    sheets_el = wb.find(f"{_NS}sheets")
    sheet_els = [s for s in sheets_el.iter(f"{_NS}sheet")]
    if sheet_name is None:
        sheet_name = sheet_els[0].get("name")
    rels_el = ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
    rel_map = {
        rel.get("Id"): rel.get("Target")
        for rel in rels_el.iter(f"{{http://schemas.openxmlformats.org/package/2006/relationships}}Relationship")
    }
    rid = None
    for s in sheet_els:
        if s.get("name") == sheet_name:
            rid = s.get(_REL_ID)
            break
    if rid is None:
        rid = sheet_els[0].get(_REL_ID)
    target = rel_map.get(rid, "worksheets/sheet1.xml")
    if not target.startswith("xl/"):
        sheet_path = "xl/" + target.lstrip("/")
    else:
        sheet_path = target

    root = ET.fromstring(z.read(sheet_path))
    rows: dict[int, dict[int, str]] = {}
    for c in root.iter(f"{_NS}c"):
        coord = c.get("r")
        if not coord:
            continue
        t = c.get("t")
        v = c.find(f"{_NS}v")
        val = ""
        if t == "s" and v is not None:
            val = shared[int(v.text)]
        elif t == "inlineStr":
            val = "".join(x.text or "" for x in c.iter(f"{_NS}t"))
        elif v is not None:
            val = v.text or ""
        r, col = _cell_ref(coord)
        rows.setdefault(r, {})[col] = val

    max_row = max(rows) if rows else -1
    max_col = max((c for row in rows.values() for c in row), default=-1)
    out = []
    for r in range(max_row + 1):
        out.append([rows.get(r, {}).get(c, "") for c in range(max_col + 1)])
    return out


def list_sheets(path: str | Path) -> list[str]:
    z = zipfile.ZipFile(str(path))
    wb = ET.fromstring(z.read("xl/workbook.xml"))
    return [s.get("name") for s in wb.find(f"{_NS}sheets").iter(f"{_NS}sheet")]
