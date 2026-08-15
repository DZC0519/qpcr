from __future__ import annotations
import math, re
from dataclasses import dataclass
from pathlib import Path
from typing import Any
import pandas as pd
from .models import WellRecord

ALIASES = {
 'well': {'pos', 'well', 'position', '\u5b54\u4f4d', '\u5b54', '\u5b54\u53f7'},
 'ct': {'cp', 'cq', 'ct', 'cycle', 'crossing point', '\u5b9a\u91cf\u5faa\u73af'},
 'include': {'include', 'included', '\u4f7f\u7528', '\u5305\u542b'},
 'color': {'color', 'colour', '\u67d3\u6599', '\u989c\u8272'},
 'name': {'name', 'sample name', 'sample', '\u540d\u79f0', '\u6837\u672c'},
 'concentration': {'concentration', 'conc', '\u6d53\u5ea6'},
 'standard': {'standard', 'std', '\u6807\u51c6'},
 'status': {'status', 'call', '\u72b6\u6001'}}

@dataclass(slots=True)
class ImportCandidate:
    sheet: str
    header_row: int
    columns: dict[str, str]

@dataclass(slots=True)
class ImportedPlate:
    source_path: str
    sheet: str
    header_row: int
    wells: tuple[WellRecord, ...]
    plate_format: int

def _norm(value: Any) -> str:
    return re.sub(r'[\s_\-]+', ' ', str(value).strip().lower())

def _column_map(values: list[Any]) -> dict[str, str]:
    found = {}
    for value in values:
        for canonical, aliases in ALIASES.items():
            if _norm(value) in aliases and canonical not in found:
                found[canonical] = str(value)
    return found

def find_candidates(path: str | Path, scan_rows: int = 30) -> list[ImportCandidate]:
    path = Path(path)
    sheets = ({path.stem: pd.read_csv(path, header=None, nrows=scan_rows)}
              if path.suffix.lower() == '.csv' else
              pd.read_excel(path, sheet_name=None, header=None, nrows=scan_rows))
    candidates = []
    for sheet, preview in sheets.items():
        for index, row in preview.iterrows():
            columns = _column_map(row.dropna().tolist())
            if 'well' in columns and 'ct' in columns:
                candidates.append(ImportCandidate(str(sheet), int(index) + 1, columns))
                break
    return candidates

def _clean_well(value: Any) -> str:
    text = str(value).strip().upper().replace(' ', '')
    match = re.fullmatch(r'([A-Z]+)0*(\d+)', text)
    return f'{match.group(1)}{int(match.group(2))}' if match else text

def _optional_float(value: Any) -> float | None:
    if value is None or pd.isna(value) or str(value).strip() == '': return None
    try:
        number = float(value)
        return number if math.isfinite(number) else None
    except (TypeError, ValueError): return None

def _truthy(value: Any) -> bool:
    if value is None or pd.isna(value): return True
    return _norm(value) not in {'false', 'no', '0', 'n', '\u5426'}

def _optional_text(value: Any) -> str:
    return '' if value is None or pd.isna(value) else str(value).strip()

def load_plate(path: str | Path, sheet: str | None = None,
               header_row: int | None = None) -> ImportedPlate:
    path = Path(path)
    candidates = find_candidates(path)
    if not candidates:
        raise ValueError('\u672a\u627e\u5230\u540c\u65f6\u5305\u542b\u5b54\u4f4d\u548c Cp/Cq/Ct \u7684\u8868\u5934\u3002')
    candidate = next((c for c in candidates if sheet is None or c.sheet == sheet), candidates[0])
    header_row = header_row or candidate.header_row
    frame = (pd.read_csv(path, header=header_row - 1) if path.suffix.lower() == '.csv' else
             pd.read_excel(path, sheet_name=candidate.sheet, header=header_row - 1))
    columns, records = _column_map(frame.columns.tolist()), []
    for offset, row in frame.iterrows():
        well_value = row.get(columns['well'])
        if well_value is None or pd.isna(well_value) or not str(well_value).strip(): continue
        raw = {str(k): (None if pd.isna(v) else v) for k, v in row.items()}
        records.append(WellRecord(
            well=_clean_well(well_value), ct=_optional_float(row.get(columns['ct'])),
            include=_truthy(row.get(columns.get('include', ''))),
            color=row.get(columns.get('color', ''), None),
            machine_name=_optional_text(row.get(columns.get('name', ''), '')),
            concentration=_optional_float(row.get(columns.get('concentration', ''))),
            standard=row.get(columns.get('standard', ''), None),
            status=_optional_text(row.get(columns.get('status', ''), '')),
            source_sheet=candidate.sheet, source_row=int(offset)+header_row+1, raw=raw))
    if not records:
        raise ValueError('\u8868\u5934\u5df2\u627e\u5230\uff0c\u4f46\u6ca1\u6709\u8bfb\u53d6\u5230\u5b54\u6570\u636e\u3002')
    rows = {re.match(r'[A-Z]+', r.well).group() for r in records if re.match(r'[A-Z]+', r.well)}
    cols = {int(re.search(r'\d+', r.well).group()) for r in records if re.search(r'\d+', r.well)}
    plate_format = 384 if len(rows) > 8 or max(cols, default=0) > 12 else 96
    return ImportedPlate(str(path), candidate.sheet, header_row, tuple(records), plate_format)
