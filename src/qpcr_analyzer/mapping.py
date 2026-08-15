from __future__ import annotations

import json
import re
from dataclasses import asdict
from pathlib import Path

from .models import GROUP_ORDER, PlateAssignment, WellRecord, normalize_group_order


def well_sort_key(well: str) -> tuple[int, int]:
    letters = ''.join(c for c in well.upper() if c.isalpha())
    number = int(''.join(c for c in well if c.isdigit()) or 0)
    row = 0
    for char in letters:
        row = row * 26 + ord(char) - 64
    return row, number


def plate_row(well: str) -> str:
    match = re.match(r'([A-Z]+)', well.upper())
    return match.group(1) if match else ''


def normalize_sample_axis(value: str | None, default: str = 'row') -> str:
    return value if value in {'auto', 'row', 'column'} else default


def normalize_mapping_scope(value: str | None, default: str = 'axis') -> str:
    return value if value in {'axis', 'well'} else default


def sample_axis_key(well: str, sample_axis: str) -> str:
    if normalize_sample_axis(sample_axis) == 'column':
        return str(well_sort_key(well)[1])
    return plate_row(well)


def detect_sample_axis(
    wells: list[str], purpose: str, plate_rows: int = 0, plate_columns: int = 0,
) -> str | None:
    """Infer the sample direction from a deliberate horizontal or vertical selection."""
    if not wells or purpose not in {'group', 'gene'}:
        return None
    rows: dict[str, set[int]] = {}
    columns: dict[int, set[int]] = {}
    for well in wells:
        row_index, column = well_sort_key(well)
        rows.setdefault(plate_row(well), set()).add(column)
        columns.setdefault(column, set()).add(row_index)
    def longest_contiguous_run(values: set[int]) -> int:
        ordered = sorted(values)
        longest = current = 0
        previous = None
        for value in ordered:
            current = current + 1 if previous is not None and value == previous + 1 else 1
            longest = max(longest, current)
            previous = value
        return longest

    horizontal = max((longest_contiguous_run(values) for values in rows.values()), default=0)
    vertical = max((longest_contiguous_run(values) for values in columns.values()), default=0)

    if plate_columns and horizontal == plate_columns and vertical < plate_rows:
        return 'row'
    if plate_rows and vertical == plate_rows and horizontal < plate_columns:
        return 'column'
    if horizontal == vertical:
        return None
    # A horizontal strip represents one row-oriented sample; a vertical
    # strip represents one column-oriented sample.  The same geometry is
    # valid for both group assignment and contiguous technical replicates.
    return 'row' if horizontal > vertical else 'column'


def suggest_row_layout(
    wells: tuple[WellRecord, ...], empty_rows: tuple[str, ...] = ('O', 'P')
) -> dict[str, PlateAssignment]:
    assignments: dict[str, PlateAssignment] = {}
    for record in wells:
        row_name = plate_row(record.well)
        if row_name in empty_rows:
            assignments[record.well] = PlateAssignment(
                record.well, role='empty', plate_row=row_name)
        else:
            assignments[record.well] = PlateAssignment(
                well=record.well, sample=f'样本{row_name}', role='sample', plate_row=row_name)
    return assignments


def suggest_triplicate_layout(
    wells: tuple[WellRecord, ...], empty_rows: tuple[str, ...] = ('O', 'P')
) -> dict[str, PlateAssignment]:
    assignments = suggest_row_layout(wells, empty_rows)
    for record in wells:
        assignment = assignments[record.well]
        if assignment.role != 'sample':
            continue
        column = well_sort_key(record.well)[1]
        assignment.group = 'Ctrl'
        assignment.gene = f'基因{(column - 1) // 3 + 1}'
        assignment.technical_replicate = (column - 1) % 3 + 1
    renumber_samples(assignments)
    return assignments


def renumber_samples(
    assignments: dict[str, PlateAssignment], sample_axis: str = 'row',
) -> None:
    sample_axis = normalize_sample_axis(sample_axis)
    if sample_axis == 'auto':
        return
    axis_groups: dict[str, str] = {}
    for assignment in sorted(assignments.values(), key=lambda a: well_sort_key(a.well)):
        if assignment.role == 'sample' and assignment.group:
            key = sample_axis_key(assignment.well, sample_axis)
            axis_groups.setdefault(key, assignment.group)
    counters = {group: 0 for group in GROUP_ORDER}
    sample_by_axis: dict[str, str] = {}
    for key, group in axis_groups.items():
        counters.setdefault(group, 0)
        counters[group] += 1
        sample_by_axis[key] = f'{group}-{counters[group]}'
    for assignment in assignments.values():
        assignment.plate_row = assignment.plate_row or plate_row(assignment.well)
        if assignment.role == 'sample':
            key = sample_axis_key(assignment.well, sample_axis)
            fallback = f'样本{key}' if sample_axis == 'row' else f'样本列{key}'
            assignment.sample = sample_by_axis.get(key, fallback)


def apply_group_to_rows(
    assignments: dict[str, PlateAssignment], wells: list[str], group: str
) -> set[str]:
    return apply_group_to_axis(assignments, wells, group, 'row')


def apply_group_to_selected_wells(
    assignments: dict[str, PlateAssignment], wells: list[str], group: str,
    sample_number: int,
) -> set[str]:
    """Assign one explicit biological sample to exactly the selected wells."""
    if not str(group).strip():
        raise ValueError('Group name must not be empty')
    sample_id = f'{group}-{max(1, int(sample_number))}'
    modified: set[str] = set()
    for well in wells:
        assignment = assignments.setdefault(
            well, PlateAssignment(well=well, plate_row=plate_row(well)))
        if assignment.role != 'sample':
            continue
        assignment.group = group
        assignment.sample = sample_id
        assignment.plate_row = assignment.plate_row or plate_row(well)
        modified.add(well)
    return modified


def apply_group_to_axis(
    assignments: dict[str, PlateAssignment], wells: list[str], group: str,
    sample_axis: str,
) -> set[str]:
    sample_axis = normalize_sample_axis(sample_axis)
    selected_keys = {sample_axis_key(well, sample_axis) for well in wells}
    for assignment in assignments.values():
        key = sample_axis_key(assignment.well, sample_axis)
        if key in selected_keys and assignment.role == 'sample':
            assignment.group = group
            assignment.plate_row = assignment.plate_row or plate_row(assignment.well)
    renumber_samples(assignments, sample_axis)
    return selected_keys


def apply_gene_to_wells(
    assignments: dict[str, PlateAssignment], wells: list[str], gene: str,
    is_reference: bool, sample_axis: str = 'row', group_by_sample: bool = False,
) -> None:
    sample_axis = normalize_sample_axis(sample_axis)
    by_sample_axis: dict[str, list[str]] = {}
    for well in sorted(wells, key=well_sort_key):
        assignment = assignments.get(well)
        key = (assignment.sample if group_by_sample and assignment and assignment.sample
               else sample_axis_key(well, sample_axis))
        by_sample_axis.setdefault(key, []).append(well)
    for axis_wells in by_sample_axis.values():
        position_key = ((lambda well: (well_sort_key(well)[0], well_sort_key(well)[1]))
                        if sample_axis == 'row'
                        else (lambda well: (well_sort_key(well)[1], well_sort_key(well)[0])))
        for replicate, well in enumerate(sorted(axis_wells, key=position_key), 1):
            assignment = assignments.setdefault(well, PlateAssignment(well, plate_row=plate_row(well)))
            if assignment.role != 'sample':
                continue
            assignment.gene = gene
            assignment.is_reference = is_reference
            assignment.technical_replicate = replicate


def save_template(
    path: str | Path, assignments, genes=None, sample_axis: str = 'row',
    mapping_scope: str = 'axis', group_order=None,
) -> None:
    payload = {
        'schema_version': 2,
        'assignments': {k: asdict(v) for k, v in assignments.items()},
        'genes': [asdict(v) for v in (genes or [])],
        'sample_axis': normalize_sample_axis(sample_axis),
        'mapping_scope': normalize_mapping_scope(mapping_scope),
        'group_order': normalize_group_order(group_order),
    }
    Path(path).write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding='utf-8')


def load_template(
    path: str | Path, include_mapping_scope: bool = False,
    include_group_order: bool = False,
):
    payload = json.loads(Path(path).read_text(encoding='utf-8'))
    assignments = {}
    for key, value in payload['assignments'].items():
        value.setdefault('plate_row', plate_row(key))
        assignments[key] = PlateAssignment(**value)
    base = (assignments, payload.get('genes', []), normalize_sample_axis(
        payload.get('sample_axis'), 'row'))
    if include_mapping_scope:
        result = (*base, normalize_mapping_scope(payload.get('mapping_scope'), 'axis'))
        if include_group_order:
            return (*result, normalize_group_order(payload.get('group_order')))
        return result
    if include_group_order:
        return (*base, normalize_group_order(payload.get('group_order')))
    return base
