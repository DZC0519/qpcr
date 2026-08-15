from __future__ import annotations

import json
import shutil
import tempfile
import zipfile
from pathlib import Path

from .mapping import (
    normalize_mapping_scope, normalize_sample_axis, plate_row, renumber_samples,
)
from .models import (
    GROUP_ORDER, AnalysisConfig, AuditAction, ExclusionRecord, GeneDefinition,
    PlateAssignment, PlotSelection, ProjectState, WellRecord, normalize_group_order,
)


def save_project(path: str | Path, state: ProjectState) -> None:
    path = Path(path)
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        (root / 'project.json').write_text(
            json.dumps(state.to_dict(), ensure_ascii=False, indent=2, default=str), encoding='utf-8')
        source = Path(state.source_path)
        if source.exists():
            shutil.copy2(source, root / ('source' + source.suffix.lower()))
        with zipfile.ZipFile(path, 'w', zipfile.ZIP_DEFLATED) as archive:
            for item in root.iterdir():
                archive.write(item, item.name)


def _migrate_group(group: str) -> str:
    aliases = {
        'control':'Ctrl', 'con':'Ctrl', 'ctrl':'Ctrl', '对照':'Ctrl',
        'model':'Model', '模型':'Model', 'low':'Low', '低':'Low',
        'med':'Med', 'medium':'Med', '中':'Med', 'high':'High', '高':'High',
        'pos':'Pos', 'positive':'Pos', '阳性':'Pos',
    }
    return aliases.get(str(group).strip().lower(), group if group in GROUP_ORDER else '')


def load_project(path: str | Path) -> ProjectState:
    path = Path(path)
    with zipfile.ZipFile(path) as archive:
        payload = json.loads(archive.read('project.json').decode('utf-8'))
    version = int(payload.get('schema_version', 1))
    assignments = {}
    for key, value in payload.get('assignments', {}).items():
        value.setdefault('plate_row', plate_row(key))
        if version < 2:
            value['group'] = _migrate_group(value.get('group', ''))
        assignments[key] = PlateAssignment(**value)

    genes_payload = payload.get('genes', [])
    if genes_payload:
        genes = [GeneDefinition(**value) for value in genes_payload]
    else:
        seen = []
        for assignment in assignments.values():
            if assignment.gene and assignment.gene not in seen:
                seen.append(assignment.gene)
        genes = [GeneDefinition(name, any(
            a.gene == name and a.is_reference for a in assignments.values()), index)
            for index, name in enumerate(seen)]

    excluded = payload.get('excluded_wells', {})
    exclusions_payload = payload.get('exclusions', [])
    exclusions = ([ExclusionRecord(**value) for value in exclusions_payload]
                  if exclusions_payload else
                  [ExclusionRecord(well, 'legacy', reason, True) for well, reason in excluded.items()])
    config_payload = payload.get('config', {})
    config_payload.setdefault('calibrator_group', 'Ctrl')
    config_payload.setdefault('sample_axis', 'row')
    config_payload.setdefault('mapping_scope', 'axis')
    config_payload.setdefault('group_order', list(GROUP_ORDER))
    config = AnalysisConfig(**config_payload)
    config.sample_axis = normalize_sample_axis(config.sample_axis, 'row')
    config.mapping_scope = normalize_mapping_scope(config.mapping_scope, 'axis')
    config.group_order = normalize_group_order(config.group_order)
    selections = {key: PlotSelection(**value)
                  for key, value in payload.get('plot_selections', {}).items()}
    audit = []
    for value in payload.get('audit', []):
        value.setdefault('details', {})
        audit.append(AuditAction(**value))
    state = ProjectState(
        source_path=payload.get('source_path', ''), source_sheet=payload.get('source_sheet', ''),
        header_row=payload.get('header_row', 0),
        wells=tuple(WellRecord(**value) for value in payload.get('wells', [])),
        assignments=assignments, genes=genes, config=config,
        excluded_wells=excluded, exclusions=exclusions,
        plot_selections=selections, audit=audit)
    if state.config.mapping_scope == 'axis':
        renumber_samples(state.assignments, state.config.sample_axis)
    state.sync_gene_flags()
    return state
