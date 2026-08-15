from __future__ import annotations

from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Literal

Role = Literal['sample', 'empty', 'ntc', 'no_rt', 'standard']
GROUP_ORDER = ('Ctrl', 'Model', 'Low', 'Med', 'High', 'Pos')


def normalize_group_order(value: Any = None) -> list[str]:
    """Return a stable, non-empty group list while preserving user order."""
    source = list(GROUP_ORDER) if value is None else value
    if isinstance(source, str):
        source = [source]
    result: list[str] = []
    seen: set[str] = set()
    for item in source or ():
        name = str(item).strip()
        key = name.casefold()
        if name and key not in seen:
            result.append(name)
            seen.add(key)
    return result or list(GROUP_ORDER)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec='seconds')


@dataclass(frozen=True, slots=True)
class WellRecord:
    well: str
    ct: float | None
    include: bool = True
    color: Any = None
    machine_name: str = ''
    concentration: float | None = None
    standard: Any = None
    status: str = ''
    source_sheet: str = ''
    source_row: int = 0
    raw: dict[str, Any] = field(default_factory=dict)


@dataclass(slots=True)
class GeneDefinition:
    name: str
    is_reference: bool = False
    order: int = 0


@dataclass(slots=True)
class PlateAssignment:
    well: str
    sample: str = ''
    group: str = ''
    gene: str = ''
    is_reference: bool = False
    technical_replicate: int | None = None
    role: Role = 'sample'
    plate_row: str = ''


@dataclass(slots=True)
class AnalysisConfig:
    plate_format: int = 384
    replicate_range_threshold: float = 0.5
    minimum_valid_replicates: int = 2
    minimum_calibrator_samples: int = 2
    amplification_efficiency: float = 1.0
    calibrator_group: str = 'Ctrl'
    alpha: float = 0.05
    error_bar: str = 'SEM'
    sample_axis: Literal['auto', 'row', 'column'] = 'auto'
    mapping_scope: Literal['axis', 'well'] = 'axis'
    group_order: list[str] = field(default_factory=lambda: list(GROUP_ORDER))

    def __post_init__(self) -> None:
        self.group_order = normalize_group_order(self.group_order)


@dataclass(slots=True)
class ExclusionRecord:
    well: str
    source: Literal['qc_suggestion', 'manual', 'legacy'] = 'manual'
    reason: str = ''
    active: bool = True
    timestamp: str = field(default_factory=utc_now)


@dataclass(slots=True)
class PlotSelection:
    gene: str
    samples: list[str] = field(default_factory=list)
    exploratory: bool = False


@dataclass(slots=True)
class QCFlag:
    code: str
    severity: Literal['info', 'warning', 'error']
    message: str
    wells: list[str] = field(default_factory=list)
    sample: str = ''
    gene: str = ''
    value: float | None = None


@dataclass(slots=True)
class ReplicateQC:
    sample: str
    group: str
    gene: str
    wells: list[str] = field(default_factory=list)
    raw_ct: list[float | None] = field(default_factory=list)
    valid_wells: list[str] = field(default_factory=list)
    excluded_wells: list[str] = field(default_factory=list)
    n_valid: int = 0
    mean_ct: float | None = None
    range_ct: float | None = None
    status: str = ''
    suggested_well: str = ''
    suggested_ct: float | None = None
    range_after_suggestion: float | None = None


@dataclass(slots=True)
class AuditAction:
    action: str
    wells: list[str] = field(default_factory=list)
    reason: str = ''
    details: dict[str, Any] = field(default_factory=dict)
    timestamp: str = field(default_factory=utc_now)


@dataclass(slots=True)
class AnalysisResult:
    replicate_summary: list[dict[str, Any]] = field(default_factory=list)
    delta_ct: list[dict[str, Any]] = field(default_factory=list)
    relative_expression: list[dict[str, Any]] = field(default_factory=list)
    statistics: list[dict[str, Any]] = field(default_factory=list)
    qc_flags: list[QCFlag] = field(default_factory=list)
    replicate_qc: list[ReplicateQC] = field(default_factory=list)
    gene_tables: dict[str, list[dict[str, Any]]] = field(default_factory=dict)
    calibrator_counts: dict[str, int] = field(default_factory=dict)
    exploratory_genes: set[str] = field(default_factory=set)


@dataclass(slots=True)
class ProjectState:
    source_path: str = ''
    source_sheet: str = ''
    header_row: int = 0
    wells: tuple[WellRecord, ...] = ()
    assignments: dict[str, PlateAssignment] = field(default_factory=dict)
    config: AnalysisConfig = field(default_factory=AnalysisConfig)
    genes: list[GeneDefinition] = field(default_factory=list)
    excluded_wells: dict[str, str] = field(default_factory=dict)
    exclusions: list[ExclusionRecord] = field(default_factory=list)
    plot_selections: dict[str, PlotSelection] = field(default_factory=dict)
    audit: list[AuditAction] = field(default_factory=list)

    def exclude(self, wells: list[str], reason: str = '', source: str = 'manual') -> None:
        for well in wells:
            self.excluded_wells[well] = reason or ('QC建议排除' if source == 'qc_suggestion' else '手动排除')
            for record in reversed(self.exclusions):
                if record.well == well and record.active:
                    record.active = False
            self.exclusions.append(ExclusionRecord(well, source, reason, True))
        self.audit.append(AuditAction('exclude', list(wells), reason, {'source': source}))

    def restore(self, wells: list[str], reason: str = '') -> None:
        for well in wells:
            self.excluded_wells.pop(well, None)
            for record in reversed(self.exclusions):
                if record.well == well and record.active:
                    record.active = False
                    break
        self.audit.append(AuditAction('restore', list(wells), reason))

    def set_plot_selection(self, gene: str, samples: list[str], all_samples: list[str]) -> None:
        chosen = sorted(set(samples))
        exploratory = set(chosen) != set(all_samples)
        self.plot_selections[gene] = PlotSelection(gene, chosen, exploratory)
        self.audit.append(AuditAction(
            'plot_filter', reason='探索性筛选' if exploratory else '恢复全部样本',
            details={'gene': gene, 'samples': chosen, 'exploratory': exploratory}))

    def reference_genes(self) -> list[str]:
        if self.genes:
            return [g.name for g in sorted(self.genes, key=lambda x: x.order) if g.is_reference]
        seen: list[str] = []
        for assignment in self.assignments.values():
            if assignment.gene and assignment.is_reference and assignment.gene not in seen:
                seen.append(assignment.gene)
        return seen

    def target_genes(self) -> list[str]:
        if self.genes:
            return [g.name for g in sorted(self.genes, key=lambda x: x.order) if not g.is_reference]
        seen: list[str] = []
        for assignment in self.assignments.values():
            if assignment.gene and not assignment.is_reference and assignment.gene not in seen:
                seen.append(assignment.gene)
        return seen

    def sync_gene_flags(self) -> None:
        refs = set(self.reference_genes())
        for assignment in self.assignments.values():
            assignment.is_reference = assignment.gene in refs

    def to_dict(self, include_wells: bool = True) -> dict[str, Any]:
        payload = {
            'schema_version': 2, 'app_version': '2.1.0',
            'source_path': self.source_path, 'source_sheet': self.source_sheet,
            'header_row': self.header_row,
            'assignments': {k: asdict(v) for k, v in self.assignments.items()},
            'genes': [asdict(v) for v in self.genes],
            'config': asdict(self.config),
            'excluded_wells': dict(self.excluded_wells),
            'exclusions': [asdict(v) for v in self.exclusions],
            'plot_selections': {k: asdict(v) for k, v in self.plot_selections.items()},
            'audit': [asdict(v) for v in self.audit],
        }
        if include_wells:
            payload['wells'] = [asdict(v) for v in self.wells]
        return payload

    @property
    def source_name(self) -> str:
        return Path(self.source_path).name if self.source_path else ''
