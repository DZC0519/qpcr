from __future__ import annotations

from collections import defaultdict

import numpy as np

from .mapping import well_sort_key
from .models import AnalysisConfig, PlateAssignment, QCFlag, ReplicateQC, WellRecord


def build_replicate_qc(
    wells: tuple[WellRecord, ...], assignments: dict[str, PlateAssignment],
    config: AnalysisConfig, excluded_wells: dict[str, str] | None = None,
) -> list[ReplicateQC]:
    excluded = excluded_wells or {}
    groups: dict[tuple[str, str, str], list[WellRecord]] = defaultdict(list)
    for record in wells:
        assignment = assignments.get(record.well)
        if not assignment or assignment.role != 'sample' or not assignment.sample or not assignment.gene:
            continue
        groups[(assignment.sample, assignment.group, assignment.gene)].append(record)

    results: list[ReplicateQC] = []
    for (sample, group, gene), records in groups.items():
        records.sort(key=lambda r: well_sort_key(r.well))
        valid = [r for r in records if r.ct is not None and r.include and r.well not in excluded]
        values = np.array([r.ct for r in valid], dtype=float)
        spread = float(np.ptp(values)) if len(values) >= 2 else None
        status = '通过'
        if len(valid) < config.minimum_valid_replicates:
            status = '有效复孔不足'
        elif spread is not None and spread > config.replicate_range_threshold:
            status = '偏离过大'

        suggested_well, suggested_ct, after = '', None, None
        if status == '偏离过大' and len(valid) >= 3:
            median = float(np.median(values))
            distances = np.abs(values - median)
            max_distance = float(distances.max())
            candidates = np.flatnonzero(np.isclose(distances, max_distance, rtol=1e-9, atol=1e-12))
            if len(candidates) == 1:
                candidate = int(candidates[0])
                remaining = np.delete(values, candidate)
                after = float(np.ptp(remaining)) if len(remaining) >= 2 else 0.0
                if (len(remaining) >= config.minimum_valid_replicates
                        and after <= config.replicate_range_threshold):
                    suggested_well = valid[candidate].well
                    suggested_ct = float(valid[candidate].ct)

        results.append(ReplicateQC(
            sample=sample, group=group, gene=gene,
            wells=[r.well for r in records], raw_ct=[r.ct for r in records],
            valid_wells=[r.well for r in valid],
            excluded_wells=[r.well for r in records if r.well in excluded],
            n_valid=len(valid), mean_ct=float(values.mean()) if len(valid) >= config.minimum_valid_replicates else None,
            range_ct=spread, status=status, suggested_well=suggested_well,
            suggested_ct=suggested_ct, range_after_suggestion=after,
        ))
    return sorted(results, key=lambda q: (well_sort_key(q.wells[0]) if q.wells else (999, 999), q.gene))


def run_qc(
    wells: tuple[WellRecord, ...], assignments: dict[str, PlateAssignment],
    config: AnalysisConfig, excluded_wells: dict[str, str] | None = None,
) -> list[QCFlag]:
    flags: list[QCFlag] = []
    for item in build_replicate_qc(wells, assignments, config, excluded_wells):
        if item.status == '有效复孔不足':
            flags.append(QCFlag('insufficient_replicates', 'error',
                f'{item.sample}/{item.gene} 仅有 {item.n_valid} 个有效技术孔',
                item.wells, item.sample, item.gene, float(item.n_valid)))
        elif item.status == '偏离过大':
            message = f'{item.sample}/{item.gene} 复孔极差 {item.range_ct:.2f} Cp'
            if item.suggested_well:
                message += f'；建议排除 {item.suggested_well}'
            flags.append(QCFlag('replicate_spread', 'warning', message,
                item.valid_wells, item.sample, item.gene, item.range_ct))
    return flags
