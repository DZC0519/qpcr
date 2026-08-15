from pathlib import Path
import math

import pytest
from openpyxl import Workbook

from qpcr_analyzer.analysis import analyze, bh_adjust
from qpcr_analyzer.importer import find_candidates, load_plate
from qpcr_analyzer.mapping import (
    apply_gene_to_wells, apply_group_to_axis, apply_group_to_selected_wells,
    detect_sample_axis, load_template, plate_row, save_template,
    suggest_triplicate_layout,
)
from qpcr_analyzer.models import AnalysisConfig, PlateAssignment, WellRecord
from qpcr_analyzer.qc import build_replicate_qc, run_qc


def make_book(path: Path):
    wb = Workbook(); ws = wb.active; ws.title = 'Run1'
    ws.append(['machine export']); ws.append(['Include', 'Pos', 'Name', 'Cp', 'Status'])
    values = {'A1':20.0, 'A2':20.1, 'A3':20.8, 'A4':18.0,
              'A5':18.1, 'A6':18.2, 'B1':None}
    for well, ct in values.items():
        ws.append([True, well, 'sample', ct, ''])
    ws.append([True, 'H12', 'blank', 22.0, '? - uncertain']); wb.save(path)


def test_import_and_replicate_only_qc(tmp_path):
    path = tmp_path/'input.xlsx'; make_book(path)
    assert find_candidates(path)[0].header_row == 2
    plate = load_plate(path)
    assert plate.plate_format == 96 and len(plate.wells) == 8
    assignments = suggest_triplicate_layout(plate.wells, empty_rows=('H',))
    flags = run_qc(plate.wells, assignments, AnalysisConfig())
    spread = next(flag for flag in flags if flag.code == 'replicate_spread')
    assert spread.wells == ['A1', 'A2', 'A3']
    assert not any(flag.code == 'control_amplification' for flag in flags)
    qc = build_replicate_qc(plate.wells, assignments, AnalysisConfig())
    assert next(row for row in qc if row.gene == '基因1').suggested_well == 'A3'
    assert next(well for well in plate.wells if well.well == 'B1').ct is None


@pytest.mark.parametrize('values, expected', [
    ([20.0, 20.8], ''),
    ([20.0, 20.1, 21.0], 'A3'),
    ([20.0, 20.5, 21.0], ''),
    ([20.0, 20.1, 20.2, 21.0], 'A4'),
    ([20.0, 20.1, 20.2, 20.15, 21.0], 'A5'),
])
def test_dynamic_replicate_outlier_suggestion(values, expected):
    wells = tuple(WellRecord(f'A{index+1}', value) for index, value in enumerate(values))
    assignments = {well.well: PlateAssignment(
        well.well, 'Ctrl-1', 'Ctrl', 'Target', False, index+1, plate_row='A')
        for index, well in enumerate(wells)}
    result = build_replicate_qc(wells, assignments, AnalysisConfig())[0]
    assert result.suggested_well == expected


def _analysis_data():
    wells, assignments = [], {}
    rows = [('A', 'Ctrl-1', 'Ctrl', 20.0), ('B', 'Ctrl-2', 'Ctrl', 20.2),
            ('C', 'Model-1', 'Model', 18.0), ('D', 'Model-2', 'Model', 18.2)]
    for plate_row, sample, group, target in rows:
        for index, ct in enumerate([target, target+.1, target-.1], 1):
            well = plate_row+str(index); wells.append(WellRecord(well, ct))
            assignments[well] = PlateAssignment(well, sample, group, 'Target', False, index, plate_row=plate_row)
        for index, ct in enumerate([17.0, 17.1, 16.9], 4):
            well = plate_row+str(index); wells.append(WellRecord(well, ct))
            assignments[well] = PlateAssignment(well, sample, group, 'REF', True, index-3, plate_row=plate_row)
    return tuple(wells), assignments


def test_delta_ct_fold_change_statistics_and_filter_recalculation():
    wells, assignments = _analysis_data()
    result = analyze(wells, assignments, AnalysisConfig())
    samples = [row for row in result.relative_expression if row.get('sample') != '组汇总']
    ctrl = [row for row in samples if row['group'] == 'Ctrl']
    model = [row for row in samples if row['group'] == 'Model']
    assert math.isclose(sum(row['fold_change'] for row in ctrl)/2, 1.0024, rel_tol=.02)
    assert all(row['fold_change'] > 3.5 for row in model)
    assert result.statistics[0]['test'] == 'Welch t-test'

    filtered = analyze(wells, assignments, AnalysisConfig(), sample_filters={
        'Target': {'Ctrl-1', 'Ctrl-2', 'Model-1'}})
    assert 'Target' in filtered.exploratory_genes
    assert not filtered.statistics
    assert next(row for row in filtered.gene_tables['Target'] if row['sample'] == 'Model-2')[
        'included_for_analysis'] is False


def test_custom_group_order_drives_analysis_and_statistics():
    wells, assignments = _analysis_data()
    for assignment in assignments.values():
        if assignment.group == 'Ctrl':
            assignment.group = 'Control'
            assignment.sample = assignment.sample.replace('Ctrl-', 'Control-')
        elif assignment.group == 'Model':
            assignment.group = 'Drug'
            assignment.sample = assignment.sample.replace('Model-', 'Drug-')
    config = AnalysisConfig(
        group_order=['Control', 'Drug'], calibrator_group='Control')
    result = analyze(wells, assignments, config)
    assert [row['group'] for row in result.gene_tables['Target']] == [
        'Control', 'Control', 'Drug', 'Drug']
    assert result.statistics[0]['comparison'] == 'Control vs Drug'
    assert result.calibrator_counts['Target'] == 2


def test_missing_reference_and_ctrl_are_not_computable():
    wells = (WellRecord('A1', 20), WellRecord('A2', 20.1), WellRecord('A3', 19.9))
    assignments = {well.well: PlateAssignment(
        well.well, 'Ctrl-1', 'Ctrl', 'Target', False, index+1, plate_row='A')
        for index, well in enumerate(wells)}
    result = analyze(wells, assignments, AnalysisConfig())
    assert result.delta_ct[0]['computable'] is False
    assert math.isnan(result.delta_ct[0]['delta_ct'])
    assert math.isnan(result.relative_expression[0]['fold_change'])


def test_bh_adjust_monotonic():
    adjusted = bh_adjust([.01, .04, .03, .2])
    assert all(0 <= value <= 1 for value in adjusted)
    assert adjusted[0] <= adjusted[2] <= adjusted[3]


def test_default_error_bar_is_sem():
    assert AnalysisConfig().error_bar == 'SEM'
    assert AnalysisConfig().sample_axis == 'auto'
    assert AnalysisConfig().mapping_scope == 'axis'


def test_sample_axis_is_inferred_from_action_and_selection_direction():
    vertical = ['A1', 'B1', 'C1']
    horizontal = ['A1', 'A2', 'A3']
    square = ['A1', 'A2', 'B1', 'B2']
    assert detect_sample_axis(vertical, 'group', 16, 24) == 'column'
    assert detect_sample_axis(vertical, 'gene', 16, 24) == 'column'
    assert detect_sample_axis(horizontal, 'group', 16, 24) == 'row'
    assert detect_sample_axis(horizontal, 'gene', 16, 24) == 'row'
    assert detect_sample_axis(square, 'group', 16, 24) is None


def test_discontiguous_cross_row_blocks_use_contiguous_run_direction():
    horizontal_blocks = [
        f'{row}{column}' for row in ('B', 'D', 'F', 'H') for column in range(10, 13)
    ]
    vertical_blocks = [
        f'{row}{column}' for column in (10, 12, 14) for row in ('B', 'C', 'D')
    ]
    assert detect_sample_axis(horizontal_blocks, 'gene', 16, 24) == 'row'
    assert detect_sample_axis(horizontal_blocks, 'group', 16, 24) == 'row'
    assert detect_sample_axis(vertical_blocks, 'gene', 16, 24) == 'column'


def test_exact_group_assignment_changes_only_selected_wells():
    wells = [f'{row}{column}' for row in ('A', 'B', 'C') for column in range(1, 5)]
    assignments = {well: PlateAssignment(well, plate_row=plate_row(well)) for well in wells}
    selected = ['A1', 'A2', 'C1', 'C2']
    modified = apply_group_to_selected_wells(assignments, selected, 'Model', 4)
    assert modified == set(selected)
    assert {assignments[well].sample for well in selected} == {'Model-4'}
    assert all(assignments[well].group == 'Model' for well in selected)
    assert all(assignments[well].group == '' for well in set(wells) - set(selected))


def test_exact_gene_assignment_numbers_replicates_per_explicit_sample():
    selected = [f'{row}{column}' for row in ('B', 'D') for column in range(10, 13)]
    assignments = {
        well: PlateAssignment(
            well, 'Ctrl-1' if well.startswith('B') else 'Ctrl-2', 'Ctrl',
            plate_row=plate_row(well),
        )
        for well in selected
    }
    apply_gene_to_wells(
        assignments, selected, 'Target', False, 'row', group_by_sample=True)
    assert [assignments[f'B{column}'].technical_replicate for column in range(10, 13)] == [1, 2, 3]
    assert [assignments[f'D{column}'].technical_replicate for column in range(10, 13)] == [1, 2, 3]


def test_column_sample_mapping_expands_columns_and_numbers_vertical_replicates(tmp_path):
    wells = [f'{row}{column}' for row in ('A', 'B', 'C') for column in range(1, 4)]
    assignments = {
        well: PlateAssignment(well, plate_row=well[0]) for well in wells
    }
    targets = apply_group_to_axis(assignments, ['A1', 'A2'], 'Ctrl', 'column')
    assert targets == {'1', '2'}
    assert all(assignments[f'{row}1'].group == 'Ctrl' for row in ('A', 'B', 'C'))
    assert all(assignments[f'{row}2'].group == 'Ctrl' for row in ('A', 'B', 'C'))
    assert all(assignments[f'{row}3'].group == '' for row in ('A', 'B', 'C'))
    assert {assignments[f'{row}1'].sample for row in ('A', 'B', 'C')} == {'Ctrl-1'}
    assert {assignments[f'{row}2'].sample for row in ('A', 'B', 'C')} == {'Ctrl-2'}

    apply_gene_to_wells(assignments, ['A1', 'B1', 'C1'], 'Target', False, 'column')
    assert [assignments[f'{row}1'].technical_replicate for row in ('A', 'B', 'C')] == [1, 2, 3]

    template = tmp_path / 'column-template.json'
    save_template(template, assignments, sample_axis='column')
    restored, _, sample_axis = load_template(template)
    assert sample_axis == 'column'
    assert restored['C1'].technical_replicate == 3


def test_exact_mapping_scope_survives_template_round_trip(tmp_path):
    assignments = {
        well: PlateAssignment(well, 'Low-2', 'Low', 'Target', False, index, plate_row='B')
        for index, well in enumerate(('B10', 'B11', 'B12'), 1)
    }
    template = tmp_path / 'exact-template.json'
    save_template(template, assignments, sample_axis='row', mapping_scope='well')
    restored, _, sample_axis, mapping_scope = load_template(
        template, include_mapping_scope=True)
    assert sample_axis == 'row'
    assert mapping_scope == 'well'
    assert [restored[well].sample for well in ('B10', 'B11', 'B12')] == ['Low-2'] * 3


def test_custom_group_order_survives_template_round_trip(tmp_path):
    assignments = {'A1': PlateAssignment('A1', 'Control-1', 'Control', 'Target')}
    template = tmp_path / 'custom-groups-template.json'
    save_template(
        template, assignments, sample_axis='row', mapping_scope='well',
        group_order=['Control', 'Drug'])
    restored, _, _, _, groups = load_template(
        template, include_mapping_scope=True, include_group_order=True)
    assert restored['A1'].group == 'Control'
    assert groups == ['Control', 'Drug']


def test_column_layout_runs_complete_relative_expression_analysis():
    wells = tuple(
        WellRecord(f'{row}{column}', ct)
        for column, target in ((1, 20.0), (2, 20.2), (3, 18.0), (4, 18.2))
        for row, ct in zip(
            ('A', 'B', 'C', 'D', 'E', 'F'),
            (target, target + .1, target - .1, 17.0, 17.1, 16.9),
        )
    )
    assignments = {
        record.well: PlateAssignment(record.well, plate_row=plate_row(record.well))
        for record in wells
    }
    apply_group_to_axis(assignments, ['A1', 'A2'], 'Ctrl', 'column')
    apply_group_to_axis(assignments, ['A3', 'A4'], 'Model', 'column')
    apply_gene_to_wells(
        assignments,
        [f'{row}{column}' for row in ('A', 'B', 'C') for column in range(1, 5)],
        'Target', False, 'column')
    apply_gene_to_wells(
        assignments,
        [f'{row}{column}' for row in ('D', 'E', 'F') for column in range(1, 5)],
        'REF', True, 'column')

    result = analyze(wells, assignments, AnalysisConfig(384, sample_axis='column'))
    assert [row['sample'] for row in result.gene_tables['Target']] == [
        'Ctrl-1', 'Ctrl-2', 'Model-1', 'Model-2']
    assert result.statistics[0]['test'] == 'Welch t-test'
    assert len(result.replicate_qc) == 8
