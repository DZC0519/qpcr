from PySide6.QtCore import Qt
from PySide6.QtWidgets import QAbstractItemView, QApplication, QPushButton, QTabWidget

from qpcr_analyzer.analysis import analyze
from qpcr_analyzer.exports import _reference_significance, make_gene_figure
from qpcr_analyzer.main import MainWindow
from qpcr_analyzer.models import (
    GROUP_ORDER, AnalysisConfig, GeneDefinition, PlateAssignment, ProjectState, WellRecord,
)
from qpcr_analyzer.widgets import DropLabel, GeneResultPage, PlateTable


def analysis_state():
    wells, assignments = [], {}
    layout = [('A', 'Ctrl-1', 'Ctrl', 20.0), ('B', 'Ctrl-2', 'Ctrl', 20.2),
              ('C', 'Model-1', 'Model', 18.0), ('D', 'Model-2', 'Model', 18.2)]
    for plate_row, sample, group, target in layout:
        for index, ct in enumerate([target, target + .1, target - .1], 1):
            well = f'{plate_row}{index}'; wells.append(WellRecord(well, ct))
            assignments[well] = PlateAssignment(
                well, sample, group, 'Target', False, index, plate_row=plate_row)
        for index, ct in enumerate([17.0, 17.1, 16.9], 4):
            well = f'{plate_row}{index}'; wells.append(WellRecord(well, ct))
            assignments[well] = PlateAssignment(
                well, sample, group, 'REF', True, index - 3, plate_row=plate_row)
    return ProjectState(
        wells=tuple(wells), assignments=assignments, config=AnalysisConfig(96),
        genes=[GeneDefinition('Target', False, 0), GeneDefinition('REF', True, 1)],
    )


def test_v2_ui_fixed_groups_and_bottom_tabs(qtbot):
    window = MainWindow(); qtbot.addWidget(window)
    assert window.stack.count() == 5
    assert tuple(window.group_combo.itemText(index) for index in range(window.group_combo.count())) == GROUP_ORDER
    assert window.group_combo.isEditable() is False
    assert [window.group_table.item(row, 0).text()
            for row in range(window.group_table.rowCount())] == list(GROUP_ORDER)
    assert window.gene_tabs.tabPosition() == QTabWidget.TabPosition.South
    button_texts = {button.text() for button in window.findChildren(QPushButton)}
    assert '一键导出全部靶基因 SVG' not in button_texts
    assert '导出清理后计算表 Excel' in button_texts


def test_group_names_are_single_click_editable_and_rename_assignments(qtbot):
    window = MainWindow(); qtbot.addWidget(window)
    window.state.wells = (WellRecord('A1', 20.0),)
    window.state.assignments = {
        'A1': PlateAssignment('A1', 'Model-1', 'Model', 'Target', plate_row='A')}
    item = window.group_table.item(1, 0)
    qtbot.mouseClick(
        window.group_table.viewport(), Qt.MouseButton.LeftButton,
        pos=window.group_table.visualItemRect(item).center(),
    )
    qtbot.waitUntil(
        lambda: window.group_table.state() == QAbstractItemView.State.EditingState)
    item.setText('Disease')
    assert window.commit_group_definitions()
    assert window.state.config.group_order[1] == 'Disease'
    assert window.group_combo.itemText(1) == 'Disease'
    assert window.state.assignments['A1'].group == 'Disease'
    assert window.state.assignments['A1'].sample == 'Disease-1'


def test_plate_selection_is_described_and_group_applies_only_to_selected_wells(qtbot):
    window = MainWindow(); qtbot.addWidget(window)
    window.state.config.sample_axis = 'auto'
    window._sync_sample_axis_ui()
    window.state.wells = tuple(WellRecord(f'{row}{column}', 20.0 + column/10)
                               for row in ('A', 'B') for column in range(1, 4))
    window.state.assignments = {record.well: PlateAssignment(
        record.well, sample=f'样本{record.well[0]}', plate_row=record.well[0])
        for record in window.state.wells}
    window._refresh_plate()
    window.plate_table.item(0, 0).setSelected(True)
    window.plate_table.item(0, 1).setSelected(True)
    qtbot.waitUntil(lambda: '已选择 2 孔' in window.selection_label.text())
    window.group_combo.setCurrentText('Model'); window.sample_number.setValue(2)
    window.apply_group()
    assert {window.state.assignments[f'A{column}'].sample for column in range(1, 3)} == {'Model-2'}
    assert window.state.assignments['A3'].group == ''
    assert all(window.state.assignments[f'B{column}'].group == '' for column in range(1, 4))
    assert window.state.config.mapping_scope == 'well'


def test_legacy_axis_group_button_still_expands_selected_row(qtbot):
    window = MainWindow(); qtbot.addWidget(window)
    window.state.wells = tuple(
        WellRecord(f'{row}{column}', 20.0) for row in ('A', 'B') for column in range(1, 4))
    window.state.assignments = {
        record.well: PlateAssignment(record.well, plate_row=record.well[0])
        for record in window.state.wells
    }
    window.state.config.sample_axis = 'auto'
    window._sync_sample_axis_ui(); window._refresh_plate()
    window.plate_table.item(0, 0).setSelected(True)
    window.plate_table.item(0, 1).setSelected(True)
    window.group_combo.setCurrentText('Ctrl'); window.apply_group_by_axis()
    assert all(window.state.assignments[f'A{column}'].group == 'Ctrl' for column in range(1, 4))
    assert all(window.state.assignments[f'B{column}'].group == '' for column in range(1, 4))
    assert window.state.config.mapping_scope == 'axis'


def test_group_selection_auto_detects_columns_without_touching_other_columns(qtbot):
    window = MainWindow(); qtbot.addWidget(window)
    window.state.wells = tuple(
        WellRecord(f'{row}{column}', 20.0 + column / 10)
        for row in ('A', 'B', 'C') for column in range(1, 5)
    )
    window.state.assignments = {
        record.well: PlateAssignment(record.well, plate_row=record.well[0])
        for record in window.state.wells
    }
    window.state.config.sample_axis = 'auto'
    window._sync_sample_axis_ui(); window._refresh_plate()
    for row in range(3):
        window.plate_table.item(row, 0).setSelected(True)
    window.group_combo.setCurrentText('Ctrl'); window.apply_group()

    assert window.state.config.sample_axis == 'column'
    assert '复孔纵向' in window.sample_axis_status.text()
    for row in range(3):
        for column in (1, 2):
            window.plate_table.item(row, column).setSelected(True)
    window.apply_group()
    for row in ('A', 'B', 'C'):
        assert all(window.state.assignments[f'{row}{column}'].group == 'Ctrl'
                   for column in range(1, 4))
        assert window.state.assignments[f'{row}4'].group == ''
    assert {window.state.assignments[f'{row}1'].sample for row in ('A', 'B', 'C')} == {'Ctrl-1'}


def test_cross_row_ctrl_selection_assigns_exact_blocks_and_gene_replicates(qtbot):
    window = MainWindow(); qtbot.addWidget(window)
    window.state.wells = tuple(
        WellRecord(f'{row}{column}', 20.0)
        for row in ('A', 'B', 'C', 'D') for column in range(9, 14)
    )
    window.state.assignments = {
        record.well: PlateAssignment(record.well, plate_row=record.well[0])
        for record in window.state.wells
    }
    window.state.genes = [GeneDefinition('Target', False, 0)]
    window.state.config.sample_axis = 'auto'
    window._load_gene_table_from_state(); window._refresh_gene_combo(); window._refresh_plate()
    for row in (1, 3):
        for column in (9, 10, 11):
            window.plate_table.item(row, column).setSelected(True)
    window.group_combo.setCurrentText('Low'); window.sample_number.setValue(3)
    window.apply_group()
    selected = [f'{row}{column}' for row in ('B', 'D') for column in range(10, 13)]
    assert {window.state.assignments[well].sample for well in selected} == {'Low-3'}
    assert window.state.assignments['A10'].group == ''
    assert window.state.config.sample_axis == 'row'
    window.gene_combo.setCurrentText('Target'); window.apply_gene()
    assert [window.state.assignments[f'{row}{column}'].technical_replicate
            for row in ('B', 'D') for column in range(10, 13)] == list(range(1, 7))


def test_undo_restores_mapping_scope(qtbot):
    window = MainWindow(); qtbot.addWidget(window)
    window.state.wells = tuple(WellRecord(f'A{column}', 20.0) for column in range(1, 4))
    window.state.assignments = {
        record.well: PlateAssignment(record.well, plate_row='A') for record in window.state.wells
    }
    window.state.config.sample_axis = 'auto'; window._refresh_plate()
    window.plate_table.item(0, 0).setSelected(True)
    window.plate_table.item(0, 1).setSelected(True)
    window.apply_group()
    assert window.state.config.mapping_scope == 'well'
    window.undo()
    assert window.state.config.mapping_scope == 'axis'
    assert window.state.config.sample_axis == 'auto'


def test_gene_selection_auto_detects_vertical_replicates(qtbot):
    window = MainWindow(); qtbot.addWidget(window)
    window.state.wells = tuple(WellRecord(f'{row}1', 20.0) for row in ('A', 'B', 'C'))
    window.state.assignments = {
        record.well: PlateAssignment(record.well, plate_row=record.well[0])
        for record in window.state.wells
    }
    window.state.genes = [GeneDefinition('Target', False, 0)]
    window.state.config.sample_axis = 'auto'
    window._refresh_gene_combo(); window._sync_sample_axis_ui(); window._refresh_plate()
    for row in range(3):
        window.plate_table.item(row, 0).setSelected(True)
    window.gene_combo.setCurrentText('Target'); window.apply_gene()

    assert window.state.config.sample_axis == 'column'
    assert [window.state.assignments[f'{row}1'].technical_replicate
            for row in ('A', 'B', 'C')] == [1, 2, 3]


def test_drop_zone_is_clickable(qtbot):
    drop = DropLabel(); qtbot.addWidget(drop); drop.show()
    clicks = []
    drop.clicked.connect(lambda: clicks.append(True))
    qtbot.mouseClick(drop, Qt.MouseButton.LeftButton)
    assert clicks == [True]


def test_384_plate_keeps_all_rows_columns_and_exact_selection(qtbot):
    table = PlateTable(); qtbot.addWidget(table)
    rows = tuple(chr(ord('A') + index) for index in range(16))
    wells = tuple(
        WellRecord(f'{row}{column}', 20.0)
        for row in rows for column in range(1, 25)
    )
    assignments = {
        record.well: PlateAssignment(record.well, plate_row=record.well[0])
        for record in wells
    }
    table.set_plate(wells, assignments)
    selected = [f'{row}{column}' for row in ('B', 'D', 'F', 'H') for column in range(10, 13)]
    table.select_wells(selected)
    assert table.rowCount() == 16
    assert table.columnCount() == 24
    assert table.selected_wells() == selected


def test_gene_name_uses_single_click_editing_and_new_row_focus(qtbot):
    window = MainWindow(); qtbot.addWidget(window); window.show()
    window.add_gene_row('Target')
    qtbot.mouseClick(
        window.gene_table.viewport(), Qt.MouseButton.LeftButton,
        pos=window.gene_table.visualItemRect(window.gene_table.item(0, 0)).center(),
    )
    qtbot.waitUntil(
        lambda: window.gene_table.state() == QAbstractItemView.State.EditingState)

    window.gene_table.closePersistentEditor(window.gene_table.item(0, 0))
    window.add_gene_row(start_edit=True)
    assert window.gene_table.currentRow() == 1
    qtbot.waitUntil(
        lambda: window.gene_table.state() == QAbstractItemView.State.EditingState)


def test_analysis_table_can_exclude_and_restore_selected_sample(qtbot):
    state = analysis_state()
    result = analyze(state.wells, state.assignments, state.config)
    page = GeneResultPage('Target'); qtbot.addWidget(page)
    page.update_data(state, result)
    changes = []
    page.filterChanged.connect(lambda gene, samples: changes.append((gene, samples)))

    page.table.selectRow(3)
    page._exclude_selected_rows()
    assert changes[-1][0] == 'Target'
    assert 'Model-2' not in changes[-1][1]

    filtered_samples = changes[-1][1]
    state.set_plot_selection('Target', filtered_samples, page._all_samples)
    filtered = analyze(
        state.wells, state.assignments, state.config,
        sample_filters={'Target': set(filtered_samples)},
    )
    page.update_data(state, filtered)
    page.table.selectRow(3)
    page._restore_selected_rows()
    assert 'Model-2' in changes[-1][1]


def test_analysis_page_prioritizes_table_and_collapses_chart(qtbot):
    state = analysis_state(); result = analyze(state.wells, state.assignments, state.config)
    page = GeneResultPage('Target'); qtbot.addWidget(page); page.update_data(state, result)
    assert not hasattr(page, 'filter_tree')
    assert page.table.minimumHeight() >= 320
    assert page.table.maximumHeight() > 10000
    assert page.canvas.isHidden()
    page.chart_toggle.setChecked(True)
    assert not page.canvas.isHidden()
    assert page.canvas.minimumWidth() >= 480
    assert page.canvas.minimumHeight() >= 460


def test_gene_figure_marks_significant_model_reference_comparison():
    state = analysis_state(); result = analyze(state.wells, state.assignments, state.config)
    figure = make_gene_figure(result, 'Target', 'SEM')
    labels = [text.get_text() for text in figure.axes[0].texts]
    assert '**' in labels


def test_figure_stars_use_raw_p_and_model_as_reference():
    state = analysis_state(); result = analyze(state.wells, state.assignments, state.config)
    result.statistics = [
        {'gene': 'Target', 'group_a': 'Ctrl', 'group_b': 'Model',
         'p_value': .04, 'p_adjust_bh': .4},
        {'gene': 'Target', 'group_a': 'Ctrl', 'group_b': 'Low',
         'p_value': .001, 'p_adjust_bh': .001},
        {'gene': 'Target', 'group_a': 'Model', 'group_b': 'Low',
         'p_value': .009, 'p_adjust_bh': .5},
    ]
    reference, comparisons = _reference_significance(
        result, 'Target', ['Ctrl', 'Model', 'Low'])
    assert reference == 'Model'
    assert comparisons == [('Ctrl', '*'), ('Low', '**')]


def test_chart_controls_explain_model_reference_and_raw_p(qtbot):
    state = analysis_state(); result = analyze(state.wells, state.assignments, state.config)
    page = GeneResultPage('Target'); qtbot.addWidget(page)
    page.update_data(state, result)
    assert '原始 P' in page.significance_hint.text()
    assert 'Model' in page.significance_hint.text()
    assert page.stats_table.columnCount() == 8


def test_copy_chart_places_image_on_clipboard(qtbot):
    state = analysis_state()
    result = analyze(state.wells, state.assignments, state.config)
    page = GeneResultPage('Target'); qtbot.addWidget(page)
    page.update_data(state, result)
    button = QPushButton('复制图片')
    page._copy_figure(button)
    assert not QApplication.clipboard().image().isNull()
