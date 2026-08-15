from __future__ import annotations

from collections import defaultdict
from io import BytesIO
from pathlib import Path

import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from PySide6.QtCore import QItemSelection, QItemSelectionModel, QTimer, Qt, Signal
from PySide6.QtGui import QColor, QFont, QImage, QPainter, QPen
from PySide6.QtWidgets import (
    QAbstractItemView, QApplication, QComboBox, QFrame, QHBoxLayout, QHeaderView, QLabel,
    QMessageBox, QPushButton, QSizePolicy, QTableWidget, QTableWidgetItem,
    QVBoxLayout, QWidget,
)

from .exports import group_color, make_gene_figure, significance_stars
from .mapping import plate_row, well_sort_key
from .models import AnalysisResult, ProjectState, WellRecord


class DropLabel(QLabel):
    fileDropped = Signal(str)
    clicked = Signal()

    def __init__(self) -> None:
        super().__init__('将 .xlsx / .csv 文件拖到这里\n或点击此区域选择机器数据')
        self.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.setAcceptDrops(True)
        self.setCursor(Qt.CursorShape.PointingHandCursor)
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self.setToolTip('点击整块区域即可选择 Excel 或 CSV 文件')
        self.setMinimumHeight(150)
        self.setObjectName('dropZone')

    def mousePressEvent(self, event) -> None:
        if event.button() == Qt.MouseButton.LeftButton:
            self.clicked.emit()
            event.accept()
            return
        super().mousePressEvent(event)

    def keyPressEvent(self, event) -> None:
        if event.key() in {Qt.Key.Key_Return, Qt.Key.Key_Enter, Qt.Key.Key_Space}:
            self.clicked.emit()
            event.accept()
            return
        super().keyPressEvent(event)

    def dragEnterEvent(self, event) -> None:
        if event.mimeData().hasUrls():
            event.acceptProposedAction()

    def dropEvent(self, event) -> None:
        path = event.mimeData().urls()[0].toLocalFile()
        if Path(path).suffix.lower() in {'.xlsx', '.xlsm', '.csv'}:
            self.fileDropped.emit(path)


class PlateTable(QTableWidget):
    selectionDescriptionChanged = Signal(str)

    def __init__(self) -> None:
        super().__init__()
        self.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        self.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectItems)
        self.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.setAlternatingRowColors(False)
        self.setStyleSheet(
            'QTableWidget::item:selected { background: rgba(37, 99, 235, 95); '
            'color: #0F172A; border: 1px solid #2563EB; }')
        self.itemSelectionChanged.connect(self._selection_changed)
        self.horizontalHeader().sectionClicked.connect(self._select_column)
        self.verticalHeader().sectionClicked.connect(self._select_row)
        self._row_names: list[str] = []
        self._columns = 0
        self._zoom = 85

    def set_plate(
        self, wells: tuple[WellRecord, ...], assignments, excluded=None,
        group_order=None,
    ) -> None:
        excluded = excluded or {}
        by_well = {record.well: record for record in wells}
        self._row_names = sorted({plate_row(record.well) for record in wells},
                                 key=lambda value: well_sort_key(value + '1')[0])
        self._columns = max((well_sort_key(record.well)[1] for record in wells), default=0)
        self.clear()
        self.setRowCount(len(self._row_names))
        self.setColumnCount(self._columns)
        self.setVerticalHeaderLabels(self._row_names)
        self.setHorizontalHeaderLabels([str(value) for value in range(1, self._columns + 1)])
        for row_index, row_name in enumerate(self._row_names):
            for column in range(1, self._columns + 1):
                well = f'{row_name}{column}'
                record = by_well.get(well)
                if not record:
                    continue
                assignment = assignments.get(well)
                cp = 'Undetermined' if record.ct is None else f'{record.ct:.2f}'
                lines = [cp, assignment.gene] if assignment and assignment.gene else [well, cp]
                item = QTableWidgetItem('\n'.join(lines))
                item.setData(Qt.ItemDataRole.UserRole, well)
                if assignment:
                    item.setToolTip(
                        f'{well}｜Cp: {cp}｜样本: {assignment.sample or "未分组"}｜'
                        f'组别: {assignment.group or "未分组"}｜基因: {assignment.gene or "未分配"}')
                item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
                if well in excluded:
                    font = item.font(); font.setStrikeOut(True); item.setFont(font)
                    item.setForeground(QColor('#94A3B8')); item.setBackground(QColor('#F1F5F9'))
                elif record.ct is None:
                    item.setForeground(QColor('#64748B')); item.setBackground(QColor('#F8FAFC'))
                elif assignment and assignment.role != 'sample':
                    item.setBackground(QColor('#F1F5F9'))
                elif assignment and assignment.group:
                    color = QColor(group_color(assignment.group, group_order))
                    color.setAlpha(38); item.setBackground(color)
                self.setItem(row_index, column - 1, item)
        self.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Fixed)
        self.verticalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Fixed)
        self.set_zoom(self._zoom)
        self.setMinimumHeight(470)
        self._selection_changed()

    def set_zoom(self, percent: int) -> None:
        self._zoom = max(55, min(150, int(percent)))
        width = max(42, round(72 * self._zoom / 100))
        height = max(34, round(54 * self._zoom / 100))
        self.horizontalHeader().setDefaultSectionSize(width)
        self.verticalHeader().setDefaultSectionSize(height)

    def fit_to_view(self) -> int:
        if not self.rowCount() or not self.columnCount():
            return self._zoom
        available_width = max(320, self.viewport().width() - 8)
        available_height = max(260, self.viewport().height() - 8)
        fitted_width = max(42, min(72, available_width // self.columnCount()))
        fitted_height = max(34, min(54, available_height // self.rowCount()))
        width_percent = fitted_width / 72 * 100
        height_percent = fitted_height / 54 * 100
        fitted = max(55, min(100, int(min(width_percent, height_percent))))
        self._zoom = fitted
        self.horizontalHeader().setDefaultSectionSize(fitted_width)
        self.verticalHeader().setDefaultSectionSize(fitted_height)
        return fitted

    def selected_wells(self) -> list[str]:
        wells = [item.data(Qt.ItemDataRole.UserRole) for item in self.selectedItems()
                 if item.data(Qt.ItemDataRole.UserRole)]
        return sorted(set(wells), key=well_sort_key)

    def select_wells(self, wells: list[str]) -> None:
        wanted = set(wells)
        self.clearSelection()
        for row in range(self.rowCount()):
            for column in range(self.columnCount()):
                item = self.item(row, column)
                if item and item.data(Qt.ItemDataRole.UserRole) in wanted:
                    item.setSelected(True)

    def _select_column(self, column: int) -> None:
        if not self.rowCount():
            return
        selection = QItemSelection(self.model().index(0, column),
                                   self.model().index(self.rowCount() - 1, column))
        flag = (QItemSelectionModel.SelectionFlag.Select
                if QApplication.keyboardModifiers() & Qt.KeyboardModifier.ControlModifier
                else QItemSelectionModel.SelectionFlag.ClearAndSelect)
        self.selectionModel().select(selection, flag)

    def _select_row(self, row: int) -> None:
        if not self.columnCount():
            return
        selection = QItemSelection(self.model().index(row, 0),
                                   self.model().index(row, self.columnCount() - 1))
        flag = (QItemSelectionModel.SelectionFlag.Select
                if QApplication.keyboardModifiers() & Qt.KeyboardModifier.ControlModifier
                else QItemSelectionModel.SelectionFlag.ClearAndSelect)
        self.selectionModel().select(selection, flag)

    def _selection_changed(self) -> None:
        wells = self.selected_wells()
        selected_rows = {plate_row(well) for well in wells}
        selected_columns = {well_sort_key(well)[1] for well in wells}
        for index, name in enumerate(self._row_names):
            header = self.verticalHeaderItem(index)
            if header:
                header.setBackground(QColor('#BFDBFE') if name in selected_rows else QColor('#FFFFFF'))
                header.setForeground(QColor('#1D4ED8') if name in selected_rows else QColor('#334155'))
        for index in range(self._columns):
            header = self.horizontalHeaderItem(index)
            if header:
                active = index + 1 in selected_columns
                header.setBackground(QColor('#BFDBFE') if active else QColor('#FFFFFF'))
                header.setForeground(QColor('#1D4ED8') if active else QColor('#334155'))
        if not wells:
            description = '未选择孔位'
        else:
            rows = ','.join(sorted(selected_rows, key=lambda r: well_sort_key(r+'1')[0]))
            cols = sorted(selected_columns)
            column_text = str(cols[0]) if len(cols) == 1 else f'{cols[0]}–{cols[-1]}'
            description = f'已选择 {len(wells)} 孔｜行 {rows}｜列 {column_text}'
        self.selectionDescriptionChanged.emit(description)
        self.viewport().update()

    def paintEvent(self, event) -> None:
        super().paintEvent(event)
        painter = QPainter(self.viewport())
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, False)
        painter.setPen(QPen(QColor('#1D4ED8'), 4))
        for selection_range in self.selectedRanges():
            top_left = self.item(selection_range.topRow(), selection_range.leftColumn())
            bottom_right = self.item(selection_range.bottomRow(), selection_range.rightColumn())
            if not top_left or not bottom_right:
                continue
            rect = self.visualItemRect(top_left).united(self.visualItemRect(bottom_right))
            painter.drawRect(rect.adjusted(1, 1, -2, -2))


def display_number(value, digits: int = 3) -> str:
    try:
        return f'{float(value):.{digits}f}' if np.isfinite(float(value)) else '不可计算'
    except (TypeError, ValueError):
        return '不可计算'


class GeneResultPage(QWidget):
    filterChanged = Signal(str, list)
    errorBarChanged = Signal(str)
    exportRequested = Signal(str)

    def __init__(self, gene: str) -> None:
        super().__init__()
        self.gene = gene
        self._updating = False
        self._state: ProjectState | None = None
        self._result: AnalysisResult | None = None
        self._all_samples: list[str] = []
        self._selected_samples: set[str] = set()
        root = QVBoxLayout(self)
        root.setContentsMargins(8, 8, 8, 8)
        root.setSpacing(7)
        top = QHBoxLayout()
        self.process = QLabel('原始 Cp  →  排除复核  →  复孔均值  →  ΔCt  →  校准  →  ΔΔCt  →  Fold change')
        self.process.setObjectName('processBar')
        self.process.setWordWrap(True)
        self.badge = QLabel('正式分析｜全部样本')
        self.badge.setObjectName('officialBadge')
        self.badge.setMaximumWidth(340)
        top.addWidget(self.process, 1); top.addWidget(self.badge)
        root.addLayout(top)

        self.table = QTableWidget()
        self.table.setColumnCount(8)
        self.table.setHorizontalHeaderLabels([
            '组别', '样本', f'{gene} Ct', '内参 Ct', 'ΔCt',
            'Ctrl 组平均 ΔCt', 'ΔΔCt', 'Fold change',
        ])
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setAlternatingRowColors(True)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.table.verticalHeader().setDefaultSectionSize(30)
        self.table.cellClicked.connect(self._show_source)
        self.table.setMinimumHeight(320)
        table_actions = QHBoxLayout()
        self.selection_hint = QLabel('当前使用 0 个样本｜已排除 0 个')
        self.selection_hint.setObjectName('hint')
        exclude_rows = QPushButton('排除选中')
        exclude_rows.setObjectName('warningButton')
        exclude_rows.clicked.connect(self._exclude_selected_rows)
        restore_rows = QPushButton('恢复选中')
        restore_rows.clicked.connect(self._restore_selected_rows)
        restore_all = QPushButton('全部恢复')
        restore_all.setObjectName('ghostButton')
        restore_all.clicked.connect(self._restore_all_samples)
        table_actions.addWidget(self.selection_hint, 1)
        table_actions.addWidget(exclude_rows)
        table_actions.addWidget(restore_rows)
        table_actions.addWidget(restore_all)
        root.addLayout(table_actions)
        root.addWidget(self.table, 10)

        chart_frame = QFrame(); chart_frame.setObjectName('card')
        chart_layout = QVBoxLayout(chart_frame)
        controls = QHBoxLayout()
        self.chart_toggle = QPushButton('显示图表')
        self.chart_toggle.setCheckable(True)
        self.chart_toggle.setObjectName('secondaryButton')
        self.chart_toggle.toggled.connect(self._toggle_chart)
        controls.addWidget(self.chart_toggle)
        controls.addWidget(QLabel('误差线'))
        self.error_combo = QComboBox(); self.error_combo.addItems(['SEM', 'SD', '95% CI', '无'])
        self.error_combo.currentTextChanged.connect(self.errorBarChanged.emit)
        self.significance_hint = QLabel('星号：原始 P')
        self.significance_hint.setObjectName('hint')
        copy_button = QPushButton('复制图片')
        copy_button.setObjectName('secondaryButton')
        copy_button.clicked.connect(lambda: self._copy_figure(copy_button))
        export_button = QPushButton('导出本基因 SVG')
        export_button.clicked.connect(lambda: self.exportRequested.emit(self.gene))
        controls.addWidget(self.error_combo); controls.addWidget(self.significance_hint)
        controls.addStretch()
        controls.addWidget(copy_button); controls.addWidget(export_button)
        chart_layout.addLayout(controls)
        self.canvas = FigureCanvasQTAgg()
        self.canvas.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Preferred)
        self.canvas.setMinimumSize(480, 460)
        self.canvas.setMaximumSize(660, 600)
        self.canvas.setVisible(False)
        chart_layout.addWidget(self.canvas, 0, Qt.AlignmentFlag.AlignHCenter)
        root.addWidget(chart_frame)

        self.stats_button = QPushButton('展开统计结果')
        self.stats_button.setCheckable(True)
        self.stats_button.toggled.connect(self._toggle_stats)
        self.stats_table = QTableWidget(); self.stats_table.setVisible(False)
        self.stats_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        root.addWidget(self.stats_button); root.addWidget(self.stats_table)

    def update_data(self, state: ProjectState, result: AnalysisResult) -> None:
        self._updating = True
        self._state, self._result = state, result
        rows = result.gene_tables.get(self.gene, [])
        refs = ', '.join(state.reference_genes()) or '未设置'
        calibrator = state.config.calibrator_group
        self.process.setText(
            f'原始 Cp  →  排除复核  →  复孔均值  →  ΔCt  →  {calibrator} 校准  →  ΔΔCt  →  Fold change')
        self.table.setHorizontalHeaderItem(3, QTableWidgetItem(f'内参 Ct ({refs})'))
        self.table.setHorizontalHeaderItem(5, QTableWidgetItem(f'{calibrator} 组平均 ΔCt'))
        self.table.clearSpans()
        self.table.setRowCount(len(rows))
        self._all_samples = [row['sample'] for row in rows]
        self._selected_samples = {
            row['sample'] for row in rows if row.get('included_for_analysis')
        }
        summary_map = {(row['sample'], row['gene']): row for row in result.replicate_summary}
        for row_index, row in enumerate(rows):
            sample_name = row['sample']
            values = [row['group'], row['sample'], display_number(row.get('target_ct')),
                      display_number(row.get('reference_ct')), display_number(row.get('delta_ct')),
                      display_number(row.get('calibrator_mean_delta_ct')),
                      display_number(row.get('delta_delta_ct')), display_number(row.get('fold_change'), 4)]
            for column, value in enumerate(values):
                item = QTableWidgetItem(value)
                item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
                item.setData(Qt.ItemDataRole.UserRole, row['sample'])
                if column == 0:
                    item.setBackground(QColor(group_color(row['group'], state.config.group_order)))
                    item.setForeground(QColor('#FFFFFF'))
                    font = item.font(); font.setBold(True); item.setFont(font)
                if column == 2:
                    item.setForeground(QColor('#1D4ED8'))
                    item.setToolTip(self._summary_text(summary_map.get((row['sample'], self.gene))))
                elif column == 3:
                    item.setForeground(QColor('#1D4ED8'))
                    text = [f'{ref}\n{self._summary_text(summary_map.get((row["sample"], ref)))}'
                            for ref in state.reference_genes()]
                    item.setToolTip('\n\n'.join(text))
                if not row.get('included_for_analysis'):
                    font = item.font(); font.setStrikeOut(True); item.setFont(font)
                    item.setForeground(QColor('#94A3B8'))
                    item.setBackground(QColor('#F1F5F9'))
                    if column == 1:
                        item.setText(f'{sample_name}（已排除）')
                self.table.setItem(row_index, column, item)
        group_rows: dict[str, list[int]] = defaultdict(list)
        for row_index, row in enumerate(rows):
            group_rows[row['group']].append(row_index)
        for indices in group_rows.values():
            if len(indices) > 1 and indices == list(range(indices[0], indices[-1] + 1)):
                self.table.setSpan(indices[0], 0, len(indices), 1)
        selected = [row['sample'] for row in rows if row.get('included_for_analysis')]
        excluded_count = max(0, len(rows) - len(selected))
        self.selection_hint.setText(
            f'当前使用 {len(selected)} 个样本｜已排除 {excluded_count} 个｜排除或恢复后立即重新计算')
        if self.gene in result.exploratory_genes:
            self.badge.setText(f'探索性筛选｜n={len(selected)}｜{", ".join(selected)}')
            self.badge.setObjectName('exploratoryBadge')
        else:
            self.badge.setText(f'正式分析｜全部样本｜n={len(selected)}')
            self.badge.setObjectName('officialBadge')
        if self.gene in result.exploratory_genes:
            self.badge.setText(f'已按排除样本重算｜n={len(selected)}')
            self.badge.setObjectName('exploratoryBadge')
        else:
            self.badge.setText(f'正式分析｜n={len(selected)}')
            self.badge.setObjectName('officialBadge')
        self.badge.setToolTip('图形、统计和导出均只使用计算表中未排除的样本')
        active_groups = {
            row.get('group') for row in rows if row.get('included_for_analysis')
        }
        calibrator = state.config.calibrator_group
        comparison_reference = (
            'Model' if 'Model' in active_groups
            else calibrator if calibrator in active_groups
            else next((group for group in state.config.group_order if group in active_groups), None)
        )
        if comparison_reference == 'Model':
            comparison_text = f'Model vs {calibrator}；其余组 vs Model'
        elif comparison_reference:
            comparison_text = f'各组 vs {comparison_reference}'
        else:
            comparison_text = '暂无可用对照组'
        self.significance_hint.setText(f'星号：原始 P｜{comparison_text}')
        self.badge.style().unpolish(self.badge); self.badge.style().polish(self.badge)
        self.error_combo.setCurrentText(state.config.error_bar)
        self._draw_figure(); self._populate_statistics(result)
        self._updating = False

    def _selected_table_samples(self) -> set[str]:
        samples: set[str] = set()
        for index in self.table.selectionModel().selectedRows(1):
            item = self.table.item(index.row(), 1)
            if item:
                sample = item.data(Qt.ItemDataRole.UserRole)
                if sample:
                    samples.add(str(sample))
        return samples

    def _exclude_selected_rows(self) -> None:
        samples = self._selected_table_samples()
        if not samples:
            self.selection_hint.setText('请先在计算表中选中一个或多个样本行。')
            return
        self.filterChanged.emit(self.gene, sorted(self._selected_samples - samples))

    def _restore_selected_rows(self) -> None:
        samples = self._selected_table_samples()
        if not samples:
            self.selection_hint.setText('请先在计算表中选中要恢复的样本行。')
            return
        self.filterChanged.emit(self.gene, sorted(self._selected_samples | samples))

    def _restore_all_samples(self) -> None:
        self.filterChanged.emit(self.gene, list(self._all_samples))

    def _copy_figure(self, button: QPushButton) -> None:
        if not self.canvas.figure:
            return
        buffer = BytesIO()
        self.canvas.figure.savefig(
            buffer, format='png', dpi=220, bbox_inches='tight', facecolor='white')
        image = QImage.fromData(buffer.getvalue(), 'PNG')
        if image.isNull():
            QMessageBox.warning(self, '复制失败', '图片未能写入剪贴板，请改用 SVG 导出。')
            return
        QApplication.clipboard().setImage(image)
        original = button.text()
        button.setText('已复制到剪贴板')
        QTimer.singleShot(1600, lambda: button.setText(original))

    def _toggle_chart(self, checked: bool) -> None:
        self.canvas.setVisible(checked)
        self.chart_toggle.setText('隐藏图表' if checked else '显示图表')
        if checked:
            self.canvas.draw_idle()

    def _draw_figure(self) -> None:
        if not self._result:
            return
        old_figure = self.canvas.figure
        figure = make_gene_figure(
            self._result, self.gene, self.error_combo.currentText(),
            self._state.config.group_order if self._state else None,
            self._state.config.calibrator_group if self._state else 'Ctrl')
        self.canvas.figure = figure
        figure.set_canvas(self.canvas)
        self.canvas.draw_idle()
        try:
            old_figure.clear()
        except Exception:
            pass

    def _populate_statistics(self, result: AnalysisResult) -> None:
        rows = [row for row in result.statistics if row.get('gene') == self.gene]
        keys = [
            'test', 'comparison', 'statistic', 'df', 'p_value', 'p_adjust_bh',
            'plot_stars', 'significant',
        ]
        self.stats_table.clear()
        self.stats_table.setColumnCount(len(keys)); self.stats_table.setRowCount(len(rows))
        self.stats_table.setHorizontalHeaderLabels(
            ['检验', '比较', '统计量', '自由度', 'P 值', 'BH 校正 P',
             '图中星号', 'BH 显著'])
        for row_index, row in enumerate(rows):
            for column_index, key in enumerate(keys):
                value = row.get(key, '')
                if key in {'statistic', 'p_value', 'p_adjust_bh'}:
                    value = display_number(value, 4)
                elif key == 'plot_stars':
                    value = significance_stars(row.get('p_value')) or '—'
                elif key == 'significant':
                    value = '是' if value else '否'
                self.stats_table.setItem(row_index, column_index, QTableWidgetItem(str(value)))
        if not rows:
            self.stats_button.setText('统计结果：样本量不足，仅显示描述统计')
            self.stats_button.setEnabled(False)
        else:
            self.stats_button.setEnabled(True)
            self.stats_button.setText('收起统计结果' if self.stats_button.isChecked() else '展开统计结果')
        self.stats_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)

    def _toggle_stats(self, checked: bool) -> None:
        self.stats_table.setVisible(checked)
        if self.stats_button.isEnabled():
            self.stats_button.setText('收起统计结果' if checked else '展开统计结果')

    @staticmethod
    def _summary_text(summary: dict | None) -> str:
        if not summary:
            return '无复孔来源'
        kept = ', '.join(f'{well}={value:.3f}' for well, value in zip(
            summary.get('retained_wells', []), summary.get('retained_ct', []))) or '无'
        removed = ', '.join(f'{well}={value:.3f}' for well, value in zip(
            summary.get('excluded_wells', []), summary.get('excluded_ct', []))) or '无'
        return f'参与均值：{kept}\n已排除：{removed}'

    def _show_source(self, row: int, column: int) -> None:
        if column not in {2, 3}:
            return
        item = self.table.item(row, column)
        if item:
            QMessageBox.information(self, 'Ct 均值来源', item.toolTip() or '无复孔来源')
