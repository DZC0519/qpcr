from __future__ import annotations

import json
import sys
from copy import deepcopy
from pathlib import Path

from PySide6.QtCore import QSettings, QTimer, Qt
from PySide6.QtGui import QColor, QFont, QFontDatabase
from PySide6.QtWidgets import (
    QAbstractItemView, QApplication, QCheckBox, QComboBox, QDoubleSpinBox,
    QFileDialog, QFormLayout, QFrame, QHBoxLayout, QHeaderView, QLabel,
    QInputDialog, QListWidget, QMainWindow, QMessageBox, QPushButton, QScrollArea,
    QSpinBox, QStackedWidget, QTabWidget, QTableWidget, QTableWidgetItem,
    QVBoxLayout, QWidget,
)

from .analysis import analyze
from .exports import export_all_svg, export_clean_excel, export_figure
from .importer import find_candidates, load_plate
from .mapping import (
    apply_gene_to_wells, apply_group_to_axis, apply_group_to_selected_wells,
    detect_sample_axis, load_template, normalize_mapping_scope,
    normalize_sample_axis, plate_row, renumber_samples, sample_axis_key,
    save_template, suggest_row_layout, well_sort_key,
)
from .models import (
    GROUP_ORDER, AnalysisConfig, GeneDefinition, ProjectState, normalize_group_order,
)
from .project import load_project, save_project
from .qc import build_replicate_qc
from .widgets import DropLabel, GeneResultPage, PlateTable, display_number


STEPS = ['1  导入与实验设置', '2  板图分配', '3  复孔 QC', '4  分析与作图', '5  导出与工程']


def install_chinese_font() -> None:
    app = QApplication.instance()
    if not app:
        return
    for candidate in ('C:/Windows/Fonts/msyh.ttc', 'C:/Windows/Fonts/simhei.ttf'):
        if not Path(candidate).exists():
            continue
        font_id = QFontDatabase.addApplicationFont(candidate)
        families = QFontDatabase.applicationFontFamilies(font_id)
        if families:
            app.setFont(QFont(families[0], 10)); return


class MainWindow(QMainWindow):
    def __init__(self) -> None:
        super().__init__()
        install_chinese_font()
        self.state = ProjectState(config=AnalysisConfig())
        self.result = None
        self.settings = QSettings('LocalLabTools', 'qPCRAnalyzerV2')
        self.undo_stack: list[dict] = []
        self.redo_stack: list[dict] = []
        self.result_pages: dict[str, GeneResultPage] = {}
        self.setWindowTitle('qPCR 分析助手 2.1')
        self.resize(1280, 780)
        self.setMinimumSize(1024, 650)
        self._build_ui(); self._apply_style()

    def _page(self, title: str, subtitle: str):
        page = QScrollArea(); page.setWidgetResizable(True); page.setFrameShape(QFrame.Shape.NoFrame)
        content = QWidget(); layout = QVBoxLayout(content); page.setWidget(content)
        layout.setContentsMargins(14, 10, 14, 14)
        layout.setSpacing(8)
        heading = QLabel(title); heading.setObjectName('heading')
        note = QLabel(subtitle); note.setObjectName('subtitle'); note.setWordWrap(True)
        layout.addWidget(heading); layout.addWidget(note)
        return page, layout

    def _build_ui(self) -> None:
        root = QWidget(); outer = QHBoxLayout(root); outer.setContentsMargins(0, 0, 0, 0)
        self.nav = QListWidget(); self.nav.addItems(STEPS); self.nav.setFixedWidth(195)
        self.nav.setObjectName('nav'); self.nav.setCurrentRow(0)
        self.stack = QStackedWidget()
        outer.addWidget(self.nav); outer.addWidget(self.stack, 1)
        self.setCentralWidget(root)
        self.nav.currentRowChanged.connect(self.stack.setCurrentIndex)
        self._build_import_page(); self._build_mapping_page(); self._build_qc_page()
        self._build_analysis_page(); self._build_export_page()

    def _build_import_page(self) -> None:
        page, layout = self._page(
            '导入与实验设置',
            '先导入机器数据并录入本次实验基因；板图支持 Ctrl 跨行/跨列多选，并把组别、样本号和基因精确应用到所选孔。')
        self.drop = DropLabel(); self.drop.fileDropped.connect(self.import_file)
        self.drop.clicked.connect(self.choose_import)
        layout.addWidget(self.drop)
        buttons = QHBoxLayout()
        choose = QPushButton('选择机器数据'); choose.clicked.connect(self.choose_import)
        open_project = QPushButton('打开分析工程'); open_project.clicked.connect(self.open_project)
        buttons.addWidget(choose); buttons.addWidget(open_project); buttons.addStretch()
        layout.addLayout(buttons)
        self.import_info = QLabel('尚未导入数据'); self.import_info.setObjectName('statusCard')
        layout.addWidget(self.import_info)

        group_card = QFrame(); group_card.setObjectName('card')
        group_layout = QVBoxLayout(group_card)
        group_title = QHBoxLayout(); group_title.addWidget(QLabel('本次实验组别'))
        group_title.addStretch()
        add_group = QPushButton('添加组别')
        add_group.clicked.connect(lambda: self.add_group_row(start_edit=True))
        remove_group = QPushButton('删除选中'); remove_group.clicked.connect(self.remove_group_rows)
        restore_groups = QPushButton('恢复默认组别')
        restore_groups.clicked.connect(self.restore_default_groups)
        apply_groups = QPushButton('应用组别设置')
        apply_groups.setObjectName('secondaryButton')
        apply_groups.clicked.connect(self.commit_group_definitions)
        for button in (add_group, remove_group, restore_groups, apply_groups):
            group_title.addWidget(button)
        group_layout.addLayout(group_title)
        self.group_table = QTableWidget(0, 1)
        self.group_table.setHorizontalHeaderLabels(['组别名称'])
        self.group_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)
        self.group_table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.group_table.setEditTriggers(
            QAbstractItemView.EditTrigger.SelectedClicked
            | QAbstractItemView.EditTrigger.EditKeyPressed
            | QAbstractItemView.EditTrigger.AnyKeyPressed
        )
        self.group_table.cellClicked.connect(self._edit_group_cell)
        group_layout.addWidget(self.group_table)
        layout.addWidget(group_card)
        self._load_group_table_from_state()

        gene_card = QFrame(); gene_card.setObjectName('card')
        gene_layout = QVBoxLayout(gene_card)
        title = QHBoxLayout(); title.addWidget(QLabel('本次实验基因'))
        title.addStretch()
        add_gene = QPushButton('添加基因')
        add_gene.clicked.connect(lambda: self.add_gene_row(start_edit=True))
        remove_gene = QPushButton('删除选中'); remove_gene.clicked.connect(self.remove_gene_rows)
        save_common = QPushButton('保存为常用列表'); save_common.clicked.connect(self.save_common_genes)
        load_common = QPushButton('载入常用列表'); load_common.clicked.connect(self.load_common_genes)
        for button in (add_gene, remove_gene, save_common, load_common): title.addWidget(button)
        gene_layout.addLayout(title)
        self.gene_table = QTableWidget(0, 2)
        self.gene_table.setHorizontalHeaderLabels(['基因名称', '设为内参'])
        self.gene_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)
        self.gene_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeMode.ResizeToContents)
        self.gene_table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.gene_table.setEditTriggers(
            QAbstractItemView.EditTrigger.SelectedClicked
            | QAbstractItemView.EditTrigger.EditKeyPressed
            | QAbstractItemView.EditTrigger.AnyKeyPressed
        )
        self.gene_table.cellClicked.connect(self._edit_gene_cell)
        gene_layout.addWidget(self.gene_table)
        confirm = QPushButton('确认基因设置并进入板图分配')
        confirm.setObjectName('primaryButton'); confirm.clicked.connect(self.commit_genes)
        gene_layout.addWidget(confirm)
        layout.addWidget(gene_card, 1)
        self.stack.addWidget(page)

    def add_group_row(self, name: str = '', start_edit: bool = False) -> None:
        row = self.group_table.rowCount()
        self.group_table.insertRow(row)
        item = QTableWidgetItem(name)
        item.setData(Qt.ItemDataRole.UserRole, name)
        self.group_table.setItem(row, 0, item)
        if start_edit:
            self.group_table.setCurrentCell(row, 0)
            QTimer.singleShot(0, lambda current=item: self.group_table.editItem(current))

    def _edit_group_cell(self, row: int, column: int) -> None:
        if column != 0:
            return
        item = self.group_table.item(row, column)
        if item:
            self.group_table.editItem(item)

    def _group_rows(self) -> list[str]:
        names: list[str] = []
        for row in range(self.group_table.rowCount()):
            item = self.group_table.item(row, 0)
            name = item.text().strip() if item else ''
            if name:
                names.append(name)
        return names

    def remove_group_rows(self) -> None:
        for row in sorted({index.row() for index in self.group_table.selectedIndexes()}, reverse=True):
            self.group_table.removeRow(row)

    def restore_default_groups(self) -> None:
        self.group_table.setRowCount(0)
        for name in GROUP_ORDER:
            self.add_group_row(name)

    def _load_group_table_from_state(self) -> None:
        self.group_table.setRowCount(0)
        for name in normalize_group_order(self.state.config.group_order):
            self.add_group_row(name)

    def _refresh_group_combo(self) -> None:
        if not hasattr(self, 'group_combo'):
            return
        current = self.group_combo.currentText()
        groups = normalize_group_order(self.state.config.group_order)
        self.group_combo.blockSignals(True)
        self.group_combo.clear(); self.group_combo.addItems(groups)
        self.group_combo.setCurrentText(current if current in groups else (groups[0] if groups else ''))
        self.group_combo.blockSignals(False)

    def commit_group_definitions(self, silent: bool = False) -> bool:
        previous_order = normalize_group_order(self.state.config.group_order)
        names: list[str] = []
        original_names: list[str] = []
        for row in range(self.group_table.rowCount()):
            item = self.group_table.item(row, 0)
            if not item:
                continue
            name = item.text().strip()
            if not name:
                continue
            names.append(name)
            original_names.append(str(item.data(Qt.ItemDataRole.UserRole) or name).strip())
        if not names:
            if not silent:
                QMessageBox.warning(self, '需要组别', '请至少保留一个实验组别。')
            return False
        if len({name.casefold() for name in names}) != len(names):
            if not silent:
                QMessageBox.warning(self, '名称重复', '组别名称不能重复。')
            return False
        new_names = set(names)
        self._push_undo()
        renames = {
            old: new for old, new in zip(original_names, names)
            if old and old != new
        }
        for assignment in self.state.assignments.values():
            old_group = assignment.group
            if old_group in renames:
                new_group = renames[old_group]
            elif old_group and old_group not in new_names and old_group in previous_order:
                new_group = ''
            else:
                new_group = old_group
            if new_group != old_group:
                assignment.group = new_group
                if assignment.sample and old_group:
                    prefix = f'{old_group}-'
                    if assignment.sample.startswith(prefix):
                        assignment.sample = f'{new_group}-{assignment.sample[len(prefix):]}' if new_group else assignment.sample
        calibrator = self.state.config.calibrator_group
        if calibrator in renames:
            self.state.config.calibrator_group = renames[calibrator]
        elif calibrator not in new_names:
            self.state.config.calibrator_group = names[0]
        self.state.config.group_order = normalize_group_order(names)
        for row in range(self.group_table.rowCount()):
            item = self.group_table.item(row, 0)
            if item and item.text().strip():
                item.setData(Qt.ItemDataRole.UserRole, item.text().strip())
        self._refresh_group_combo()
        self._refresh_plate()
        if not silent:
            self.statusBar().showMessage(f'已应用 {len(names)} 个实验组别', 3500)
        return True

    def add_gene_row(
        self, name: str = '', is_reference: bool = False, start_edit: bool = False,
    ) -> None:
        row = self.gene_table.rowCount(); self.gene_table.insertRow(row)
        item = QTableWidgetItem(name)
        self.gene_table.setItem(row, 0, item)
        checkbox = QCheckBox(); checkbox.setChecked(is_reference)
        holder = QWidget(); box = QHBoxLayout(holder); box.setContentsMargins(0, 0, 0, 0)
        box.setAlignment(Qt.AlignmentFlag.AlignCenter); box.addWidget(checkbox)
        self.gene_table.setCellWidget(row, 1, holder)
        if start_edit:
            self.gene_table.setCurrentCell(row, 0)
            QTimer.singleShot(0, lambda current=item: self.gene_table.editItem(current))

    def _edit_gene_cell(self, row: int, column: int) -> None:
        if column != 0:
            return
        item = self.gene_table.item(row, column)
        if item:
            self.gene_table.editItem(item)

    def _gene_rows(self) -> list[GeneDefinition]:
        genes = []
        for row in range(self.gene_table.rowCount()):
            name = (self.gene_table.item(row, 0).text() if self.gene_table.item(row, 0) else '').strip()
            holder = self.gene_table.cellWidget(row, 1)
            checkbox = holder.findChild(QCheckBox) if holder else None
            if name:
                genes.append(GeneDefinition(name, bool(checkbox and checkbox.isChecked()), len(genes)))
        return genes

    def remove_gene_rows(self) -> None:
        for row in sorted({index.row() for index in self.gene_table.selectedIndexes()}, reverse=True):
            self.gene_table.removeRow(row)

    def save_common_genes(self) -> None:
        genes = self._gene_rows()
        self.settings.setValue('common_genes', json.dumps([
            {'name': gene.name, 'is_reference': gene.is_reference} for gene in genes], ensure_ascii=False))
        QMessageBox.information(self, '已保存', '常用基因列表已保存在本机。')

    def load_common_genes(self) -> None:
        try:
            payload = json.loads(self.settings.value('common_genes', '[]'))
        except (TypeError, json.JSONDecodeError):
            payload = []
        if not payload:
            QMessageBox.information(self, '没有列表', '本机尚未保存常用基因列表。'); return
        self.gene_table.setRowCount(0)
        for gene in payload:
            self.add_gene_row(gene.get('name', ''), bool(gene.get('is_reference')))

    def commit_genes(self) -> None:
        if not self.commit_group_definitions(silent=True):
            QMessageBox.warning(self, '需要组别', '请先在上方保留至少一个有效组别。')
            return
        genes = self._gene_rows()
        if not genes:
            QMessageBox.warning(self, '需要基因', '请至少录入一个目标基因和一个内参基因。'); return
        if len({gene.name.casefold() for gene in genes}) != len(genes):
            QMessageBox.warning(self, '名称重复', '基因名称不能重复。'); return
        if not any(gene.is_reference for gene in genes):
            QMessageBox.warning(self, '需要内参', '请至少勾选一个内参基因。'); return
        if all(gene.is_reference for gene in genes):
            QMessageBox.warning(self, '需要目标基因', '请至少保留一个非内参的目标基因。'); return
        old = [gene.name for gene in self.state.genes]
        rename = {old[index]: gene.name for index, gene in enumerate(genes) if index < len(old)}
        for assignment in self.state.assignments.values():
            if assignment.gene in rename:
                assignment.gene = rename[assignment.gene]
        self.state.genes = genes; self.state.sync_gene_flags()
        self._refresh_gene_combo(); self._refresh_plate(); self.nav.setCurrentRow(1)

    def choose_import(self) -> None:
        path, _ = QFileDialog.getOpenFileName(self, '选择机器数据', '', '数据文件 (*.xlsx *.xlsm *.csv)')
        if path:
            self.import_file(path)

    def import_file(self, path: str) -> None:
        try:
            candidates = find_candidates(path)
            if not candidates:
                raise ValueError('没有找到同时包含孔位和 Cp/Cq/Ct 的表头。')
            selected = candidates[0]
            if len(candidates) > 1:
                labels = [f'{item.sheet}（表头第 {item.header_row} 行）' for item in candidates]
                label, ok = QInputDialog.getItem(self, '选择工作表', '检测到多个候选表：', labels, 0, False)
                if not ok:
                    return
                selected = candidates[labels.index(label)]
            plate = load_plate(path, selected.sheet, selected.header_row)
        except Exception as exc:
            QMessageBox.critical(self, '导入失败', str(exc)); return
        self.state = ProjectState(
            source_path=plate.source_path, source_sheet=plate.sheet, header_row=plate.header_row,
            wells=plate.wells, assignments=suggest_row_layout(plate.wells, empty_rows=()),
            config=AnalysisConfig(plate_format=plate.plate_format, sample_axis='auto'))
        self.result = None; self.undo_stack.clear(); self.redo_stack.clear()
        self._load_group_table_from_state(); self._refresh_group_combo()
        detected = []
        for record in plate.wells:
            name = record.machine_name.strip()
            if name and name.lower() not in {'sample', 'unknown'} and name not in detected:
                detected.append(name)
        self.gene_table.setRowCount(0)
        if 1 < len(detected) <= 24:
            for index, name in enumerate(detected): self.add_gene_row(name, index == len(detected)-1)
        else:
            self.add_gene_row('目标基因1', False); self.add_gene_row('内参基因1', True)
        valid = sum(record.ct is not None for record in plate.wells)
        undetermined = len(plate.wells) - valid
        self.import_info.setText(
            f'{Path(path).name}｜{plate.sheet}｜表头第 {plate.header_row} 行｜'
            f'{plate.plate_format} 孔｜有效 Cp {valid}｜Undetermined {undetermined}')
        self._sync_sample_axis_ui(); self._refresh_plate(); self.nav.setCurrentRow(0)

    def _build_mapping_page(self) -> None:
        page, layout = self._page(
            '板图分配',
            '框选或按住 Ctrl 选择任意孔块后，可将组别、样本号和基因精确应用到所选孔；技术复孔方向自动识别。')
        toolbar_card = QFrame(); toolbar_card.setObjectName('card')
        toolbar_layout = QVBoxLayout(toolbar_card)

        orientation_row = QHBoxLayout()
        orientation_row.addWidget(QLabel('样本方向'))
        self.sample_axis_status = QLabel('等待自动识别')
        self.sample_axis_status.setObjectName('axisBadge')
        orientation_row.addWidget(self.sample_axis_status)
        self.sample_axis_hint = QLabel('首次选择横向或纵向连续孔段后，软件会自动判断技术复孔方向。')
        self.sample_axis_hint.setObjectName('hint')
        orientation_row.addWidget(self.sample_axis_hint, 1)
        reset_axis = QPushButton('重新识别方向')
        reset_axis.setObjectName('ghostButton')
        reset_axis.clicked.connect(self.reset_sample_axis_detection)
        orientation_row.addWidget(reset_axis)
        toolbar_layout.addLayout(orientation_row)

        assignment_row = QHBoxLayout()
        assignment_row.addWidget(QLabel('固定组别'))
        self.group_combo = QComboBox(); self.group_combo.addItems(GROUP_ORDER)
        assignment_row.addWidget(self.group_combo)
        assignment_row.addWidget(QLabel('样本号'))
        self.sample_number = QSpinBox(); self.sample_number.setRange(1, 999)
        self.sample_number.setValue(1); self.sample_number.setFixedWidth(72)
        self.sample_number.setToolTip('生成样本名称，例如 Ctrl-1、Model-2')
        assignment_row.addWidget(self.sample_number)
        self.apply_group_button = QPushButton('应用组别/样本到所选孔')
        self.apply_group_button.clicked.connect(self.apply_group)
        self.apply_group_button.setObjectName('secondaryButton')
        assignment_row.addWidget(self.apply_group_button)
        self.axis_group_button = QPushButton('整行/整列分组')
        self.axis_group_button.clicked.connect(self.apply_group_by_axis)
        self.axis_group_button.setObjectName('ghostButton')
        self.axis_group_button.setToolTip('兼容旧用法：把组别扩展到选区所在的完整板行或板列')
        assignment_row.addWidget(self.axis_group_button)
        assignment_row.addStretch()
        toolbar_layout.addLayout(assignment_row)

        gene_row = QHBoxLayout(); gene_row.addWidget(QLabel('基因'))
        self.gene_combo = QComboBox(); gene_row.addWidget(self.gene_combo)
        apply_gene = QPushButton('应用基因到所选孔'); apply_gene.clicked.connect(self.apply_gene)
        apply_gene.setObjectName('secondaryButton')
        gene_row.addWidget(apply_gene); gene_row.addStretch()
        toolbar_layout.addLayout(gene_row)

        action_row = QHBoxLayout(); action_row.addWidget(QLabel('板图缩放'))
        self.plate_zoom = QSpinBox(); self.plate_zoom.setRange(55, 150); self.plate_zoom.setValue(85)
        self.plate_zoom.setSuffix('%'); self.plate_zoom.valueChanged.connect(self.plate_table_zoom_changed)
        action_row.addWidget(self.plate_zoom)
        fit = QPushButton('适应窗口'); fit.setObjectName('ghostButton')
        fit.clicked.connect(self.fit_plate_to_window); action_row.addWidget(fit)
        undo = QPushButton('撤销'); undo.clicked.connect(self.undo)
        redo = QPushButton('重做'); redo.clicked.connect(self.redo)
        save = QPushButton('保存模板'); save.clicked.connect(self.save_mapping_template)
        load = QPushButton('载入模板'); load.clicked.connect(self.load_mapping_template)
        for button in (undo, redo, save, load): action_row.addWidget(button)
        action_row.addStretch()
        next_button = QPushButton('进入复孔 QC'); next_button.setObjectName('primaryButton')
        next_button.clicked.connect(lambda: (self.refresh_qc(), self.nav.setCurrentRow(2)))
        action_row.addWidget(next_button)
        toolbar_layout.addLayout(action_row)
        layout.addWidget(toolbar_card)
        self.selection_label = QLabel('未选择孔位'); self.selection_label.setObjectName('selectionStatus')
        layout.addWidget(self.selection_label)
        self.plate_table = PlateTable()
        self.plate_table.selectionDescriptionChanged.connect(self.selection_label.setText)
        layout.addWidget(self.plate_table, 1)
        self.stack.addWidget(page)
        self._sync_sample_axis_ui()

    def plate_table_zoom_changed(self, value: int) -> None:
        if hasattr(self, 'plate_table'):
            self.plate_table.set_zoom(value)

    def fit_plate_to_window(self) -> None:
        if not hasattr(self, 'plate_table'):
            return
        fitted = self.plate_table.fit_to_view()
        self.plate_zoom.blockSignals(True)
        self.plate_zoom.setValue(fitted)
        self.plate_zoom.blockSignals(False)
        self.statusBar().showMessage(f'板图已适应当前窗口：{fitted}%', 3500)

    def _snapshot(self) -> dict:
        return {
            'assignments': deepcopy(self.state.assignments),
            'excluded_wells': deepcopy(self.state.excluded_wells),
            'exclusions': deepcopy(self.state.exclusions),
            'sample_axis': self.state.config.sample_axis,
            'mapping_scope': self.state.config.mapping_scope,
            'group_order': list(self.state.config.group_order),
            'plot_selections': deepcopy(self.state.plot_selections),
        }

    def _push_undo(self) -> None:
        self.undo_stack.append(self._snapshot())
        if len(self.undo_stack) > 100:
            self.undo_stack.pop(0)
        self.redo_stack.clear()

    def _restore_snapshot(self, payload: dict) -> None:
        self.state.assignments = deepcopy(payload['assignments'])
        self.state.excluded_wells = deepcopy(payload['excluded_wells'])
        self.state.exclusions = deepcopy(payload['exclusions'])
        self.state.config.sample_axis = payload.get('sample_axis', 'row')
        self.state.config.mapping_scope = payload.get('mapping_scope', 'axis')
        self.state.config.group_order = normalize_group_order(
            payload.get('group_order', self.state.config.group_order))
        self.state.plot_selections = deepcopy(payload.get('plot_selections', {}))
        self.result = None
        self._load_group_table_from_state(); self._refresh_group_combo()
        self._sync_sample_axis_ui()
        self._refresh_plate(); self.refresh_qc()

    def undo(self) -> None:
        if not self.undo_stack:
            return
        self.redo_stack.append(self._snapshot())
        self._restore_snapshot(self.undo_stack.pop())

    def redo(self) -> None:
        if not self.redo_stack:
            return
        self.undo_stack.append(self._snapshot())
        self._restore_snapshot(self.redo_stack.pop())

    def _refresh_gene_combo(self) -> None:
        current = self.gene_combo.currentText() if hasattr(self, 'gene_combo') else ''
        if not hasattr(self, 'gene_combo'):
            return
        self.gene_combo.clear(); self.gene_combo.addItems([gene.name for gene in self.state.genes])
        if current:
            self.gene_combo.setCurrentText(current)

    def _refresh_plate(self) -> None:
        if not hasattr(self, 'plate_table'):
            return
        self.plate_table.set_plate(
            self.state.wells, self.state.assignments, self.state.excluded_wells,
            self.state.config.group_order)

    def _sync_sample_axis_ui(self) -> None:
        if not hasattr(self, 'sample_axis_status'):
            return
        axis = normalize_sample_axis(self.state.config.sample_axis, 'auto')
        scope = normalize_mapping_scope(self.state.config.mapping_scope, 'axis')
        if scope == 'well' and axis == 'row':
            self.sample_axis_status.setText('精确选孔｜复孔横向')
            self.sample_axis_hint.setText('组别和样本只写入所选孔；每个样本的技术复孔从左到右编号。')
        elif scope == 'well' and axis == 'column':
            self.sample_axis_status.setText('精确选孔｜复孔纵向')
            self.sample_axis_hint.setText('组别和样本只写入所选孔；每个样本的技术复孔从上到下编号。')
        elif scope == 'well':
            self.sample_axis_status.setText('精确选孔｜等待识别复孔方向')
            self.sample_axis_hint.setText('请选择一段横向或纵向连续孔位，软件会自动识别复孔方向。')
        elif axis == 'row':
            self.sample_axis_status.setText('已识别：样本按行')
            self.sample_axis_hint.setText('旧式整行分组可用；基因复孔按从左到右编号。')
        elif axis == 'column':
            self.sample_axis_status.setText('已识别：样本按列')
            self.sample_axis_hint.setText('旧式整列分组可用；基因复孔按从上到下编号。')
        else:
            self.sample_axis_status.setText('等待自动识别')
            self.sample_axis_hint.setText('首次选择横向或纵向连续孔段后，软件会自动判断复孔方向。')
        self.sample_axis_status.style().unpolish(self.sample_axis_status)
        self.sample_axis_status.style().polish(self.sample_axis_status)

    def _axis_for_selection(self, wells: list[str], purpose: str) -> str | None:
        current = normalize_sample_axis(self.state.config.sample_axis, 'auto')
        if current in {'row', 'column'}:
            return current
        candidate = detect_sample_axis(
            wells, purpose, self.plate_table.rowCount(), self.plate_table.columnCount())
        if candidate:
            return candidate
        action = ('选择一条纵向样本列表或横向样本列表'
                  if purpose == 'group' else '选择一条横向或纵向连续复孔')
        QMessageBox.information(
            self, '暂时无法判断方向',
            f'当前选区同时跨越多行和多列，无法可靠判断样本方向。请先{action}完成第一次应用；'
            '方向识别后即可正常框选矩形区域。')
        return None

    def _commit_sample_axis(self, axis: str, mapping_scope: str | None = None) -> None:
        scope = normalize_mapping_scope(mapping_scope, self.state.config.mapping_scope)
        if self.state.config.sample_axis == axis and self.state.config.mapping_scope == scope:
            return
        self.state.config.sample_axis = axis
        self.state.config.mapping_scope = scope
        if scope == 'axis':
            renumber_samples(self.state.assignments, axis)
        self._sync_sample_axis_ui()

    def reset_sample_axis_detection(self) -> None:
        mapped = any(
            assignment.role == 'sample' and (assignment.group or assignment.gene)
            for assignment in self.state.assignments.values())
        if mapped:
            answer = QMessageBox.question(
                self, '重新识别样本方向',
                '重新识别会清除当前板图的组别、基因和复孔编号，但不会修改原始机器数据。继续吗？',
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                QMessageBox.StandardButton.No)
            if answer != QMessageBox.StandardButton.Yes:
                return
        self._push_undo()
        for assignment in self.state.assignments.values():
            if assignment.role != 'sample':
                continue
            assignment.group = ''
            assignment.gene = ''
            assignment.is_reference = False
            assignment.technical_replicate = None
            assignment.sample = f'样本{plate_row(assignment.well)}'
        self.state.config.sample_axis = 'auto'
        self.state.config.mapping_scope = 'axis'
        self.state.plot_selections.clear()
        self.result = None
        self._sync_sample_axis_ui(); self._refresh_plate()
        self.statusBar().showMessage('样本方向已重置，将根据下一次组别或基因选区自动识别。', 5000)

    def apply_group(self) -> None:
        wells = self.plate_table.selected_wells()
        if not wells:
            QMessageBox.information(self, '未选择', '请先框选孔位，或按住 Ctrl 选择多个孔块。'); return
        axis = self._axis_for_selection(wells, 'group')
        if not axis:
            return
        self._push_undo()
        self._commit_sample_axis(axis, 'well')
        targets = apply_group_to_selected_wells(
            self.state.assignments, wells, self.group_combo.currentText(),
            self.sample_number.value())
        self._refresh_plate()
        self.plate_table.select_wells(wells)
        sample_id = f'{self.group_combo.currentText()}-{self.sample_number.value()}'
        self.statusBar().showMessage(
            f'已将 {sample_id} 精确应用到 {len(targets)} 个所选孔，未选孔保持不变', 5000)

    def apply_group_by_axis(self) -> None:
        wells = self.plate_table.selected_wells()
        if not wells:
            QMessageBox.information(self, '未选择', '请先选择要整行或整列分组的孔位。'); return
        axis = self._axis_for_selection(wells, 'group')
        if not axis:
            return
        self._push_undo()
        self._commit_sample_axis(axis, 'axis')
        targets = apply_group_to_axis(
            self.state.assignments, wells, self.group_combo.currentText(), axis)
        self._refresh_plate()
        self.plate_table.select_wells(wells)
        ordered = sorted(targets, key=(lambda value: int(value)) if axis == 'column' else None)
        unit = '列' if axis == 'column' else '行'
        self.statusBar().showMessage(
            f'已按旧式布局将 {self.group_combo.currentText()} 应用到{unit} {", ".join(ordered)}',
            5000)

    def apply_gene(self) -> None:
        wells = self.plate_table.selected_wells(); gene = self.gene_combo.currentText()
        if not wells or not gene:
            QMessageBox.information(self, '未选择', '请先选择基因并框选连续孔位。'); return
        axis = self._axis_for_selection(wells, 'gene')
        if not axis:
            return
        positions_by_axis: dict[str, list[int]] = {}
        for well in wells:
            row_index, column = well_sort_key(well)
            key = plate_row(well) if axis == 'row' else str(column)
            position = column if axis == 'row' else row_index
            positions_by_axis.setdefault(key, []).append(position)
        if any(sorted(positions) != list(range(min(positions), max(positions) + 1))
               for positions in positions_by_axis.values()):
            direction = '横向' if axis == 'row' else '纵向'
            QMessageBox.warning(
                self, '孔位不连续', f'同一样本内的技术复孔必须为{direction}连续孔，请重新框选。')
            return
        self._push_undo()
        scope = normalize_mapping_scope(self.state.config.mapping_scope, 'axis')
        self._commit_sample_axis(axis, scope)
        definition = next((item for item in self.state.genes if item.name == gene), None)
        apply_gene_to_wells(
            self.state.assignments, wells, gene,
            bool(definition and definition.is_reference), axis,
            group_by_sample=scope == 'well')
        self._refresh_plate()
        self.plate_table.select_wells(wells)
        direction = '从左到右' if axis == 'row' else '从上到下'
        self.statusBar().showMessage(
            f'已将 {gene} 应用到 {len(wells)} 个孔，并按{direction}自动编号复孔', 5000)

    def save_mapping_template(self) -> None:
        path, _ = QFileDialog.getSaveFileName(self, '保存板图模板', 'plate_template.json', '模板 (*.json)')
        if path:
            save_template(
                path, self.state.assignments, self.state.genes,
                self.state.config.sample_axis, self.state.config.mapping_scope,
                self.state.config.group_order)

    def load_mapping_template(self) -> None:
        path, _ = QFileDialog.getOpenFileName(self, '载入板图模板', '', '模板 (*.json)')
        if not path:
            return
        try:
            assignments, genes_payload, sample_axis, mapping_scope, group_order = load_template(
                path, include_mapping_scope=True, include_group_order=True)
        except Exception as exc:
            QMessageBox.critical(self, '模板载入失败', str(exc)); return
        self._push_undo(); self.state.assignments = assignments
        self.state.config.sample_axis = sample_axis
        self.state.config.mapping_scope = mapping_scope
        self.state.config.group_order = normalize_group_order(group_order)
        self._load_group_table_from_state(); self._refresh_group_combo()
        if genes_payload:
            self.state.genes = [GeneDefinition(**value) for value in genes_payload]
            self._load_gene_table_from_state(); self._refresh_gene_combo()
        if mapping_scope == 'axis':
            renumber_samples(self.state.assignments, sample_axis)
        self.state.sync_gene_flags(); self._sync_sample_axis_ui(); self._refresh_plate()

    def _load_gene_table_from_state(self) -> None:
        self.gene_table.setRowCount(0)
        for gene in sorted(self.state.genes, key=lambda item: item.order):
            self.add_gene_row(gene.name, gene.is_reference)

    def _build_qc_page(self) -> None:
        page, layout = self._page(
            '技术复孔 QC 与排除值',
            '这里只检查同一样本、同一基因的 Cp 复孔。橙色为可自动建议的唯一异常值；灰色删除线表示已排除，双击 Cp 可手动排除或恢复。')
        controls = QHBoxLayout(); controls.addWidget(QLabel('极差阈值'))
        self.qc_threshold = QDoubleSpinBox(); self.qc_threshold.setRange(0.05, 5.0)
        self.qc_threshold.setSingleStep(0.05); self.qc_threshold.setDecimals(2); self.qc_threshold.setValue(0.5)
        controls.addWidget(self.qc_threshold); controls.addWidget(QLabel('Cp'))
        controls.addSpacing(16); controls.addWidget(QLabel('排除后至少保留'))
        self.qc_minimum = QSpinBox(); self.qc_minimum.setRange(2, 12); self.qc_minimum.setValue(2)
        controls.addWidget(self.qc_minimum); controls.addWidget(QLabel('个有效复孔'))
        refresh = QPushButton('重新检查'); refresh.clicked.connect(self.refresh_qc)
        exclude_all = QPushButton('一键剔除全部建议值'); exclude_all.setObjectName('warningButton')
        exclude_all.clicked.connect(self.exclude_suggestions)
        toggle = QPushButton('排除/恢复选中 Cp'); toggle.clicked.connect(self.toggle_selected_qc)
        undo = QPushButton('撤销'); undo.clicked.connect(self.undo)
        controls.addWidget(refresh); controls.addWidget(exclude_all); controls.addWidget(toggle)
        controls.addWidget(undo); controls.addStretch()
        layout.addLayout(controls)
        self.qc_summary = QLabel('尚未运行 QC'); self.qc_summary.setObjectName('statusCard')
        layout.addWidget(self.qc_summary)
        self.qc_table = QTableWidget(); self.qc_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.qc_table.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        self.qc_table.itemDoubleClicked.connect(lambda *_: self.toggle_selected_qc())
        layout.addWidget(self.qc_table, 1)
        analyze_button = QPushButton('生成分析表与图'); analyze_button.setObjectName('primaryButton')
        analyze_button.clicked.connect(lambda: (self.refresh_analysis(), self.nav.setCurrentRow(3)))
        layout.addWidget(analyze_button)
        self.stack.addWidget(page)

    def refresh_qc(self) -> None:
        if not hasattr(self, 'qc_table'):
            return
        self.state.config.replicate_range_threshold = self.qc_threshold.value()
        self.state.config.minimum_valid_replicates = self.qc_minimum.value()
        rows = build_replicate_qc(
            self.state.wells, self.state.assignments, self.state.config, self.state.excluded_wells)
        max_replicates = max((len(row.wells) for row in rows), default=0)
        headers = ['样本', '组别', '基因'] + [f'Cp{i+1}' for i in range(max_replicates)] + [
            '有效数', '均值', '极差', 'QC 状态', '建议排除值']
        self.qc_table.clear(); self.qc_table.setRowCount(len(rows)); self.qc_table.setColumnCount(len(headers))
        self.qc_table.setHorizontalHeaderLabels(headers)
        by_well = {record.well: record for record in self.state.wells}
        suggestions = 0; warnings = 0
        for row_index, item in enumerate(rows):
            for column, value in enumerate((item.sample, item.group, item.gene)):
                self.qc_table.setItem(row_index, column, QTableWidgetItem(value))
            for replicate, well in enumerate(item.wells):
                record = by_well.get(well); value = record.ct if record else None
                cell = QTableWidgetItem('Undetermined' if value is None else f'{value:.3f}')
                cell.setData(Qt.ItemDataRole.UserRole, well)
                cell.setToolTip(f'孔位 {well}｜双击排除或恢复')
                if well in self.state.excluded_wells:
                    font = cell.font(); font.setStrikeOut(True); cell.setFont(font)
                    cell.setForeground(QColor('#94A3B8')); cell.setBackground(QColor('#E2E8F0'))
                elif well == item.suggested_well:
                    cell.setBackground(QColor('#FED7AA')); cell.setForeground(QColor('#9A3412'))
                self.qc_table.setItem(row_index, 3 + replicate, cell)
            offset = 3 + max_replicates
            if item.n_valid < self.state.config.minimum_valid_replicates:
                status = '有效复孔不足'; warnings += 1
            elif item.range_ct is not None and item.range_ct > self.state.config.replicate_range_threshold:
                status = '偏离过大'; warnings += 1
            else:
                status = '通过'
            if item.suggested_well: suggestions += 1
            values = [str(item.n_valid), display_number(item.mean_ct), display_number(item.range_ct),
                      status, item.suggested_well or '—']
            for column, value in enumerate(values, start=offset):
                cell = QTableWidgetItem(value)
                if status != '通过' and column == offset + 3:
                    cell.setBackground(QColor('#FEE2E2')); cell.setForeground(QColor('#B91C1C'))
                self.qc_table.setItem(row_index, column, cell)
        self.qc_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        self.qc_summary.setText(
            f'共 {len(rows)} 组复孔｜需复核 {warnings} 组｜可一键排除的唯一建议值 {suggestions} 个｜'
            f'当前已排除 {len(self.state.excluded_wells)} 孔')

    def exclude_suggestions(self) -> None:
        rows = build_replicate_qc(
            self.state.wells, self.state.assignments, self.state.config, self.state.excluded_wells)
        wells = [row.suggested_well for row in rows if row.suggested_well]
        if not wells:
            QMessageBox.information(self, '没有建议值', '当前没有满足自动排除规则的唯一异常值。'); return
        self._push_undo(); self.state.exclude(wells, '极差超阈值且移除后恢复合格', 'qc_suggestion')
        self._refresh_plate(); self.refresh_qc()

    def toggle_selected_qc(self) -> None:
        wells = sorted({item.data(Qt.ItemDataRole.UserRole) for item in self.qc_table.selectedItems()
                        if item.data(Qt.ItemDataRole.UserRole)}, key=well_sort_key)
        if not wells:
            return
        self._push_undo()
        to_restore = [well for well in wells if well in self.state.excluded_wells]
        to_exclude = [well for well in wells if well not in self.state.excluded_wells]
        if to_restore: self.state.restore(to_restore, '手动恢复')
        if to_exclude: self.state.exclude(to_exclude, '', 'manual')
        self._refresh_plate(); self.refresh_qc()

    def _build_analysis_page(self) -> None:
        page, layout = self._page(
            'Excel / Prism 式分析',
            '每个靶基因一个工作表标签。表格显示完整计算链，图中柱为组均值、点为独立生物样本；点击 Ct 均值可查看来源孔。')
        self.analysis_notice = QLabel('请先完成板图分配并设置 Ctrl、目标基因和内参。')
        self.analysis_notice.setObjectName('statusCard'); layout.addWidget(self.analysis_notice)
        self.gene_tabs = QTabWidget(); self.gene_tabs.setTabPosition(QTabWidget.TabPosition.South)
        self.gene_tabs.setDocumentMode(True); layout.addWidget(self.gene_tabs, 1)
        self.stack.addWidget(page)

    def _sample_filters(self) -> dict[str, set[str]]:
        return {gene: set(selection.samples) for gene, selection in self.state.plot_selections.items()}

    def refresh_analysis(self) -> None:
        self.state.sync_gene_flags()
        self.result = analyze(
            self.state.wells, self.state.assignments, self.state.config,
            self.state.excluded_wells, self._sample_filters())
        targets = self.state.target_genes()
        for gene in list(self.result_pages):
            if gene not in targets:
                page = self.result_pages.pop(gene)
                self.gene_tabs.removeTab(self.gene_tabs.indexOf(page)); page.deleteLater()
        for gene in targets:
            if gene not in self.result_pages:
                page = GeneResultPage(gene)
                page.filterChanged.connect(self.apply_plot_filter)
                page.errorBarChanged.connect(self.change_error_bar)
                page.exportRequested.connect(self.export_single_svg)
                self.result_pages[gene] = page; self.gene_tabs.addTab(page, gene)
            self.result_pages[gene].update_data(self.state, self.result)
        if not targets:
            self.analysis_notice.setText('没有靶基因：请回到“导入与实验设置”添加非内参基因。')
            return
        no_ctrl = [gene for gene in targets if self.result.calibrator_counts.get(gene, 0)
                   < self.state.config.minimum_calibrator_samples]
        if no_ctrl:
            self.analysis_notice.setText(
                f'以下基因的有效 {self.state.config.calibrator_group} 样本少于 '
                f'{self.state.config.minimum_calibrator_samples} 个，'
                f'Fold change 与统计已停止：{", ".join(no_ctrl)}')
        elif self.result.exploratory_genes:
            self.analysis_notice.setText(
                f'探索性筛选正在生效：{", ".join(sorted(self.result.exploratory_genes))}。'
                'Ctrl 校准均值、Fold change 和统计均已按筛选样本重算。')
        else:
            self.analysis_notice.setText('正式分析：全部未排除生物样本均参与计算。')

    def apply_plot_filter(self, gene: str, samples: list[str]) -> None:
        if not self.result:
            return
        all_samples = [row['sample'] for row in self.result.gene_tables.get(gene, [])]
        self.state.set_plot_selection(gene, samples, all_samples)
        self.refresh_analysis()

    def change_error_bar(self, value: str) -> None:
        if value == self.state.config.error_bar:
            return
        self.state.config.error_bar = value
        self.state.audit.append(self._audit('error_bar', value))
        self.refresh_analysis()

    @staticmethod
    def _audit(action: str, reason: str = ''):
        from .models import AuditAction
        return AuditAction(action=action, reason=reason)

    def _build_export_page(self) -> None:
        page, layout = self._page(
            '导出与工程',
            '主流程只保留工程保存和清理后计算表 Excel；单张 SVG 请在对应基因的分析页导出。原始机器数据始终保持不变。')
        project_card = QFrame(); project_card.setObjectName('card')
        project_layout = QVBoxLayout(project_card); project_layout.addWidget(QLabel('分析工程'))
        row = QHBoxLayout()
        save_button = QPushButton('保存分析工程'); save_button.clicked.connect(self.save_project)
        open_button = QPushButton('打开分析工程'); open_button.clicked.connect(self.open_project)
        row.addWidget(save_button); row.addWidget(open_button); row.addStretch(); project_layout.addLayout(row)
        layout.addWidget(project_card)

        export_card = QFrame(); export_card.setObjectName('card')
        export_layout = QVBoxLayout(export_card); export_layout.addWidget(QLabel('常用导出'))
        row = QHBoxLayout()
        excel = QPushButton('导出清理后计算表 Excel'); excel.setObjectName('primaryButton')
        excel.clicked.connect(self.export_excel)
        row.addWidget(excel); row.addStretch(); export_layout.addLayout(row)
        help_text = QLabel(
            'Excel：每个靶基因一个工作表，ΔCt/ΔΔCt/Fold change 使用可见公式，Ct 批注明确列出保留孔和排除孔；'
            '最后仅附“排除记录”。每个基因的分析页可直接复制图片或单张导出可编辑 SVG。')
        help_text.setWordWrap(True); help_text.setObjectName('hint'); export_layout.addWidget(help_text)
        layout.addWidget(export_card)
        self.export_status = QLabel('尚未导出'); self.export_status.setObjectName('statusCard')
        layout.addWidget(self.export_status); layout.addStretch()
        self.stack.addWidget(page)

    def save_project(self) -> None:
        if not self.state.wells:
            QMessageBox.information(self, '没有数据', '请先导入机器数据。'); return
        path, _ = QFileDialog.getSaveFileName(self, '保存分析工程', 'analysis.qpcrproj', 'qPCR 工程 (*.qpcrproj)')
        if not path:
            return
        if not path.lower().endswith('.qpcrproj'):
            path += '.qpcrproj'
        try:
            save_project(path, self.state)
            self.export_status.setText(f'工程已保存：{path}')
        except Exception as exc:
            QMessageBox.critical(self, '保存失败', str(exc))

    def open_project(self) -> None:
        path, _ = QFileDialog.getOpenFileName(self, '打开分析工程', '', 'qPCR 工程 (*.qpcrproj)')
        if not path:
            return
        try:
            self.state = load_project(path)
        except Exception as exc:
            QMessageBox.critical(self, '打开失败', str(exc)); return
        self.undo_stack.clear(); self.redo_stack.clear(); self.result = None
        self._load_group_table_from_state(); self._refresh_group_combo()
        self._load_gene_table_from_state(); self._refresh_gene_combo()
        self._sync_sample_axis_ui(); self._refresh_plate()
        self.qc_threshold.setValue(self.state.config.replicate_range_threshold)
        self.qc_minimum.setValue(self.state.config.minimum_valid_replicates)
        valid = sum(record.ct is not None for record in self.state.wells)
        self.import_info.setText(
            f'已打开工程｜{self.state.source_name or Path(path).name}｜{len(self.state.wells)} 孔｜'
            f'有效 Cp {valid}｜已排除 {len(self.state.excluded_wells)} 孔')
        self.refresh_qc(); self.refresh_analysis(); self.nav.setCurrentRow(1)
        axis = normalize_sample_axis(self.state.config.sample_axis, 'row')
        scope = normalize_mapping_scope(self.state.config.mapping_scope, 'axis')
        if scope == 'well':
            unassigned_wells = sorted(
                (assignment.well for assignment in self.state.assignments.values()
                 if assignment.role == 'sample' and assignment.gene and not assignment.group),
                key=well_sort_key)
            if unassigned_wells:
                QMessageBox.warning(
                    self,
                    '需要重新选择组别',
                    '工程中以下孔位已有基因但尚未分组：'
                    f'{", ".join(unassigned_wells)}。\n\n'
                    '请在“板图分配”页选中这些孔，并应用组别和样本号。',
                )
            return
        unassigned_axes = {
            sample_axis_key(assignment.well, axis)
            for assignment in self.state.assignments.values()
            if assignment.role == 'sample' and assignment.gene and not assignment.group
        }
        if unassigned_axes:
            ordered_axes = sorted(
                unassigned_axes,
                key=(lambda value: int(value)) if axis == 'column' else None)
            axis_text = ', '.join(ordered_axes)
            unit = '列' if axis == 'column' else '行'
            QMessageBox.warning(
                self,
                '需要重新选择组别',
                f'工程中以下板{unit}尚未使用内置组别：'
                f'{axis_text}。\n\n'
                '这通常来自旧版自由文本组别无法自动映射，或工程保存时尚未完成分组。'
                    f'请在“板图分配”页选中相应{unit}，并重新应用当前实验组别。',
            )

    def export_excel(self) -> None:
        if not self.state.wells:
            QMessageBox.information(self, '没有数据', '请先导入并完成分析。'); return
        path, _ = QFileDialog.getSaveFileName(self, '导出清理后计算表', 'qPCR_cleaned.xlsx', 'Excel (*.xlsx)')
        if not path:
            return
        if not path.lower().endswith('.xlsx'):
            path += '.xlsx'
        try:
            export_clean_excel(path, self.state)
            self.export_status.setText(f'清理后计算表已导出：{path}')
        except Exception as exc:
            QMessageBox.critical(self, '导出失败', str(exc))

    def export_svgs(self) -> None:
        if not self.result:
            self.refresh_analysis()
        folder = QFileDialog.getExistingDirectory(self, '选择 SVG 输出文件夹')
        if not folder:
            return
        try:
            paths = export_all_svg(folder, self.state, self.result, self.state.config.error_bar)
            self.export_status.setText(f'已导出 {len(paths)} 张 SVG：{folder}')
        except Exception as exc:
            QMessageBox.critical(self, '导出失败', str(exc))

    def export_single_svg(self, gene: str) -> None:
        if not self.result:
            return
        suffix = '_filtered' if gene in self.result.exploratory_genes else ''
        path, _ = QFileDialog.getSaveFileName(self, '导出 SVG', f'{gene}{suffix}.svg', 'SVG (*.svg)')
        if not path:
            return
        if not path.lower().endswith('.svg'):
            path += '.svg'
        export_figure(
            path, self.result, gene, self.state.config.error_bar,
            group_order=self.state.config.group_order,
            calibrator_group=self.state.config.calibrator_group)

    def _apply_style(self) -> None:
        self.setStyleSheet('''
            QMainWindow, QWidget { background: #F8FAFC; color: #0F172A; }
            #nav { background: #0F172A; color: #CBD5E1; border: 0; padding: 18px 8px; }
            #nav::item { padding: 14px 12px; margin: 3px; border-radius: 7px; }
            #nav::item:selected { background: #1D4ED8; color: white; }
            #heading { font-size: 24px; font-weight: 700; margin: 14px 12px 2px 12px; }
            #subtitle { color: #64748B; margin: 0 12px 10px 12px; }
            #dropZone { border: 2px dashed #60A5FA; border-radius: 12px; background: #EFF6FF;
                        color: #1D4ED8; font-size: 15px; margin: 4px 12px; }
            #card, #statusCard { background: white; border: 1px solid #E2E8F0; border-radius: 9px;
                                padding: 10px; margin: 4px 12px; }
            #selectionStatus { background: #DBEAFE; color: #1D4ED8; border: 1px solid #60A5FA;
                               border-radius: 6px; padding: 8px; font-weight: 700; }
            #axisBadge { background: #E0F2FE; color: #075985; border: 1px solid #7DD3FC;
                         border-radius: 6px; padding: 6px 10px; font-weight: 700; }
            #processBar { background: #EFF6FF; color: #1D4ED8; padding: 8px; border-radius: 6px; }
            #officialBadge { background: #DCFCE7; color: #166534; padding: 7px; border-radius: 6px; }
            #exploratoryBadge { background: #FEF3C7; color: #92400E; padding: 7px; border-radius: 6px; }
            #hint { color: #64748B; font-size: 9pt; }
            QPushButton { background: white; border: 1px solid #CBD5E1; border-radius: 8px;
                          padding: 8px 13px; min-height: 18px; }
            QPushButton:hover { background: #EFF6FF; border-color: #60A5FA; color: #1D4ED8; }
            QPushButton:pressed { background: #DBEAFE; border-color: #2563EB; }
            QPushButton:disabled { background: #F1F5F9; color: #94A3B8; border-color: #E2E8F0; }
            #primaryButton { background: #2563EB; color: white; border-color: #2563EB;
                             font-weight: 700; padding: 9px 16px; }
            #primaryButton:hover { background: #1D4ED8; color: white; border-color: #1D4ED8; }
            #secondaryButton { background: #EFF6FF; color: #1D4ED8; border-color: #93C5FD;
                               font-weight: 600; }
            #secondaryButton:hover { background: #DBEAFE; border-color: #60A5FA; }
            #ghostButton { background: transparent; color: #475569; border-color: #E2E8F0; }
            #warningButton { background: #FFF7ED; color: #C2410C; border-color: #FDBA74; }
            #warningButton:hover { background: #FFEDD5; border-color: #FB923C; }
            QComboBox, QSpinBox, QDoubleSpinBox { background: white; border: 1px solid #CBD5E1;
                                                 border-radius: 7px; padding: 6px 9px; min-height: 20px; }
            QComboBox:hover, QSpinBox:hover, QDoubleSpinBox:hover { border-color: #60A5FA; }
            QTableWidget { background: white; alternate-background-color: #F8FAFC;
                           gridline-color: #E2E8F0; border: 1px solid #CBD5E1; }
            QHeaderView::section { background: #E2E8F0; color: #334155; padding: 6px;
                                   border: 0; border-right: 1px solid #CBD5E1; font-weight: 700; }
            QTabBar::tab { background: #E2E8F0; padding: 8px 18px; border: 1px solid #CBD5E1; }
            QTabBar::tab:selected { background: white; color: #1D4ED8; font-weight: 700; }
            QScrollBar:vertical { background: white; width: 13px; border: 1px solid #E2E8F0; }
            QScrollBar:horizontal { background: white; height: 13px; border: 1px solid #E2E8F0; }
            QScrollBar::handle:vertical, QScrollBar::handle:horizontal {
                background: #CBD5E1; min-height: 34px; min-width: 34px; border-radius: 4px; margin: 1px; }
            QScrollBar::handle:vertical:hover, QScrollBar::handle:horizontal:hover { background: #94A3B8; }
        ''')


def main() -> int:
    app = QApplication.instance() or QApplication(sys.argv)
    app.setApplicationName('qPCR 分析助手')
    window = MainWindow(); window.showMaximized()
    return app.exec()


def run() -> None:
    raise SystemExit(main())


if __name__ == '__main__':
    raise SystemExit(main())
