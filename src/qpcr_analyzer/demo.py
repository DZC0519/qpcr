from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path

from openpyxl import Workbook

from .exports import export_clean_excel
from .models import AnalysisConfig, GeneDefinition, PlateAssignment, ProjectState, WellRecord


def generate_demo(folder: str | Path) -> tuple[Path, Path, Path]:
    folder = Path(folder); folder.mkdir(parents=True, exist_ok=True)
    machine_path = folder/'synthetic_qpcr_example.xlsx'
    calculation_path = folder/'synthetic_qpcr_calculation.xlsx'
    template_path = folder/'synthetic_plate_template.json'
    workbook = Workbook(); sheet = workbook.active; sheet.title = 'Synthetic_Run'
    sheet.append(['Synthetic qPCR export - generated data only'])
    sheet.append(['Include', 'Color', 'Pos', 'Name', 'Cp', 'Concentration', 'Standard', 'Status'])

    group_rows = [('A', 'Ctrl', 1), ('B', 'Ctrl', 2), ('C', 'Model', 1), ('D', 'Model', 2),
                  ('E', 'Low', 1), ('F', 'Low', 2), ('G', 'Med', 1), ('H', 'Med', 2),
                  ('I', 'High', 1), ('J', 'High', 2), ('K', 'Pos', 1), ('L', 'Pos', 2)]
    effects_a = {'Ctrl':0.0, 'Model':2.0, 'Low':1.5, 'Med':1.0, 'High':0.3, 'Pos':2.6}
    effects_b = {'Ctrl':0.0, 'Model':-1.1, 'Low':-.8, 'Med':-.4, 'High':-.1, 'Pos':-1.5}
    genes = [GeneDefinition('Target_A', False, 0), GeneDefinition('Target_B', False, 1),
             GeneDefinition('REF', True, 2)]
    offsets = (-.08, .03, .05)
    assignments = {}; records = []
    row_info = {row:(group, number) for row, group, number in group_rows}
    for row in 'ABCDEFGHIJKLMNOP':
        for column in range(1, 25):
            well = f'{row}{column}'
            if row in row_info and column <= 9:
                group, number = row_info[row]; sample = f'{group}-{number}'
                gene_index = (column - 1)//3; gene = genes[gene_index].name
                base = (20.2-effects_a[group], 22.5-effects_b[group], 17.0)[gene_index]
                ct = base + (number-1)*.16 + offsets[(column-1)%3]
                if well == 'C3': ct += 1.15
                ct = round(ct, 3)
                assignment = PlateAssignment(well, sample, group, gene, gene == 'REF',
                                             (column-1)%3+1, plate_row=row)
                sheet.append([True, gene_index+1, well, sample, ct, None, 0, ''])
                records.append(WellRecord(well, ct, machine_name=sample, source_sheet=sheet.title))
            else:
                assignment = PlateAssignment(well, role='empty', plate_row=row)
                sheet.append([True, 0, well, 'Blank', None, None, 0, ''])
                records.append(WellRecord(well, None, machine_name='Blank', source_sheet=sheet.title))
            assignments[well] = assignment
    workbook.save(machine_path)
    state = ProjectState(source_path=str(machine_path), source_sheet=sheet.title, header_row=2,
        wells=tuple(records), assignments=assignments, config=AnalysisConfig(384), genes=genes)
    export_clean_excel(calculation_path, state)
    template_path.write_text(json.dumps({'schema_version':2,
        'assignments':{key:asdict(value) for key, value in assignments.items()},
        'genes':[asdict(value) for value in genes]}, ensure_ascii=False, indent=2), encoding='utf-8')
    return machine_path, calculation_path, template_path


if __name__ == '__main__':
    generate_demo(Path(__file__).resolve().parents[2]/'assets')
