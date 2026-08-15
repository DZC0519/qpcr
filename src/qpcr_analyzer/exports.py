from __future__ import annotations

import re
from dataclasses import asdict
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
import pandas as pd

from .analysis import analyze
from .models import GROUP_ORDER, AnalysisResult, ProjectState, normalize_group_order

CHINESE_FONT = Path('C:/Windows/Fonts/msyh.ttc')
FONT_STACK = ['Arial', 'Helvetica', 'DejaVu Sans']
if CHINESE_FONT.exists():
    font_manager.fontManager.addfont(str(CHINESE_FONT))
    FONT_STACK.insert(2, font_manager.FontProperties(fname=str(CHINESE_FONT)).get_name())
    plt.rcParams['axes.unicode_minus'] = False

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = FONT_STACK

GROUP_COLORS = {
    'Ctrl':'#4C78A8', 'Model':'#E45756', 'Low':'#72B7B2',
    'Med':'#F2CF5B', 'High':'#B279A2', 'Pos':'#59A14F',
}

# Keep the plate-map colors unchanged while using a quieter publication palette
# for exported figures.
PLOT_GROUP_COLORS = {
    'Ctrl': '#C8D5E6', 'Model': '#E8C2BF', 'Low': '#C9DECF',
    'Med': '#F0DAB4', 'High': '#D7C8E2', 'Pos': '#D0D2D4',
}

CUSTOM_GROUP_COLORS = [
    '#7B9ACC', '#D48782', '#70A99C', '#D7B65A', '#9A80A8', '#6FA36D',
    '#B58B68', '#6D9AA8', '#C47D9D', '#81965C',
]
CUSTOM_PLOT_COLORS = [
    '#D2DDEC', '#EACBC8', '#D2E5D8', '#F0E1B9', '#DFD4E8', '#D3E1D1',
    '#E4D4C5', '#D1E0E5', '#E8D2DC', '#DCE3C8',
]


def group_color(group: str, group_order: list[str] | tuple[str, ...] | None = None) -> str:
    if group in GROUP_COLORS:
        return GROUP_COLORS[group]
    order = normalize_group_order(group_order)
    index = order.index(group) if group in order else len(order)
    return CUSTOM_GROUP_COLORS[index % len(CUSTOM_GROUP_COLORS)]


def plot_group_color(group: str, group_order: list[str] | tuple[str, ...] | None = None) -> str:
    if group in PLOT_GROUP_COLORS:
        return PLOT_GROUP_COLORS[group]
    order = normalize_group_order(group_order)
    index = order.index(group) if group in order else len(order)
    return CUSTOM_PLOT_COLORS[index % len(CUSTOM_PLOT_COLORS)]

plt.rcParams.update({
    'svg.fonttype': 'none',
    'pdf.fonttype': 42,
    'axes.linewidth': 1.0,
})


def _safe_float(value):
    try:
        return float(value) if np.isfinite(value) else None
    except (TypeError, ValueError):
        return None


def _safe_sheet_name(name: str, used: set[str]) -> str:
    base = re.sub(r'[\\/*?:\[\]]', '_', name).strip() or 'Gene'
    base = base[:31]
    candidate, index = base, 2
    while candidate.lower() in used:
        suffix = f'_{index}'
        candidate = base[:31-len(suffix)] + suffix
        index += 1
    used.add(candidate.lower())
    return candidate


def _safe_filename(name: str) -> str:
    value = re.sub(r'[<>:"/\\|?*]+', '_', name).strip(' .')
    return value or 'gene'


def _group_summary(result: AnalysisResult, gene: str) -> dict[str, dict]:
    return {row['group']: row for row in result.relative_expression
            if row.get('gene') == gene and row.get('sample') == '组汇总'}


def significance_stars(p_value) -> str:
    value = _safe_float(p_value)
    if value is None:
        return ''
    if value < .0001:
        return '****'
    if value < .001:
        return '***'
    if value < .01:
        return '**'
    if value < .05:
        return '*'
    return ''


def _reference_significance(
    result: AnalysisResult, gene: str, groups: list[str],
    group_order: list[str] | tuple[str, ...] | None = None,
    calibrator_group: str = 'Ctrl',
) -> tuple[str | None, list[tuple[str, str]]]:
    reference = (
        'Model' if 'Model' in groups
        else calibrator_group if calibrator_group in groups
        else 'Ctrl' if 'Ctrl' in groups else None
    )
    if reference is None:
        return None, []
    comparisons = []
    for row in result.statistics:
        if row.get('gene') != gene:
            continue
        group_a, group_b = row.get('group_a'), row.get('group_b')
        if reference not in {group_a, group_b}:
            continue
        other = group_b if group_a == reference else group_a
        if other not in groups:
            continue
        stars = significance_stars(row.get('p_value'))
        if stars:
            comparisons.append((other, stars))
    order = {group: index for index, group in enumerate(normalize_group_order(group_order))}
    reference_index = groups.index(reference)
    comparisons.sort(key=lambda item: (
        abs(groups.index(item[0]) - reference_index), order.get(item[0], 999)))
    return reference, comparisons


def make_gene_figure(
    result: AnalysisResult, gene: str, error_bar: str = 'SEM',
    group_order: list[str] | tuple[str, ...] | None = None,
    calibrator_group: str = 'Ctrl',
):
    samples = [row for row in result.relative_expression
               if row.get('gene') == gene and row.get('sample') != '组汇总'
               and row.get('included_for_analysis') and _safe_float(row.get('fold_change')) is not None]
    summaries = _group_summary(result, gene)
    configured_groups = normalize_group_order(group_order)
    present_groups = {row['group'] for row in samples if row.get('group')}
    ordered_groups = configured_groups + sorted(present_groups - set(configured_groups))
    groups = [group for group in ordered_groups if any(row['group'] == group for row in samples)]
    figure_width = max(3.7, min(5.4, 2.6 + len(groups) * .55))
    fig, ax = plt.subplots(figsize=(figure_width, 4.2), constrained_layout=True)
    if not groups:
        ax.text(.5, .5, f'没有可计算数据\n请检查 {calibrator_group}、内参和有效复孔数量',
                ha='center', va='center', transform=ax.transAxes, color='#64748B')
        ax.set_axis_off()
        return fig

    x = np.arange(len(groups), dtype=float)
    heights, low_errors, high_errors = [], [], []
    for group in groups:
        summary = summaries.get(group, {})
        height = _safe_float(summary.get('fold_change'))
        if height is None:
            values = [row['fold_change'] for row in samples if row['group'] == group]
            height = float(np.mean(values)) if values else np.nan
        heights.append(height)
        low = high = 0.0
        if error_bar == '95% CI':
            ci_low = _safe_float(summary.get('fold_change_ci_low'))
            ci_high = _safe_float(summary.get('fold_change_ci_high'))
            if ci_low is not None and ci_high is not None:
                low, high = max(0.0, height-ci_low), max(0.0, ci_high-height)
        elif error_bar in {'SD', 'SEM'}:
            key = 'delta_ct_sd' if error_bar == 'SD' else 'delta_ct_sem'
            spread = _safe_float(summary.get(key))
            if spread is not None:
                lower_fc, upper_fc = height * 2**(-spread), height * 2**spread
                low, high = height-lower_fc, upper_fc-height
        low_errors.append(low)
        high_errors.append(high)

    colors = [plot_group_color(group, ordered_groups) for group in groups]
    ax.bar(x, heights, width=.50, color=colors,
           edgecolor='#4A4A4A', linewidth=.8)
    if error_bar != '无':
        ax.errorbar(x, heights, yerr=np.array([low_errors, high_errors]), fmt='none',
                    ecolor='#4A4A4A', elinewidth=1.05, capsize=3.5,
                    capthick=1.05, zorder=3)
    rng = np.random.default_rng(20260724)
    for index, group in enumerate(groups):
        values = [row['fold_change'] for row in samples if row['group'] == group]
        jitter = rng.uniform(-.13, .13, len(values)) if len(values) > 1 else np.zeros(len(values))
        ax.scatter(index+jitter, values, s=34, facecolor='white', edgecolor='#4A4A4A',
                   linewidth=1.0, zorder=4)
    reference, comparisons = _reference_significance(
        result, gene, groups, ordered_groups, calibrator_group)
    if comparisons and reference in groups:
        group_tops = []
        for index, group in enumerate(groups):
            values = [float(row['fold_change']) for row in samples
                      if row['group'] == group and np.isfinite(float(row['fold_change']))]
            local_top = heights[index] + high_errors[index]
            if values:
                local_top = max(local_top, max(values))
            group_tops.append(local_top)
        base_top = max(group_tops, default=1.0)
        padding = max(base_top * .065, .06)
        annotation_tops = []
        for other, stars in comparisons:
            target = 'Model' if reference == 'Model' and other == 'Ctrl' else other
            target_index = groups.index(target)
            y = group_tops[target_index] + padding * .58
            ax.plot([target_index - .20, target_index + .20], [y, y],
                    color='#4A4A4A', linewidth=1.0, solid_capstyle='butt',
                    clip_on=False)
            ax.text(target_index, y + padding * .06, stars,
                    ha='center', va='bottom', fontsize=11, weight='bold',
                    color='#1F1F1F')
            annotation_tops.append(y + padding * .78)
        current_bottom, current_top = ax.get_ylim()
        needed_top = max(annotation_tops, default=base_top + padding)
        ax.set_ylim(0, max(current_top, needed_top))
    else:
        current_bottom, current_top = ax.get_ylim()
        ax.set_ylim(0, current_top)

    rotate_labels = max((len(group) for group in groups), default=0) > 8
    ax.set_xticks(x, groups, rotation=38 if rotate_labels else 0,
                  ha='right' if rotate_labels else 'center')
    ax.set_ylabel(f'Relative {gene} mRNA' if gene else 'Relative mRNA',
                  fontsize=10.5)
    ax.set_xlim(-.55, len(groups) - .45)
    ax.tick_params(axis='both', direction='out', length=4, width=1.0,
                   labelsize=9.5, color='#1F1F1F')
    ax.spines[['top','right']].set_visible(False)
    ax.spines['left'].set_color('#1F1F1F')
    ax.spines['bottom'].set_color('#1F1F1F')
    ax.spines['left'].set_linewidth(1.1)
    ax.spines['bottom'].set_linewidth(1.1)
    ax.grid(False)
    return fig


def make_fold_change_figure(
    result: AnalysisResult, group_order: list[str] | tuple[str, ...] | None = None,
):
    genes = sorted({row.get('gene') for row in result.relative_expression if row.get('gene')})
    return (make_gene_figure(result, genes[0], group_order=group_order)
            if genes else make_gene_figure(result, '', group_order=group_order))


def export_figure(path: str | Path, result: AnalysisResult, gene: str | None = None,
                  error_bar: str = 'SEM', dpi: int = 300,
                  group_order: list[str] | tuple[str, ...] | None = None,
                  calibrator_group: str = 'Ctrl') -> None:
    path = Path(path)
    if gene is None:
        gene = next(iter(result.gene_tables), '')
    figure = make_gene_figure(result, gene, error_bar, group_order, calibrator_group)
    figure.savefig(path, dpi=dpi, bbox_inches='tight')
    plt.close(figure)


def export_all_svg(folder: str | Path, state: ProjectState, result: AnalysisResult,
                   error_bar: str = 'SEM') -> list[Path]:
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    written = []
    for gene in state.target_genes() or list(result.gene_tables):
        suffix = '_filtered' if gene in result.exploratory_genes else ''
        path = folder / f'{_safe_filename(gene)}{suffix}.svg'
        export_figure(
            path, result, gene, error_bar, group_order=state.config.group_order,
            calibrator_group=state.config.calibrator_group)
        written.append(path)
    return written


def _source_comment(summary: dict | None) -> str:
    if not summary:
        return '无可用复孔信息'
    kept = ', '.join(f'{w}={v:.3f}' for w, v in zip(
        summary.get('retained_wells', []), summary.get('retained_ct', []))) or '无'
    removed = ', '.join(f'{w}={v:.3f}' for w, v in zip(
        summary.get('excluded_wells', []), summary.get('excluded_ct', []))) or '无'
    return f'参与均值：{kept}\n已排除：{removed}'


def export_clean_excel(path: str | Path, state: ProjectState) -> None:
    sample_filters = {
        gene: set(selection.samples)
        for gene, selection in state.plot_selections.items()
    }
    official = analyze(
        state.wells, state.assignments, state.config, state.excluded_wells,
        sample_filters=sample_filters,
    )
    summary_map = {(r['sample'], r['gene']): r for r in official.replicate_summary}
    used: set[str] = set()
    with pd.ExcelWriter(path, engine='xlsxwriter') as writer:
        workbook = writer.book
        header = workbook.add_format({'bold':True, 'bg_color':'#DCE6F1', 'border':1,
                                      'align':'center', 'valign':'vcenter', 'text_wrap':True})
        cell = workbook.add_format({'border':1, 'align':'center', 'num_format':'0.000'})
        text_cell = workbook.add_format({'border':1, 'align':'center'})
        formula = workbook.add_format({'border':1, 'align':'center', 'num_format':'0.0000'})
        group_order = normalize_group_order(getattr(state.config, 'group_order', None))
        group_formats = {g:workbook.add_format({'border':1, 'align':'center', 'valign':'vcenter',
            'bg_color':group_color(g, group_order), 'font_color':'#FFFFFF', 'bold':True})
            for g in group_order}
        refs = state.reference_genes()
        for gene in state.target_genes() or list(official.gene_tables):
            rows = [
                row for row in official.gene_tables.get(gene, [])
                if row.get('included_for_analysis')
            ]
            sheet_name = _safe_sheet_name(gene, used)
            worksheet = workbook.add_worksheet(sheet_name)
            writer.sheets[sheet_name] = worksheet
            headers = ['组别', '样本', f'{gene} Ct', f'内参 Ct ({", ".join(refs)})',
                       'ΔCt', f'{state.config.calibrator_group}组平均ΔCt',
                       'ΔΔCt', 'Fold change=2^(-ΔΔCt)']
            worksheet.write_row(0, 0, headers, header)
            worksheet.set_row(0, 34)
            worksheet.set_column('A:B', 13)
            worksheet.set_column('C:H', 20)
            ctrl_excel_rows = []
            group_ranges: dict[str, list[int]] = {}
            for index, row in enumerate(rows, start=1):
                excel_row = index + 1
                group_ranges.setdefault(row['group'], []).append(index)
                if row['group'] == state.config.calibrator_group and row.get('computable'):
                    ctrl_excel_rows.append(excel_row)
                worksheet.write(index, 1, row['sample'], text_cell)
                target = _safe_float(row.get('target_ct'))
                reference = _safe_float(row.get('reference_ct'))
                if target is None:
                    worksheet.write(index, 2, '不可计算', text_cell)
                else:
                    worksheet.write_number(index, 2, target, cell)
                if reference is None:
                    worksheet.write(index, 3, '不可计算', text_cell)
                else:
                    worksheet.write_number(index, 3, reference, cell)
                target_summary = summary_map.get((row['sample'], gene))
                worksheet.write_comment(index, 2, _source_comment(target_summary))
                ref_lines = [_source_comment(summary_map.get((row['sample'], ref))) for ref in refs]
                worksheet.write_comment(index, 3, '\n\n'.join(
                    f'{ref}:\n{text}' for ref, text in zip(refs, ref_lines)))
                if target is not None and reference is not None:
                    worksheet.write_formula(index, 4, f'=C{excel_row}-D{excel_row}', formula,
                                            _safe_float(row.get('delta_ct')) or 0)
                else:
                    worksheet.write(index, 4, '不可计算', text_cell)
            ctrl_ready = (len(ctrl_excel_rows) >= state.config.minimum_calibrator_samples
                          and official.calibrator_counts.get(gene, 0)
                          >= state.config.minimum_calibrator_samples)
            if ctrl_ready:
                first_ctrl, last_ctrl = min(ctrl_excel_rows), max(ctrl_excel_rows)
                ctrl_formula = f'=AVERAGE(E{first_ctrl}:E{last_ctrl})'
                ctrl_value = next((_safe_float(r.get('calibrator_mean_delta_ct'))
                                   for r in rows if _safe_float(r.get('calibrator_mean_delta_ct')) is not None), None)
                for index, row in enumerate(rows, start=1):
                    excel_row = index + 1
                    worksheet.write_formula(index, 5, ctrl_formula, formula, ctrl_value)
                    if row.get('computable') and ctrl_value is not None:
                        worksheet.write_formula(index, 6, f'=E{excel_row}-F{excel_row}', formula,
                                                _safe_float(row.get('delta_delta_ct')) or 0)
                        worksheet.write_formula(index, 7, f'=POWER(2,-G{excel_row})', formula,
                                                _safe_float(row.get('fold_change')) or 0)
                    else:
                        worksheet.write(index, 6, '不可计算', text_cell)
                        worksheet.write(index, 7, '不可计算', text_cell)
            else:
                for index in range(1, len(rows) + 1):
                    worksheet.write(index, 5, '不可计算', text_cell)
                    worksheet.write(index, 6, '不可计算', text_cell)
                    worksheet.write(index, 7, '不可计算', text_cell)
            for group, indices in group_ranges.items():
                start, end = min(indices), max(indices)
                fmt = group_formats.get(group, text_cell)
                if start == end:
                    worksheet.write(start, 0, group, fmt)
                else:
                    worksheet.merge_range(start, 0, end, 0, group, fmt)
            worksheet.freeze_panes(1, 2)

        log_name = _safe_sheet_name('排除记录', used)
        log = workbook.add_worksheet(log_name)
        writer.sheets[log_name] = log
        log_headers = ['孔位', 'Cp', '样本', '组别', '基因', '来源', '理由', '当前状态', '时间']
        log.write_row(0, 0, log_headers, header)
        by_well = {record.well: record for record in state.wells}
        for index, record in enumerate(state.exclusions, start=1):
            assignment = state.assignments.get(record.well)
            well_record = by_well.get(record.well)
            log.write_row(index, 0, [
                record.well, well_record.ct if well_record else '',
                assignment.sample if assignment else '', assignment.group if assignment else '',
                assignment.gene if assignment else '',
                'QC建议' if record.source == 'qc_suggestion' else ('手动' if record.source == 'manual' else '旧版迁移'),
                record.reason, '已排除' if record.active and record.well in state.excluded_wells else '已恢复',
                record.timestamp,
            ], text_cell)
        log.set_column('A:I', 18)


def export_excel(path: str | Path, state: ProjectState, result: AnalysisResult | None = None) -> None:
    export_clean_excel(path, state)


def export_csv_folder(folder: str | Path, result: AnalysisResult) -> None:
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(result.delta_ct).to_csv(folder/'delta_ct.csv', index=False, encoding='utf-8-sig')
    pd.DataFrame(result.relative_expression).to_csv(
        folder/'fold_change.csv', index=False, encoding='utf-8-sig')


def export_pdf(path: str | Path, state: ProjectState, result: AnalysisResult) -> None:
    """Compatibility export retained for v1 callers; it is no longer shown in the main UI."""
    from matplotlib.backends.backend_pdf import PdfPages
    with PdfPages(path) as pdf:
        for gene in state.target_genes() or list(result.gene_tables):
            figure = make_gene_figure(
                result, gene, state.config.error_bar, state.config.group_order,
                state.config.calibrator_group)
            pdf.savefig(figure, bbox_inches='tight')
            plt.close(figure)
