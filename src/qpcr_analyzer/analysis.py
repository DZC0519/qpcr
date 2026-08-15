from __future__ import annotations
from collections import defaultdict
from itertools import combinations
from math import sqrt
import numpy as np
from scipy import stats
from .models import (
    GROUP_ORDER, AnalysisConfig, AnalysisResult, PlateAssignment, WellRecord,
    normalize_group_order,
)
from .qc import build_replicate_qc, run_qc

def bh_adjust(pvalues):
    values = np.asarray(list(pvalues), float)
    if not len(values): return []
    order = np.argsort(values); ranked = values[order]
    adjusted = np.minimum.accumulate((ranked*len(values)/np.arange(1,len(values)+1))[::-1])[::-1]
    result = np.empty_like(adjusted); result[order] = np.minimum(adjusted, 1)
    return result.tolist()

def welch_anova(groups):
    groups = [np.asarray(g,float) for g in groups if len(g)>=2]; k = len(groups)
    if k < 2: return np.nan, np.nan, np.nan, np.nan
    n = np.array([len(g) for g in groups],float)
    means = np.array([g.mean() for g in groups]); variances = np.array([g.var(ddof=1) for g in groups])
    if np.any(variances <= 0): return np.nan, np.nan, np.nan, np.nan
    weights = n/variances; mean = np.sum(weights*means)/np.sum(weights)
    term = np.sum(((1-weights/np.sum(weights))**2)/(n-1))
    f_value = (np.sum(weights*(means-mean)**2)/(k-1))/(1+2*(k-2)*term/(k**2-1))
    df1, df2 = k-1, (k**2-1)/(3*term)
    return float(f_value), float(stats.f.sf(f_value,df1,df2)), float(df1), float(df2)

def games_howell(groups):
    rows, names, k = [], [n for n,v in groups.items() if len(v)>=2], len(groups)
    for a,b in combinations(names,2):
        x,y = np.asarray(groups[a],float), np.asarray(groups[b],float)
        vx,vy = x.var(ddof=1), y.var(ddof=1); se2 = vx/len(x)+vy/len(y)
        if se2 <= 0: pvalue=df=q=np.nan
        else:
            df = se2**2/((vx/len(x))**2/(len(x)-1)+(vy/len(y))**2/(len(y)-1))
            q = abs(x.mean()-y.mean())/sqrt(se2/2)
            pvalue = float(stats.studentized_range.sf(q,k,df))
        rows.append({'comparison': f'{a} vs {b}', 'group_a': a, 'group_b': b,
            'mean_diff_delta_ct': float(x.mean()-y.mean()), 'statistic': float(q),
            'df': float(df), 'p_value': pvalue, 'test': 'Games-Howell'})
    return rows

def analyze(wells: tuple[WellRecord,...], assignments: dict[str,PlateAssignment],
            config: AnalysisConfig, excluded_wells: dict[str,str] | None = None,
            sample_filters: dict[str, set[str]] | None = None):
    excluded = excluded_wells or {}
    filters = sample_filters or {}
    qc = run_qc(wells, assignments, config, excluded)
    replicate_qc = build_replicate_qc(wells, assignments, config, excluded)
    grouped = defaultdict(list)
    for record in wells:
        a = assignments.get(record.well)
        if a and a.role == 'sample' and a.sample and a.group and a.gene:
            grouped[(a.sample, a.group, a.gene, a.is_reference)].append(record)
    summaries = []
    for (sample,group,gene,is_reference), records in grouped.items():
        records = sorted(records, key=lambda r: r.well)
        valid = [r for r in records if r.ct is not None and r.include and r.well not in excluded]
        values = np.array([r.ct for r in valid],float)
        enough = len(values) >= config.minimum_valid_replicates
        summaries.append({'sample':sample,'group':group,'gene':gene,'is_reference':is_reference,
            'wells':', '.join(r.well for r in records), 'well_list':[r.well for r in records],
            'retained_wells':[r.well for r in valid],
            'excluded_wells':[r.well for r in records if r.well in excluded],
            'retained_ct':[float(r.ct) for r in valid],
            'excluded_ct':[float(r.ct) for r in records if r.well in excluded and r.ct is not None],
            'n_total':len(records),'n_valid':len(valid),
            'mean_ct':float(values.mean()) if enough else np.nan,
            'sd_ct':float(values.std(ddof=1)) if len(values)>=2 else np.nan,
            'range_ct':float(np.ptp(values)) if len(values)>=2 else np.nan})

    reference_genes = {r['gene'] for r in summaries if r['is_reference']}
    references = defaultdict(dict)
    for row in summaries:
        if row['is_reference'] and np.isfinite(row['mean_ct']):
            references[row['sample']][row['gene']] = row['mean_ct']
    delta_rows = []
    for row in summaries:
        if row['is_reference']:
            continue
        refs = references.get(row['sample'], {})
        complete_refs = bool(reference_genes) and reference_genes.issubset(refs)
        reference_ct = float(np.mean(list(refs.values()))) if complete_refs else np.nan
        ok = np.isfinite(row['mean_ct']) and complete_refs
        selected = row['sample'] in filters.get(row['gene'], {row['sample']})
        delta = row['mean_ct'] - reference_ct if ok else np.nan
        delta_rows.append({'sample':row['sample'],'group':row['group'],'gene':row['gene'],
            'target_ct':row['mean_ct'],'reference_ct':reference_ct,
            'reference_genes':', '.join(sorted(reference_genes)),
            'delta_ct':float(delta) if ok else np.nan,'computable':bool(ok),
            'included_for_analysis': selected,
            'target_wells':row['retained_wells'], 'target_excluded':row['excluded_wells'],
            'reference_details':refs})

    present_groups = {r['group'] for r in delta_rows if r.get('group')}
    configured_groups = normalize_group_order(getattr(config, 'group_order', None))
    group_order = configured_groups + sorted(present_groups - set(configured_groups))

    calibrator = config.calibrator_group or 'Ctrl'
    genes = sorted({r['gene'] for r in delta_rows})
    calibrator_means, calibrator_counts = {}, {}
    for gene in genes:
        vals = [r['delta_ct'] for r in delta_rows if r['gene']==gene and r['group']==calibrator
                and r['included_for_analysis'] and np.isfinite(r['delta_ct'])]
        calibrator_counts[gene] = len(vals)
        if len(vals) >= config.minimum_calibrator_samples:
            calibrator_means[gene] = float(np.mean(vals))
    expression_rows, base = [], 1 + config.amplification_efficiency
    for row in delta_rows:
        cal = calibrator_means.get(row['gene'], np.nan)
        usable = row['included_for_analysis'] and np.isfinite(row['delta_ct']) and np.isfinite(cal)
        ddct = row['delta_ct'] - cal if usable else np.nan
        expression_rows.append({**row,'calibrator_group':calibrator,
            'calibrator_mean_delta_ct':cal,'delta_delta_ct':float(ddct) if np.isfinite(ddct) else np.nan,
            'fold_change':float(base**(-ddct)) if np.isfinite(ddct) else np.nan})
    statistics_rows = []
    for gene in genes:
        names = [g for g in group_order if any(r['gene'] == gene and r['group'] == g
                                                for r in delta_rows)]
        by_group = {g:np.array([r['delta_ct'] for r in delta_rows if r['gene']==gene
            and r['group']==g and r['included_for_analysis']
            and np.isfinite(r['delta_ct'])],float) for g in names}
        eligible = {k:v for k,v in by_group.items() if len(v)>=2}
        if len(eligible)==2:
            a,b = list(eligible); test = stats.ttest_ind(eligible[a],eligible[b],equal_var=False)
            statistics_rows.append({'gene':gene,'comparison':f'{a} vs {b}',
                'group_a':a,'group_b':b,'test':'Welch t-test','statistic':float(test.statistic),
                'df':float(test.df),'p_value':float(test.pvalue)})
        elif len(eligible)>=3:
            f_value,pvalue,df1,df2 = welch_anova(list(eligible.values()))
            statistics_rows.append({'gene':gene,'comparison':'\u603b\u4f53','group_a':'',
                'group_b':'','test':'Welch ANOVA','statistic':f_value,
                'df':f'{df1:.2f}, {df2:.2f}','p_value':pvalue})
            statistics_rows.extend({'gene':gene,**row} for row in games_howell(eligible))
    finite = [r for r in statistics_rows if np.isfinite(r.get('p_value',np.nan))]
    for row,p_adj in zip(finite,bh_adjust([r['p_value'] for r in finite])):
        row['p_adjust_bh'],row['significant'] = p_adj,p_adj<config.alpha
    for row in statistics_rows:
        row.setdefault('p_adjust_bh',np.nan); row.setdefault('significant',False)
    for gene in genes:
        for group in group_order:
            vals = np.array([r['delta_ct'] for r in delta_rows if r['gene']==gene and
                r['group']==group and r['included_for_analysis'] and np.isfinite(r['delta_ct'])],float)
            cal = calibrator_means.get(gene,np.nan)
            if not len(vals) or not np.isfinite(cal): continue
            mean = vals.mean(); low=high=np.nan
            sd = float(vals.std(ddof=1)) if len(vals)>=2 else np.nan
            sem = float(stats.sem(vals)) if len(vals)>=2 else np.nan
            if len(vals)>=2:
                if np.isfinite(sem) and sem > 0:
                    low,high = stats.t.interval(.95,len(vals)-1,loc=mean,scale=sem)
                else:
                    low=high=float(mean)
            expression_rows.append({'sample':'\u7ec4\u6c47\u603b','group':group,'gene':gene,
                'delta_ct':float(mean),'delta_delta_ct':float(mean-cal),
                'fold_change':float(base**(-(mean-cal))),
                'fold_change_ci_low':float(base**(-(high-cal))) if np.isfinite(high) else np.nan,
                'fold_change_ci_high':float(base**(-(low-cal))) if np.isfinite(low) else np.nan,
                'delta_ct_sd':sd, 'delta_ct_sem':sem,
                'n':len(vals),'computable':True,'calibrator_group':calibrator,
                'calibrator_mean_delta_ct':cal})

    gene_tables = {}
    group_rank = {group:index for index,group in enumerate(group_order)}
    for gene in genes:
        rows = [r for r in expression_rows if r['gene']==gene and r['sample']!='组汇总']
        rows.sort(key=lambda r: (group_rank.get(r['group'], 99), r['sample']))
        gene_tables[gene] = rows
    exploratory = set()
    for gene, selected in filters.items():
        all_samples = {r['sample'] for r in delta_rows if r['gene']==gene}
        if set(selected) != all_samples:
            exploratory.add(gene)
    return AnalysisResult(
        replicate_summary=summaries, delta_ct=delta_rows,
        relative_expression=expression_rows, statistics=statistics_rows,
        qc_flags=qc, replicate_qc=replicate_qc, gene_tables=gene_tables,
        calibrator_counts=calibrator_counts, exploratory_genes=exploratory)
