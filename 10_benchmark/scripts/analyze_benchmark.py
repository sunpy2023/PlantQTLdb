#!/usr/bin/env python3
"""Recompute benchmark summaries and paired bootstrap CIs from published TSVs."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

from validate_benchmark import FULL, KG, SCORE_COLUMNS, SYSTEM, TOTAL, validate


def save(frame, path):
    frame.to_csv(path, sep='\t', index=False, float_format='%.9f', lineterminator='\n')


def analyze(root, output, resamples=10000, seed=20261115, ablation_seed=20261116):
    validate(root)
    output.mkdir(parents=True, exist_ok=True)
    scores = pd.read_csv(root / 'scores/primary_scores.tsv', sep='\t')
    questions = pd.read_csv(root / 'questions.tsv', sep='\t').set_index('Question ID')
    scores = scores.merge(questions[['Species group', 'Task category']].reset_index(),
                          on='Question ID', validate='many_to_one')
    tasks = pd.read_csv(root / 'task_categories.tsv', sep='\t')['Task category'].tolist()
    matrix = scores.pivot(index='Question ID', columns=SYSTEM, values=TOTAL).sort_index()
    meta = questions.loc[matrix.index]
    values = matrix.to_numpy(float)
    rng = np.random.default_rng(seed)
    boot = np.zeros((resamples, len(matrix.columns)))
    # Keep each question's scores paired across systems; preserve all 12 stratum sizes.
    for indices in meta.groupby(['Species group', 'Task category'], sort=True).indices.values():
        boot += values[rng.choice(indices, size=(resamples, len(indices)), replace=True)].sum(axis=1) / len(matrix)
    ci = np.percentile(boot, [2.5, 97.5], axis=0)
    overall = pd.DataFrame({SYSTEM: matrix.columns, 'Questions': len(matrix),
                            'Mean total score': values.mean(axis=0),
                            'CI lower': ci[0], 'CI upper': ci[1]})
    save(overall.sort_values('Mean total score', ascending=False), output / 'overall_scores.tsv')
    task_means = scores.groupby([SYSTEM, 'Task category'], sort=True)[TOTAL].agg(['size', 'mean']).reset_index()
    save(task_means.rename(columns={'size': 'Questions', 'mean': 'Mean total score'}), output / 'task_scores.tsv')
    dims = scores.groupby(SYSTEM, sort=True)[list(SCORE_COLUMNS)].mean().reset_index()
    save(dims, output / 'dimension_scores.tsv')
    fi = matrix.columns.get_loc(FULL)
    paired = []
    for system in matrix.columns:
        if system in (FULL, KG):
            continue
        idx = matrix.columns.get_loc(system)
        lo, hi = np.percentile(boot[:, fi] - boot[:, idx], [2.5, 97.5])
        delta = values[:, fi] - values[:, idx]
        paired.append({'Comparison': f'Full Assistant minus {system}', 'Questions': len(delta),
                       'Mean paired difference': delta.mean(), 'CI lower': lo, 'CI upper': hi,
                       'Wins': int((delta > 0).sum()), 'Ties': int((delta == 0).sum()),
                       'Losses': int((delta < 0).sum())})
    save(pd.DataFrame(paired), output / 'paired_baseline_differences.tsv')

    # Reproduce Figure 3C's separate seed and task-first resampling order.
    rng = np.random.default_rng(ablation_seed)
    both = scores[scores[SYSTEM].isin([KG, FULL])].pivot(
        index=['Question ID', 'Species group', 'Task category'], columns=SYSTEM, values=TOTAL
    )[[KG, FULL]].reset_index()
    ablation = []
    pooled_boot = np.zeros((resamples, 2))

    def summarize(task, subset, sampled):
        means = subset[[KG, FULL]].mean().to_numpy()
        bounds = np.percentile(sampled, [2.5, 97.5], axis=0)
        diff_ci = np.percentile(sampled[:, 1] - sampled[:, 0], [2.5, 97.5])
        return {'task': task, 'n': len(subset), 'kg_mean': means[0], 'kg_low': bounds[0, 0],
                'kg_high': bounds[1, 0], 'full_mean': means[1], 'full_low': bounds[0, 1],
                'full_high': bounds[1, 1], 'difference': means[1] - means[0],
                'difference_low': diff_ci[0], 'difference_high': diff_ci[1]}

    for task in tasks:
        subset = both[both['Task category'] == task]
        sampled = np.zeros((resamples, 2))
        for _, group in subset.groupby('Species group', sort=True):
            vals = group[[KG, FULL]].to_numpy(float)
            sampled += vals[rng.integers(0, len(group), size=(resamples, len(group)))].sum(axis=1)
        pooled_boot += sampled
        ablation.append(summarize(task, subset, sampled / len(subset)))
    save(pd.DataFrame([summarize('Overall', both, pooled_boot / len(both)), *ablation]), output / 'resource_ablation.tsv')

    for name, round_col in [('repeated_scoring', 'Scoring round'),
                            ('repeated_generation_scores', 'Generation round')]:
        frame = pd.read_csv(root / f'scores/{name}.tsv', sep='\t')
        summaries = []
        rounds = []
        for system, group in frame.groupby(SYSTEM, sort=True):
            mat = group.pivot(index='Question ID', columns=round_col, values=TOTAL)[[1, 2, 3]]
            ranges = mat.max(axis=1) - mat.min(axis=1)
            summaries.append({SYSTEM: system, 'Questions': len(mat),
                              'Mean score range': ranges.mean(), 'Maximum score range': ranges.max(),
                              'Mean sample SD': mat.std(axis=1, ddof=1).mean(),
                              'Responses or questions with range <=10': int((ranges <= 10).sum()),
                              'Identical scores across rounds': int((ranges == 0).sum())})
            for round_id, mean in mat.mean(axis=0).items():
                rounds.append({SYSTEM: system, round_col: round_id, 'Mean total score': mean})
        round_frame = pd.DataFrame(rounds)
        round_frame['Rank among all nine configurations'] = round_frame.groupby(round_col)['Mean total score'].rank(method='min', ascending=False).astype(int)
        plotted = round_frame[round_frame[SYSTEM] != KG].copy()
        plotted['Rank among eight plotted systems'] = plotted.groupby(round_col)['Mean total score'].rank(method='min', ascending=False).astype(int)
        save(pd.DataFrame(summaries), output / f'{name}_summary.tsv')
        save(round_frame, output / f'{name}_round_means.tsv')
        save(plotted, output / f'{name}_plotted_round_means.tsv')
    inputs = [root / p for p in ('questions.tsv', 'task_categories.tsv', 'scoring_rubric.tsv',
              'responses/primary_responses.tsv', 'responses/repeated_generation_responses.tsv',
              'scores/primary_scores.tsv', 'scores/repeated_scoring.tsv', 'scores/repeated_generation_scores.tsv')]
    metadata = {'bootstrap_resamples': resamples, 'overall_and_baseline_seed': seed,
                'resource_ablation_seed': ablation_seed, 'sample_sd_ddof': 1,
                'python': sys.version.split()[0], 'numpy': np.__version__, 'pandas': pd.__version__,
                'input_sha256': {p.relative_to(root).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}}
    (output / 'analysis_metadata.json').write_text(json.dumps(metadata, indent=2) + '\n', encoding='utf-8')
    print(f'Wrote benchmark summaries to {output}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-dir', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--resamples', type=int, default=10000)
    parser.add_argument('--seed', type=int, default=20261115)
    parser.add_argument('--ablation-seed', type=int, default=20261116)
    args = parser.parse_args()
    if args.resamples < 1 or args.seed < 0 or args.ablation_seed < 0:
        parser.error('Resamples must be positive and seeds nonnegative')
    root, output = args.data_dir.resolve(), args.output_dir.resolve()
    protected = (root / 'scores', root / 'responses', root / 'scripts', root / 'tests')
    if output == root or output in root.parents or any(p == output or p in output.parents for p in protected):
        parser.error('Output directory must not contain or replace benchmark inputs/scripts')
    analyze(root, output, args.resamples, args.seed, args.ablation_seed)


if __name__ == '__main__':
    main()
