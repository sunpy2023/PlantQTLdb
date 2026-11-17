#!/usr/bin/env python3
"""Validate published benchmark records without calling any model or external service."""
import argparse
from collections import Counter
import csv
from pathlib import Path

DIMENSIONS = (
    'Factual accuracy', 'Evidence completeness', 'Source traceability',
    'Hallucination control and calibration', 'Synthesis quality',
)
SCORE_COLUMNS = tuple(f'{name} (0–20)' for name in DIMENSIONS)
TOTAL = 'Total score (0–100)'
SYSTEM = 'System configuration'
FULL = 'PlantQTLdb Assistant (DeepSeek V4.1 Flash)'
KG = 'PlantQTLdb Assistant (DeepSeek V4.1 Flash; literature KG only)'
SYSTEMS = {FULL, KG, 'GPT-5.6 Sol', 'GLM-5.3', 'DeepSeek V4 Pro',
           'DeepSeek V4.1 Flash', 'Kimi K2.5', 'MiniMax M3', 'Qwen-3.8-Flash'}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read_tsv(path):
    with Path(path).open(encoding='utf-8-sig', newline='') as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        require(reader.fieldnames and len(set(reader.fieldnames)) == len(reader.fieldnames),
                f'{path}: missing or duplicate column names')
        rows = list(reader)
    require(bool(rows), f'{path}: no records')
    require(all(None not in r and all(v is not None for v in r.values()) for r in rows),
            f'{path}: inconsistent number of fields')
    return rows


def index_rows(rows, columns, label):
    require(all(all(c in r and r[c].strip() for c in columns) for r in rows),
            f'{label}: missing key column/value {columns}')
    keys = [tuple(r[c] for c in columns) for r in rows]
    require(len(set(keys)) == len(keys), f'{label}: duplicate key {columns}')
    return dict(zip(keys, rows))


def validate_scores(rows, responses, label):
    index_rows(rows, ('Evaluation ID',), label)
    index_rows(rows, ('Question ID', SYSTEM, 'Generation round', 'Scoring round'), label)
    for r in rows:
        require(all(c in r and r[c].isdigit() and 0 <= int(r[c]) <= 20
                    for c in SCORE_COLUMNS), f'{label}: invalid dimension score')
        require(TOTAL in r and r[TOTAL].isdigit() and
                int(r[TOTAL]) == sum(int(r[c]) for c in SCORE_COLUMNS),
                f'{label}: total is not the sum of five dimension scores')
        require(all(r.get(f'{c} rationale', '').strip() for c in DIMENSIONS),
                f'{label}: missing scoring rationale')
        require(r.get('Judge model') == 'GPT-5.6 Sol', f'{label}: unexpected judge')
        rid = r['Response ID']
        require((rid,) in responses, f'{label}: unknown response ID {rid}')
        answer = responses[(rid,)]
        require(all(r[k] == answer[k] for k in
                    ('Question ID', SYSTEM, 'Generation round')),
                f'{label}: response metadata mismatch for {rid}')
        require(r['Evaluation ID'] == rid + '__s' + r['Scoring round'],
                f'{label}: evaluation ID disagrees with response/scoring round')


def validate(root):
    root = Path(root)
    questions = read_tsv(root / 'questions.tsv')
    qindex = index_rows(questions, ('Question ID',), 'questions')
    require(len(questions) == 300, 'Expected 300 questions')
    require(sorted(Counter(r['Species group'] for r in questions).values()) == [150, 150],
            'Questions must contain 150 questions per species group')
    require(sorted(Counter(r['Task category'] for r in questions).values()) == [50] * 6,
            'Questions must contain six categories of 50 questions')
    require(set(Counter((r['Species group'], r['Task category']) for r in questions).values()) == {25},
            'Expected 25 questions per species/task stratum')
    require(all(r['Repeat subset'] in ('yes', 'no') for r in questions), 'Invalid repeat flag')
    subset = {r['Question ID'] for r in questions if r['Repeat subset'] == 'yes'}
    require(len(subset) == 30, 'Expected a fixed 30-question repeat subset')
    repeated_questions = [r for r in questions if r['Question ID'] in subset]
    require(sorted(Counter(r['Species group'] for r in repeated_questions).values()) == [15, 15],
            'Repeat subset must contain 15 questions per species group')
    require(sorted(Counter(r['Task category'] for r in repeated_questions).values()) == [5] * 6,
            'Repeat subset must contain five questions per task')
    task_rows = read_tsv(root / 'task_categories.tsv')
    require(set(r['Task category'] for r in task_rows) == set(r['Task category'] for r in questions),
            'Task-label mapping disagrees with question categories')
    rubric = read_tsv(root / 'scoring_rubric.tsv')
    require(Counter(r['Dimension'] for r in rubric) == Counter({d: 5 for d in DIMENSIONS}),
            'Expected five score bands for each of five rubric dimensions')
    require(all(r['Definition'] and r['Criterion'] for r in rubric), 'Incomplete rubric row')

    primary_answers = read_tsv(root / 'responses/primary_responses.tsv')
    extra_answers = read_tsv(root / 'responses/repeated_generation_responses.tsv')
    answers = primary_answers + extra_answers
    answer_index = index_rows(answers, ('Response ID',), 'responses')
    index_rows(answers, ('Question ID', SYSTEM, 'Generation round'), 'responses')
    require(len(primary_answers) == 2700 and len(extra_answers) == 540,
            'Expected 2700 primary and 540 additional responses')
    for r in answers:
        require((r['Question ID'],) in qindex, 'Unknown question in responses')
        require(r[SYSTEM] in SYSTEMS and r['Answer'].strip(), 'Unknown system or empty answer')
        q = qindex[(r['Question ID'],)]
        require(all(r[k] == q[k] for k in ('Species group', 'Task category', 'Question')),
                'Question text or metadata mismatch')
    require({(r['Question ID'], r[SYSTEM], r['Generation round']) for r in primary_answers} ==
            {(q[0], s, '1') for q in qindex for s in SYSTEMS}, 'Primary answer coverage mismatch')
    require({(r['Question ID'], r[SYSTEM], r['Generation round']) for r in extra_answers} ==
            {(q, s, str(g)) for q in subset for s in SYSTEMS for g in (2, 3)},
            'Repeated-generation answer coverage mismatch')

    primary = read_tsv(root / 'scores/primary_scores.tsv')
    scoring = read_tsv(root / 'scores/repeated_scoring.tsv')
    generation = read_tsv(root / 'scores/repeated_generation_scores.tsv')
    for rows, label in ((primary, 'primary scores'), (scoring, 'repeated scoring'),
                        (generation, 'repeated generation scores')):
        validate_scores(rows, answer_index, label)
    require(len(primary) == 2700 and len(scoring) == len(generation) == 810,
            'Unexpected score record counts')
    key = lambda r: (r['Question ID'], r[SYSTEM], r['Generation round'], r['Scoring round'])
    require({key(r) for r in primary} == {(q[0], s, '1', '1') for q in qindex for s in SYSTEMS},
            'Primary score coverage mismatch')
    require({key(r) for r in scoring} == {(q, s, '1', str(n)) for q in subset for s in SYSTEMS for n in (1, 2, 3)},
            'Repeated-scoring coverage mismatch')
    require({key(r) for r in generation} == {(q, s, str(n), '1') for q in subset for s in SYSTEMS for n in (1, 2, 3)},
            'Repeated-generation score coverage mismatch')
    pindex = index_rows(primary, ('Evaluation ID',), 'primary scores')
    for rows in (scoring, generation):
        for r in rows:
            if r['Generation round'] == r['Scoring round'] == '1':
                require(r == pindex[(r['Evaluation ID'],)], 'Shared round-1 record differs from primary evaluation')
    result = {'questions': 300, 'systems': 9, 'primary_responses': len(primary_answers),
              'additional_responses': len(extra_answers), 'primary_scores': len(primary),
              'repeated_scoring_records': len(scoring), 'repeated_generation_records': len(generation),
              'shared_round_1_records': 270}
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-dir', type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    import json
    try:
        print(json.dumps(validate(args.data_dir), indent=2))
    except (ValueError, KeyError, OSError) as error:
        parser.exit(1, f'Validation failed: {error}\n')


if __name__ == '__main__':
    main()
