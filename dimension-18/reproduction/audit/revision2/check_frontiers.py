"""Audit joins between stored ADE certificate domains, including every surviving label.

Usage: python3 -B check_frontiers.py PATH/TO/EL18_58_search
Python standard library. Numerical and structural validity of each elimination
witness is checked by the original stage verifier. This program checks transport
between stages without executing the original programs or trusting cached joins.
"""
import json
import sys
from collections import Counter
from pathlib import Path


def require(condition, reason):
    if not condition:
        raise ValueError(reason)


def distinct(items, where):
    keys = list(items)
    require(len(keys) == len(set(keys)), f'duplicate IDs in {where}')
    return set(keys)


def positive(table):
    return {name: ids for name, ids in table.items() if ids}


def compare(left, right, stage):
    require(set(left) == set(right), f'incomplete type keys at {stage}')
    for name, expected in left.items():
        require(expected == right[name], f'wrong labels at {stage} / {name}')


def label_pairs(table):
    return {name: distinct(map(tuple, pairs), name) for name, pairs in table.items()}


def main(folder):
    def get(name):
        return json.loads((folder / 'root_profiles' / f'{name}.json').read_text())
    front9 = get('frontier9')
    current = get('frontier_current')
    prior = get('prior_lower')
    pure16 = '.'.join(['A1'] * 16)
    source = get('plane16')
    live = positive({
        name: {(int(a), int(q)) for a, rec in row.items()
               for q, cut in rec.items() if cut is None}
        for name, row in source.items()
    })
    print('rank 16 recovered types and (anchor, area) IDs')
    print('plane16', len(live), sum(map(len, live.values())))
    for name in ('plane16_next', 'plane16_root', 'refine16', 'caps16',
                 'capcycle16', 'splits16', 'kernel16', 'finish16'):
        stage = get(name)
        require(set(live) == set(stage), f'{name}: input type mismatch')
        for family, record in stage.items():
            steps = [record['cuts']] if 'cuts' in record else record['stages']
            removed = distinct((tuple(map(int, key.split(':')[-2:])) for step in steps
                                for key in step), f'{name} / {family} / cuts')
            after = distinct(map(tuple, record['live']), f'{name} / {family} / live')
            require(removed.isdisjoint(after) and removed | after == live[family],
                    f'{name} / {family}: incomplete elimination partition')
            live[family] = after
        live = positive(live)
        print(name, len(live), sum(map(len, live.values())))
        if name == 'refine16':
            live = {k: ids for k, ids in live.items() if k in prior['16'] and k != pure16}
            compare(live, label_pairs(front9['rank16']), 'frontier9')
            require(set(front9['types']['16']) == set(live), 'type list at frontier9')
            print('frontier9', len(live), sum(map(len, live.values())))
    compare(live, label_pairs(current['rank16']), 'frontier_current')
    require(set(current['types']['16']) == set(live), 'type list at frontier_current')
    labelled = {f'{family}:{a}:{q}' for family, items in live.items() for a,q in items}
    for table in ('word16', 'orbit16'):
        require(set(get(table)) == labelled, f'terminal coverage {table}')
    # A mutation that keeps every per-family count and the number of all
    # retained IDs fixed must still fail this comparison.
    mutated = {name: labels.copy() for name, labels in live.items()}
    key = min(mutated); a,q = min(mutated[key])
    mutated[key].remove((a,q))
    fake = max(j for _,j in live[key]) + 1
    mutated[key].add((a,fake))
    try:
        compare(live, mutated, 'deliberate count-preserving substitution')
    except ValueError:
        pass
    else:
        raise AssertionError('mutation was accepted')
    print('21 terminal word and orbit IDs verified; count-preserving substitution rejected')

    totals = Counter()
    code_sectors = {}
    short_codes = get('dual_spaces')
    full_codes = get('full_dual')
    for rank in range(6, 16):
        key = str(rank)
        table = get('lower_' + key)
        allowed = set(prior[key]) - {'.'.join(['A1'] * rank)}
        active = {k: distinct(v['live'], k) for k,v in table.items() if k in allowed and v['live']}
        for cached in (front9, current):
            require(distinct(cached['types'][key], 'rank '+key) == set(active),
                    'initial type intersection rank '+key)
        for basename in ('pack','energy','join','block','pencil','knap','purebranch','twonormal'):
            if not active:
                break
            stage = get(basename + '_' + key)
            require(set(stage) == set(active), basename+' / '+key)
            for family, record in stage.items():
                removed = distinct((cut[0] for cut in record['cuts']), family)
                after = distinct(record['live'], family)
                require(after.isdisjoint(removed) and after | removed == active[family],
                        f'{basename} / {rank} / {family}: incomplete label partition')
                totals[basename] += len(removed)
                active[family] = after
            active = positive(active)
        bucket = 'binary' if rank <= 12 else 'ternary'
        if active:
            stage = get(bucket+'_'+key)
            require(set(stage) == set(active), bucket+' '+key)
            for family,r in stage.items():
                cuts = distinct((cut[0] for cut in r['cuts']), family)
                after = distinct(r['live'], family)
                require(after.isdisjoint(cuts) and after|cuts == active[family], family)
                totals[bucket] += len(cuts)
                active[family] = after
        print('mixed rank', rank, 'post-code types', len(positive(active)),
              'live anchor IDs', sum(map(len, active.values())))
        active = positive(active)
        if not active:
            continue
        local_codes = get('dual_' + key)
        require(set(active) == set(local_codes), 'initial code families ' + key)
        field = 2 if rank < 13 else 3
        fixed_file = folder / 'root_profiles' / f'fixed_{rank}.json'
        fixed = get('fixed_' + key) if fixed_file.exists() else {}
        axes_file = folder / 'root_profiles' / f'fixedaxes_{rank}.json'
        axes = get('fixedaxes_' + key) if axes_file.exists() else {}
        needed_fixed = set()
        for family, anchors in active.items():
            factors = family.split('.')
            tag = f'{field}:{factors.count("A" + str(field - 1))}:{18-rank}'
            row = local_codes[family]
            require(len(row) == len(short_codes[tag]), 'initial code-list size ' + family)
            if field == 2:
                need = {str(i) for i,v in enumerate(row) if v is None}
            else:
                # The larger-dimensional ternary codes are explicitly enumerated;
                # they receive every branch of their own table.
                need = set(map(str, range(len(full_codes[tag]))))
            if not need:
                continue
            needed_fixed.add(family)
            require(set(fixed[family]) == need, 'fixed code IDs ' + family)
            live_codes = {}
            for j,entry in fixed[family].items():
                if entry.get('initial') is not None:
                    require(not entry['live'], 'impossible initial code listed as live')
                    continue
                cuts = distinct((v[0] for v in entry['cuts']), family)
                after = distinct(entry['live'], family)
                require(cuts.isdisjoint(after) and cuts|after == anchors, 'fixed code anchor partition')
                if after:
                    live_codes[j] = after
            later = axes.get(family, {})
            require(set(later) == set(live_codes), 'fixed-axis code IDs')
            for j,entry in later.items():
                cuts = distinct((v[0] for v in entry['cuts']), family)
                after = distinct(entry['live'], family)
                require(cuts.isdisjoint(after) and cuts|after == live_codes[j], 'fixed-axis anchor partition')
                if after:
                    code_sectors[family+':'+j] = after
        require(set(fixed) == needed_fixed, 'redundant or missing fixed code family')
    high = get('high_normals')
    reads = get('finish_roots')
    require(set(high) == set(code_sectors) == set(reads), 'high-support code cases')
    total_orbits = total_alive = total_assignments = 0
    for key,record in high.items():
        total_orbits += record['stats'].get('orbits', 0)
        starts = {(tuple(pts),mult) for pts,mult,unused in record['live']}
        ends = {(tuple(row[0]),row[1]) for row in reads[key]}
        require(starts == ends and len(starts) == len(record['live']) == len(reads[key]),
                'missing or duplicated support orbit ' + key)
        require(all(not row[5] for row in reads[key]), 'surviving normal Gram')
        total_alive += len(starts)
        total_assignments += sum(row[4] for row in reads[key])
    print('final code cases, tested support orbits, surviving support orbits, tested assignments:',
          len(code_sectors), total_orbits, total_alive, total_assignments)
    print('checked labeled mixed-rank elimination events', dict(totals))

    pairs = get('pairs17'); codim = get('codim1')
    expected = {f'{name}:{q}' for name,row in codim.items() for q,cut in row.items() if cut is None}
    require(set(pairs) == expected, 'rank17 codim -> pairs')
    saturation = get('saturation17')
    require(set(saturation) == {name for name,v in pairs.items() if v['surviving']},
            'rank17 pairs -> saturation')
    print('rank 17 joins',len(codim),'types',len(pairs),'type/height IDs',len(saturation),'survivors')


if __name__ == '__main__':
    main(Path(sys.argv[1]))
