#!/usr/bin/env python3
"""Run bounded, serial scientific workload measurements without changing source data."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import resource
import shutil
import signal
import subprocess
import time

GIB = 1024 ** 3


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def available_memory():
    match = re.search(r'^MemAvailable:\s+(\d+)', Path('/proc/meminfo').read_text(), re.M)
    return int(match.group(1)) * 1024


def write_json(path, data):
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(data, indent=2) + '\n')
    temporary.replace(path)


def change(text, key, value):
    expression = rf'(?m)^\s*(?://\s*)?{re.escape(key)}:.*$'
    result, count = re.subn(expression, f'{key}: {value}', text)
    if count != 1:
        raise RuntimeError(f'Expected one setting for {key}, found {count}')
    return result


def prepare(source, case, work):
    work.mkdir()
    assets = source / 'Input_files'
    (work / 'Input_files').symlink_to(assets, target_is_directory=True)
    text = (assets / '1_base/instructions.txt').read_text()
    updates = {'main_repository': '.', 'nghb_count': str(case.get('neighbors', 100))}
    if case.get('radius_scale', 1.0) != 1.0:
        updates['domain'] = 'density-box.txt'
        with (assets / '1_base/example_box.txt').open() as incoming, (work / 'density-box.txt').open('w') as outgoing:
            data = False
            for line in incoming:
                if data and line.strip():
                    fields = line.split()
                    if float(fields[1]) > 0:
                        fields[1] = f'{float(fields[1]) * case["radius_scale"]:.10f}'
                    line = '\t'.join(fields) + '\n'
                outgoing.write(line)
                if line.startswith('Index'):
                    data = True
    if case.get('water_tables', False) or case.get('features', False):
        updates['springs'] = 'springs-wt12.txt'
        updates['surf_wat_table'] = 'Input_files/1_base/example_watertable_surf1.txt Input_files/1_base/example_watertable_surf2.txt'
        shutil.copy2(source / 'KarstNSim/tests/regression/fixtures/springs_wt12.txt', work / 'springs-wt12.txt')
        (work / 'connectivity_matrix.txt').write_text('1\t0\n1\t0\n1\t0\n0\t1\n0\t1\n')
    else:
        shutil.copy2(assets / '1_base/connectivity_matrix.txt', work / 'connectivity_matrix.txt')
    if case.get('features', False):
        updates.update({
            'use_waypoints': 'true', 'use_previous_networks': 'true',
            'previous_networks': 'Input_files/3_amplification/base_0_karst.txt Input_files/3_amplification/polyphasic_0_karst.txt',
            'fraction_old_karst_perm': '0.5', 'use_deadend_points': 'true',
            'nb_deadend_points': '20', 'max_distance_of_deadend_pts': '150',
            'use_cycle_amplification': 'true', 'max_distance_amplification': '300',
            'nb_cycles': '30', 'use_noise': 'true', 'use_noise_on_all': 'false',
            'simulate_sections': 'true',
        })
    if case.get('graph_export', False):
        updates.update(create_nghb_graph='true', create_nghb_graph_property='true')
    for key, value in updates.items():
        text = change(text, key, value)
    (work / 'instructions.txt').write_text(text)
    return updates


def limits():
    resource.setrlimit(resource.RLIMIT_AS, (5 * GIB, 5 * GIB))
    resource.setrlimit(resource.RLIMIT_FSIZE, (4 * GIB, 4 * GIB))
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    os.nice(5)


def stop_child(process):
    children_file = Path(f'/proc/{process.pid}/task/{process.pid}/children')
    children = children_file.read_text().split() if children_file.exists() else []
    for child in children:
        try:
            os.kill(int(child), signal.SIGKILL)
        except ProcessLookupError:
            pass
    if not children and process.poll() is None:
        process.kill()


def run(source, binary, case, out):
    work = out / case['name']
    updates = prepare(source, case, work)
    result = dict(case=case, parameter_updates=updates, binary=str(binary),
                  binary_sha256=digest(binary), instructions_sha256=digest(work / 'instructions.txt'),
                  address_space_limit_bytes=5 * GIB, file_size_limit_bytes=4 * GIB)
    if available_memory() < 3 * GIB:
        result.update(status='not_run_host_memory', exit_code=None)
        return result
    started = time.monotonic()
    terminated = None
    print(json.dumps({'event': 'started', 'case': case['name']}), flush=True)
    with (work / 'run.log').open('w') as log:
        process = subprocess.Popen(['/usr/bin/time', '-v', '-o', str(work / 'time.txt'),
                                    str(binary), 'instructions.txt'], cwd=work,
                                   stdout=log, stderr=subprocess.STDOUT, preexec_fn=limits)
        while process.poll() is None:
            elapsed = time.monotonic() - started
            free = available_memory()
            write_json(out / 'status.json', {'case': case['name'], 'elapsed_seconds': elapsed,
                                           'host_available_bytes': free, 'pid': process.pid})
            if terminated is None and (elapsed > 600 or free < 2 * GIB):
                terminated = 'wall_time_limit' if elapsed > 600 else 'host_memory_guard'
                stop_child(process)
            time.sleep(1)
    runner_wall = time.monotonic() - started
    report = (work / 'time.txt').read_text() if (work / 'time.txt').exists() else ''
    elapsed_match = re.search(r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s*(\S+)', report)
    wall = sum(float(value) * 60 ** i for i, value in enumerate(reversed(elapsed_match.group(1).split(':')))) if elapsed_match else runner_wall
    log = (work / 'run.log').read_text(errors='replace')
    rss_match = re.search(r'Maximum resident set size \(kbytes\):\s*(\d+)', report)
    exports = {}
    points = []
    for path in sorted((work / 'outputs').glob('*')):
        if path.is_file() and not path.name.startswith('karstnsim_console_'):
            exports[path.name] = {'sha256': digest(path), 'bytes': path.stat().st_size}
            if re.fullmatch(r'base_\d+_pts\.txt', path.name):
                with path.open() as stream:
                    points.append(sum(1 for _ in stream) - 1)
    complete = process.returncode == 0 and 'Simulation completed successfully' in log
    status = 'completed' if complete else terminated or 'failed'
    if not complete and 'bad_alloc' in log:
        status = 'allocation_failed_under_5_gib_address_space_limit'
    if not complete and process.returncode == 153:
        status = 'file_size_limit'
    generated = re.search(r'Using (\d+) automatically sampled points', log)
    result.update(status=status, exit_code=process.returncode, wall_seconds=wall,
                  runner_wall_seconds=runner_wall,
                  peak_rss_bytes=int(rss_match.group(1)) * 1024 if rss_match else None,
                  support_points=points, automatically_sampled_points=int(generated.group(1)) if generated else None,
                  output_bytes=sum(item['bytes'] for item in exports.values()), exports=exports,
                  stage_timings=[{'line': line, 'seconds': float(seconds)} for line, seconds in
                      re.findall(r'(?m)^([^\n]*?)\(([\d.]+) s\)\s*$', log)])
    write_json(work / 'metrics.json', result)
    print(json.dumps({'event': 'finished', 'case': case['name'], 'status': status,
                      'peak_rss_bytes': result['peak_rss_bytes'], 'wall_seconds': wall,
                      'support_points': points}), flush=True)
    return result


def summarize_comparisons(results, out):
    indexed = {result['case']['name']: result for result in results}
    original = indexed['feature_parity_reference']
    optimized = indexed['feature_parity_optimized']
    both = original['status'] == optimized['status'] == 'completed'
    comparison = {'both_completed': both,
                  'exports_identical': both and original['exports'] == optimized['exports']}
    write_json(out / 'feature-parity.json', comparison)
    extra = {}
    if 'export_control_medium' in indexed and 'export_costs_medium' in indexed:
        control, exported = indexed['export_control_medium'], indexed['export_costs_medium']
        complete = control['status'] == exported['status'] == 'completed'
        extra = {'both_completed': complete,
                 'common_exports_identical': complete and all(
                     exported['exports'].get(name) == value for name, value in control['exports'].items()),
                 'additional_exports': sorted(set(exported.get('exports', {})) - set(control.get('exports', {})))}
        write_json(out / 'export-parity.json', extra)
    write_json(out / 'status.json', {'complete': True, 'feature_parity': comparison, 'export_parity': extra})
    print(json.dumps({'event': 'complete', 'feature_parity': comparison, 'export_parity': extra}), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    source, binary, reference, out = [p.resolve() for p in (args.source, args.binary, args.reference, args.out)]
    out.mkdir(parents=True, exist_ok=False)
    cases = [
        {'name': 'base_control'},
        {'name': 'neighbors_200', 'neighbors': 200},
        {'name': 'two_water_tables', 'water_tables': True},
        {'name': 'combined_features_full', 'features': True},
        {'name': 'denser_080', 'radius_scale': 0.8},
        {'name': 'denser_060', 'radius_scale': 0.6},
        {'name': 'denser_features_080', 'radius_scale': 0.8, 'neighbors': 150, 'features': True},
        {'name': 'full_graph_export_costs', 'graph_export': True},
        {'name': 'feature_parity_reference', 'radius_scale': 1.5, 'features': True, 'reference': True},
        {'name': 'feature_parity_optimized', 'radius_scale': 1.5, 'features': True},
        {'name': 'export_control_medium', 'radius_scale': 2.0},
        {'name': 'export_costs_medium', 'radius_scale': 2.0, 'graph_export': True},
    ]
    results = []
    manifest = {'source_commit': subprocess.check_output(['git', '-C', str(source), 'rev-parse', 'HEAD'], text=True).strip(),
                'cases': cases, 'limits': {'address_space_bytes': 5 * GIB, 'single_file_bytes': 4 * GIB,
                                         'wall_seconds': 600, 'host_memory_reserve_bytes': 2 * GIB}}
    write_json(out / 'manifest.json', manifest)
    for case in cases:
        results.append(run(source, reference if case.get('reference') else binary, case, out))
        write_json(out / 'results.json', results)
    summarize_comparisons(results, out)


if __name__ == '__main__':
    main()
