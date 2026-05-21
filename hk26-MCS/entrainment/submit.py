"""
Smart SLURM submitter for calc_entrainment.py.

Checks which chunk done files are absent for the requested model/region, then
submits only those chunks as a SLURM array job.  If the zarr store has not been
initialised yet (no init.done sentinel), --init is run inline first.

Generated files:
    slurm/tasks/<timestamp>_<model>_<region>.json   — task list for the array
    slurm/scripts/<timestamp>_<model>_<region>.sh   — SLURM batch script
    slurm/output/<job_id>_<array_idx>.{out,err}     — SLURM logs

Usage:
    python submit.py --model um_glm_n2560_RAL3p3_tuned_hk26
    python submit.py --model um_glm_n2560_RAL3p3_tuned_hk26 --dry-run
"""
import argparse
import json
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import models

CHUNK_SIZE = 10      # must match CHUNK_SIZE in calc_entrainment.py

SBATCH_OPTS = {
    'account':       'hrcm',
    'ntasks':        '1',
    'cpus-per-task': '10',
    'mem':           '8G',
    'time':          '10:0:00',
    'partition':     'standard',
    'qos':           'high',
}


def count_zarr_times(zarr_path):
    import xarray as xr
    return xr.open_zarr(zarr_path).sizes['time']


def pending_chunks(model, region):
    zarr_path = models.data_dir(model) / f'entrainment_{region}.zarr'
    n_times   = count_zarr_times(zarr_path)
    n_chunks  = (n_times + CHUNK_SIZE - 1) // CHUNK_SIZE
    return [i for i in range(n_chunks)
            if not models.chunk_donefile(model, i).exists()]


def write_task_json(model, region, chunks, name):
    path = Path('slurm') / 'tasks' / f'{name}.json'
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(
        {'model': model, 'region': region, 'tasks': [{'chunk': c} for c in chunks]},
        indent=2,
    ))
    return path


def write_slurm_script(name, json_path, n_tasks):
    script_path = Path('slurm') / 'scripts' / f'{name}.sh'
    script_path.parent.mkdir(parents=True, exist_ok=True)
    output_dir  = Path('slurm') / 'output'
    output_dir.mkdir(parents=True, exist_ok=True)

    lines = ['#!/bin/bash']
    for k, v in SBATCH_OPTS.items():
        lines.append(f'#SBATCH --{k}={v}')
    lines += [
        f'#SBATCH --job-name=entr_{name[:20]}',
        f'#SBATCH --array=0-{n_tasks - 1}',
        f'#SBATCH --output={output_dir}/%A_%a.out',
        f'#SBATCH --error={output_dir}/%A_%a.err',
        '',
        f'python calc_entrainment.py {json_path} $SLURM_ARRAY_TASK_ID',
    ]
    script_path.write_text('\n'.join(lines) + '\n')
    return script_path


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    models.add_model_arg(parser)
    models.add_region_arg(parser)
    parser.add_argument('--dry-run', action='store_true',
                        help='Print what would happen without running sbatch or --init')
    args = parser.parse_args()

    model, region = args.model, args.region

    # --- Ensure zarr is initialised ---
    init_done = models.init_donefile(model, region)
    if not init_done.exists():
        print(f'init.done not found — running --init for {model} / {region}...')
        cmd = [sys.executable, 'calc_entrainment.py', '--init',
               '--model', model, '--region', region]
        if args.dry_run:
            print('  [dry-run]', ' '.join(cmd))
        else:
            subprocess.run(cmd, check=True)
    else:
        print(f'init.done exists — skipping --init.')

    if args.dry_run and not init_done.exists():
        print('[dry-run] Cannot compute pending chunks without zarr; stopping.')
        return

    # --- Find pending chunks ---
    chunks = pending_chunks(model, region)
    if not chunks:
        print('All chunks complete — nothing to submit.')
        return

    n = len(chunks)
    preview = str(chunks[:5]) + ('...' if n > 5 else '')
    print(f'{n} pending chunks: {preview}')

    # --- Write artefacts and submit ---
    ts   = datetime.now().strftime('%Y%m%d_%H%M%S')
    name = f'{ts}_{model}_{region}'

    json_path   = write_task_json(model, region, chunks, name)
    script_path = write_slurm_script(name, json_path, n)
    print(f'Task JSON:    {json_path}')
    print(f'SLURM script: {script_path}')

    if args.dry_run:
        print(f'[dry-run] Would run: sbatch {script_path}')
        return

    result = subprocess.run(['sbatch', str(script_path)],
                            capture_output=True, text=True)
    print(result.stdout.strip())
    if result.returncode != 0:
        print(result.stderr, file=sys.stderr)
        sys.exit(result.returncode)


if __name__ == '__main__':
    main()
