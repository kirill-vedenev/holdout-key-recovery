#!/usr/bin/env python3
"""Execute the toy notebooks with fresh SageMath kernels and save their outputs."""
import argparse
import hashlib
import json
import os
import platform
from pathlib import Path
from time import perf_counter

import nbformat
from nbclient import NotebookClient
from jupyter_client.kernelspec import KernelSpecManager

HERE = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('notebooks', nargs='*', help='Notebook names; default: all four')
    parser.add_argument('--kernel', help='Registered SageMath kernel name')
    parser.add_argument('--timeout', type=int, default=600, help='Per-cell timeout in seconds')
    args = parser.parse_args()
    kernels = KernelSpecManager().find_kernel_specs()
    available = sorted(name for name in kernels if name.lower().startswith('sagemath'))
    kernel = args.kernel or (available[-1] if available else None)
    if kernel not in kernels:
        parser.error('Register a SageMath Jupyter kernel or supply --kernel NAME')
    paths = [HERE/name for name in args.notebooks] if args.notebooks else sorted(HERE.glob('[0-9][0-9]-*.ipynb'))
    summaries = []
    for path in paths:
        notebook = nbformat.read(path, as_version=4)
        for cell in notebook.cells:
            cell.metadata.pop('execution', None)
            if cell.cell_type == 'code':
                cell.outputs = []
                cell.execution_count = None
        print(f'{path.name}: running with a fresh kernel', flush=True)
        started = perf_counter()
        try:
            NotebookClient(notebook, kernel_name=kernel, timeout=args.timeout,
                           resources={'metadata': {'path': str(HERE)}},
                           record_timing=False, allow_errors=False).execute()
        except Exception:
            failure = HERE/'runs'/'failed'
            failure.mkdir(parents=True, exist_ok=True)
            nbformat.write(notebook, failure/path.name)
            raise
        nbformat.validate(notebook)
        nbformat.write(notebook, path)
        result = json.loads((HERE/'results'/f'{path.stem}.json').read_text())
        if result.get('status') != 'passed':
            raise RuntimeError(f'{path.name}: missing success record')
        duration = perf_counter()-started
        summary = dict(notebook=path.name, status='passed', kernel=kernel,
                       executed_code_cells=sum(c.cell_type=='code' for c in notebook.cells),
                       execution_seconds=duration,
                       notebook_sha256=hashlib.sha256(path.read_bytes()).hexdigest())
        summaries.append(summary)
        print(f'{path.name}: passed in {duration:.2f} seconds', flush=True)
    report = dict(status='passed', notebooks=summaries,
                  reference_sha256=hashlib.sha256((HERE/'reference.py').read_bytes()).hexdigest(),
                  python_version=platform.python_version(), platform=platform.platform(),
                  logical_cpu_count=os.cpu_count(),
                  execution_scope='selected' if args.notebooks else 'all')
    name = 'selected-execution-summary.json' if args.notebooks else 'execution-summary.json'
    (HERE/'results'/name).write_text(json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    main()
