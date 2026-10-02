#!/usr/bin/env python3
"""Verify and unpack saved TII coefficients and any archived alternate encoding.

This performs file-format conversion only. It never computes a holdout kernel,
a jet, or a recovered key. --panel needs NumPy; unpacking uses the standard library.
"""
import argparse
from array import array
import hashlib
import json
from math import comb
from pathlib import Path
import shutil
import struct
import sys


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(8 << 20), b''):
            h.update(block)
    return h.hexdigest()


def check_file(path, record):
    if path.stat().st_size != record['bytes'] or digest(path) != record['sha256']:
        raise ValueError(f'Size or SHA-256 mismatch: {path}')


def restore_panel(kernel, source, output, metadata):
    import numpy as np

    k, d, n, count = (metadata[x] for x in ('variables', 'degree', 'monomials', 'polynomials'))
    words = metadata['panel_words']
    if metadata['format'] != 'THK1' or n != comb(k, d) or count > 64 * words:
        raise ValueError('Unsupported panel dimensions or ordering')
    polys = output / 'polynomials'
    polys.mkdir()
    packed = np.memmap(kernel, dtype=np.uint8, mode='r', offset=20,
                       shape=(count, (n + 7) // 8))
    panel_path = polys / 'kernel-panel.u64le'
    panel = np.memmap(panel_path, dtype=np.uint8, mode='w+', shape=(n, words * 8))
    for start in range(0, n, 65536):
        end = min(start + 65536, n)
        rows = np.unpackbits(np.ascontiguousarray(packed[:, start // 8:(end + 7) // 8]),
                             axis=1, bitorder='little')[:, :end - start].T
        encoded = np.packbits(rows, axis=1, bitorder='little')
        panel[start:end] = 0
        panel[start:end, :encoded.shape[1]] = encoded
    panel.flush()
    del panel, packed
    # Increasing fixed-weight masks are precisely colex subset order.
    with (polys / 'monomials-colex.u64').open('wb') as stream:
        mask = (1 << d) - 1
        for start in range(0, n, 65536):
            block = array('Q')
            for _ in range(min(65536, n - start)):
                block.append(mask)
                low = mask & -mask
                above = mask + low
                mask = above | (((above ^ mask) >> 2) // low)
            if sys.byteorder != 'little':
                block.byteswap()
            stream.write(block.tobytes())
    for name, record in metadata['derived_files'].items():
        check_file(polys / name, record)
    shutil.copyfile(source / 'polynomial-index.json', polys / 'index.json')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('instance', type=Path, help='An instance directory containing storage.json')
    ap.add_argument('--output', type=Path, help='Fresh destination; omit to verify chunks only')
    ap.add_argument('--panel', action='store_true', help='Also restore the exact THK1 panel and monomial map')
    args = ap.parse_args()
    source = args.instance.resolve()
    meta = json.loads((source / 'storage.json').read_text())
    if args.panel and (not args.output or meta['format'] != 'THK1'):
        ap.error('--panel requires --output and a THK1 instance')
    combined = hashlib.sha256()
    offset = 0
    for part in meta['chunks']:
        path = source / part['file']
        if path.parent != source or part['offset'] != offset:
            raise ValueError('Invalid chunk path or offset')
        check_file(path, part)
        with path.open('rb') as stream:
            for block in iter(lambda: stream.read(8 << 20), b''):
                combined.update(block)
        offset += part['bytes']
    if offset != meta['bytes'] or combined.hexdigest() != meta['sha256']:
        raise ValueError('Combined kernel size or SHA-256 mismatch')
    if digest(source / 'public.json') != meta['public_sha256']:
        raise ValueError('Public coordinate basis SHA-256 mismatch')
    alternates = meta.get('alternate_encodings', [])
    for alternate in alternates:
        path = source / alternate['file']
        if path.parent != source:
            raise ValueError('Invalid alternate-encoding path')
        check_file(path, alternate)
        with path.open('rb') as stream:
            header = struct.unpack('<4s4I', stream.read(20))
        expected = (alternate['format'].encode(), meta['variables'], meta['degree'],
                    meta['monomials'], meta['polynomials'])
        if header != expected:
            raise ValueError('Alternate kernel header differs from storage.json')
    if args.output:
        output = args.output.resolve()
        output.mkdir(parents=True, exist_ok=False)
        target = output / ('kernel.bin' if meta['format'] == 'WHK1' else 'kernel.thk1')
        with target.open('xb') as stream:
            for part in meta['chunks']:
                with (source / part['file']).open('rb') as chunk:
                    shutil.copyfileobj(chunk, stream, 8 << 20)
        check_file(target, meta)
        with target.open('rb') as stream:
            header = struct.unpack('<4s4I', stream.read(20))
        expected = (meta['format'].encode(), meta['variables'], meta['degree'],
                    meta['monomials'], meta['polynomials'])
        if header != expected:
            raise ValueError('Kernel header differs from storage.json')
        shutil.copyfile(source / 'public.json', output / 'public.json')
        for alternate in alternates:
            target_alternate = output / alternate['file']
            if target_alternate.exists():
                raise ValueError('Alternate encoding collides with another output')
            shutil.copyfile(source / alternate['file'], target_alternate)
        if args.panel:
            restore_panel(target, source, output, meta)
    print(f"{meta['instance']}: saved coefficients verified ({meta['polynomials']} polynomials)")


if __name__ == '__main__':
    main()
