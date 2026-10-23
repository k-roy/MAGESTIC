"""Validate paired FASTQ, trim/merge with BBTools, then collapse on disk.

The counts table carries a merged consensus and its reverse complement. These
columns are not the original sequencing mates. This representation is recorded
explicitly so downstream anchor searches cannot mistake a reconstruction for
independent read evidence. Unmerged pairs are counted, but are not assigned.
"""
from __future__ import annotations

import csv
import fcntl
import gzip
import hashlib
from itertools import zip_longest
import json
from pathlib import Path
import re
import shutil
import sqlite3
import subprocess
import tempfile

from ..core.config import REDIConfig
from ..core.io import provenance, publish_tables
from magestic.utils.provenance import md5_file

TRUSEQ_ADAPTER = "AGATCGGAAGAGCGTCGTGTAGGGAAAGA"
_COMPLEMENT = str.maketrans('ACGTRYKMSWBDHVN', 'TGCAYRMKSWVHDBN')


def _integer(value, name, minimum):
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ValueError(f'{name} must be an integer >= {minimum}')
    return value


def _fastq(path, *, expected_mate=None):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', encoding='ascii') as handle:
        record = 0
        while True:
            header = handle.readline()
            if not header:
                return
            record += 1
            lines = [handle.readline() for _ in range(3)]
            if any(not line for line in lines):
                raise ValueError(f'Truncated FASTQ record {record}: {path}')
            sequence, plus, quality = [line.rstrip('\r\n') for line in lines]
            if (not header.startswith('@') or not header[1:].strip() or
                    not plus.startswith('+') or not sequence or
                    len(sequence) != len(quality) or
                    not re.fullmatch('[ACGTRYKMSWBDHVNacgtrykmswbdhvn]+', sequence) or
                    any(not 33 <= ord(char) <= 126 for char in quality)):
                raise ValueError(f'Invalid FASTQ record {record}: {path}')
            header_fields = header[1:].split()
            name = header_fields[0]
            declared_mates = []
            if name.endswith(('/1', '/2')):
                declared_mates.append(int(name[-1]))
                name = name[:-2]
            # Both old /1,/2 and CASAVA read-number annotations carry mate
            # identity. Unlabelled formats remain supported, but explicit
            # declarations cannot contradict one another or the input side.
            if len(header_fields) > 1:
                casava = re.match(r'^([12]):[YN]:\d+:', header_fields[1])
                if casava:
                    declared_mates.append(int(casava.group(1)))
            if len(set(declared_mates)) > 1 or (declared_mates and expected_mate is not None
                                               and declared_mates[0] != expected_mate):
                raise ValueError(f'FASTQ mate identity/read number mismatch in record {record}: {path}')
            yield name, sequence.upper(), quality


def _validate_pairs(first, second):
    count = 0
    for left, right in zip_longest(_fastq(first, expected_mate=1), _fastq(second, expected_mate=2)):
        if left is None or right is None:
            raise ValueError('Unequal FASTQ mate pair counts')
        if left[0] != right[0]:
            raise ValueError(f'FASTQ mate identity mismatch at pair {count + 1}')
        count += 1
    return count


def _tool_identity():
    """Fingerprint actual Java classes, not just the shell launcher/version."""
    identity = {'java': shutil.which('java'), 'tools': {}}
    if identity['java']:
        identity['java'] = str(Path(identity['java']).resolve())
        identity['java_md5'] = md5_file(identity['java'])
    installations = {}
    for name in ('bbduk.sh', 'bbmerge.sh'):
        launcher = shutil.which(name)
        if not launcher:
            identity['tools'][name] = None
            continue
        launcher = Path(launcher).resolve()
        root = launcher.parent
        if root not in installations:
            digest = hashlib.sha256(); count = size = 0
            for folder in ('current', 'jni'):
                for path in sorted((root/folder).rglob('*')):
                    if path.is_file():
                        data = path.read_bytes()
                        digest.update(str(path.relative_to(root)).encode()+b'\0'+hashlib.sha256(data).digest())
                        count += 1; size += len(data)
            installations[root] = dict(implementation_sha256=digest.hexdigest(), files=count, bytes=size)
        identity['tools'][name] = dict(launcher=str(launcher), launcher_md5=md5_file(launcher), **installations[root])
    return identity


def _run(cmd):
    # BBTools shell wrappers interpolate $@ through eval, which splits paths
    # containing spaces. Invoke the installed Java class directly, with an argv.
    launcher, java = shutil.which(cmd[0]), shutil.which('java')
    if not launcher or not java:
        raise FileNotFoundError(f'Required BBTools launcher/Java absent from PATH: {cmd[0]}')
    root = Path(launcher).resolve().parent
    classname = {'bbduk.sh': 'jgi.BBDuk', 'bbmerge.sh': 'jgi.BBMerge'}[cmd[0]]
    if not (root/'current'/Path(*classname.split('.')).with_suffix('.class')).is_file():
        raise FileNotFoundError(f'Unsupported BBTools installation: missing class {classname} under {root}/current')
    heap = [arg for arg in cmd[1:] if arg.startswith('-Xmx')]
    arguments = [arg for arg in cmd[1:] if not arg.startswith('-Xmx')]
    if len(heap) != 1:
        raise ValueError('Exactly one explicit BBTools heap limit is required')
    cmd = [str(Path(java).resolve()), '-ea', heap[0], heap[0].replace('-Xmx', '-Xms'),
           f'-Djava.library.path={root / "jni"}', '-cp', str(root/'current'), classname, *arguments]
    result = subprocess.run(cmd, check=True, capture_output=True, text=True)
    return {'argv': cmd, 'executable_md5': md5_file(cmd[0]),
            'stdout': result.stdout, 'stderr': result.stderr}


def _bbduk(in1, in2, out1, out2, *, ftl=20, overwrite=True, threads=1, memory_mb=512):
    """Legacy adapter/end trimming; no quality threshold is silently added."""
    return _run(['bbduk.sh', f'-Xmx{memory_mb}m', f'threads={threads}',
                 f'in1={in1}', f'in2={in2}', f'out1={out1}', f'out2={out2}',
                 f'literal={TRUSEQ_ADAPTER}', f'overwrite={str(overwrite).lower()}',
                 f'ftl={ftl}', 'ktrim=r', 'k=23', 'mink=11', 'hdist=1', 'tpe', 'tbo'])


def _bbmerge(in1, in2, out, *, overwrite=True, threads=1, memory_mb=512):
    return _run(['bbmerge.sh', f'-Xmx{memory_mb}m', f'threads={threads}',
                 f'in1={in1}', f'in2={in2}', f'out={out}',
                 f'overwrite={str(overwrite).lower()}'])


def _collapse(merged, collapsed, counts, sample, min_counts, database):
    """Aggregate full consensus sequences in SQLite with a 16 MiB page cache."""
    merged_count = retained = unique = below = 0
    with sqlite3.connect(database) as db:
        db.execute('PRAGMA cache_size=-16384')
        db.execute('PRAGMA temp_store=FILE')
        db.execute('CREATE TABLE sequences(seq TEXT PRIMARY KEY, quality TEXT, n INTEGER NOT NULL)')
        for _, sequence, quality in _fastq(merged):
            db.execute('INSERT INTO sequences VALUES (?,?,1) ON CONFLICT(seq) DO UPDATE SET n=n+1',
                       (sequence, quality))
            merged_count += 1
        db.commit()
        with open(collapsed, 'w') as fastq, open(counts, 'w', newline='') as table:
            writer = csv.writer(table, delimiter='\t', lineterminator='\n')
            writer.writerow(['R1_seq', 'R2_seq', 'counts', 'sample_name'])
            for sequence, quality, count in db.execute('SELECT seq,quality,n FROM sequences ORDER BY n DESC,seq'):
                if count < min_counts:
                    below += count
                    continue
                retained += count
                unique += 1
                fastq.write(f'@{count}\n{sequence}\n+\n{quality}\n')
                writer.writerow([sequence, sequence.translate(_COMPLEMENT)[::-1], count, sample])
    assert merged_count == retained + below
    return dict(merged_reads=merged_count, retained_merged_reads=retained,
                retained_unique_sequences=unique, below_minimum_count_reads=below)


def _cache_matches(prov, output, expected_outputs):
    try:
        saved = json.loads(Path(prov.sidecar_path()).read_text())
        if (saved.get('pipeline') != 'redi' or saved.get('step') != 'a_trim_counts' or
                saved.get('parameters') != prov.parameters or
                saved.get('code', {}).get('redi_source_sha256') != prov.code['redi_source_sha256']):
            return False
        expected = {(item['path'], item['md5']) for item in prov.inputs}
        if expected != {(item['path'], item['md5']) for item in saved['inputs']}:
            return False
        outputs = saved['outputs']
        expected_paths = {Path(path).resolve() for path in expected_outputs}
        recorded_paths = [Path(item['path']).resolve() for item in outputs]
        return (len(expected_paths) == 5 and len(recorded_paths) == 5 and
                set(recorded_paths) == expected_paths and Path(output).resolve() in expected_paths and
                all(md5_file(item['path']) == item['md5'] for item in outputs))
    except (OSError, ValueError, KeyError, TypeError):
        return False


def run_step_a_one_sample(prefix: str, config: REDIConfig, *,
                          general_scripts_dir: str | None = None,
                          min_counts: int = 1, bbduk_ftl: int = 20,
                          overwrite: bool = True, threads: int = 1,
                          memory_mb: int = 512) -> None:
    if (not isinstance(prefix, str) or not prefix.strip() or prefix in ('.', '..') or
            any(c in prefix for c in '/\\\x00\t\r\n')):
        raise ValueError('sample prefix must be a nonempty safe file basename')
    if general_scripts_dir is not None:
        raise ValueError('general_scripts_dir is retired; Stage A uses native, attested collapse/counts')
    _integer(min_counts, 'min_counts', 1)
    _integer(bbduk_ftl, 'bbduk_ftl', 0)
    _integer(threads, 'threads', 1)
    _integer(memory_mb, 'memory_mb', 128)
    if not isinstance(overwrite, bool):
        raise ValueError('overwrite must be boolean')
    paths = config.paths
    final = paths.counts_tables / f'{prefix}_R1_R2_counts.tsv'
    expected_outputs = [paths.trimmed_fastq/f'{prefix}_trimmed_R1.fastq',
                        paths.trimmed_fastq/f'{prefix}_trimmed_R2.fastq',
                        paths.merged_fastq/f'{prefix}_trimmed_merged.fastq',
                        paths.collapsed_fastq/f'{prefix}_collapsed.fastq', final]
    final.parent.mkdir(parents=True, exist_ok=True)
    paths.intermediate.mkdir(parents=True, exist_ok=True)
    # The run lock covers cache inspection, computation and publication; the
    # distinct publication lock also protects readers/cooperating stage writers.
    with Path(str(final)+'.run.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        prov = provenance('a_trim_counts', final)
        raw1 = paths.raw_fastq / f'{prefix}_R1.fastq.gz'
        raw2 = paths.raw_fastq / f'{prefix}_R2.fastq.gz'
        prov.add_input(raw1, role='raw_R1'); prov.add_input(raw2, role='raw_R2')
        for name in ('bbduk.sh', 'bbmerge.sh'):
            if shutil.which(name):
                prov.add_input(Path(shutil.which(name)).resolve(), role='tool_launcher')
        prov.set_parameters(sample_name=prefix, min_counts=min_counts, bbduk_ftl=bbduk_ftl,
                            threads=threads, memory_mb=memory_mb, adapter=TRUSEQ_ADAPTER,
                            tool_implementation=_tool_identity(),
                            sequence_representation='merged_consensus_and_reverse_complement')
        if not overwrite and final.exists():
            if _cache_matches(prov, final, expected_outputs):
                return
            raise FileExistsError(f'Existing counts output lacks matching complete provenance: {final}')
        raw_count = _validate_pairs(raw1, raw2)
        with tempfile.TemporaryDirectory(prefix='redi_stage_a_', dir=paths.intermediate) as scratch:
            scratch = Path(scratch)
            trim1, trim2, merged, collapsed, table = [scratch / name for name in
                ('trim_R1.fastq', 'trim_R2.fastq', 'merged.fastq', 'collapsed.fastq', 'counts.tsv')]
            if raw_count:
                trim_record = _bbduk(raw1, raw2, trim1, trim2, ftl=bbduk_ftl,
                                      overwrite=True, threads=threads, memory_mb=memory_mb)
                trimmed_count = _validate_pairs(trim1, trim2)
                if trimmed_count > raw_count:
                    raise ValueError('Trimmed pair count exceeds raw input count')
                merge_record = _bbmerge(trim1, trim2, merged, overwrite=True,
                                        threads=threads, memory_mb=memory_mb) if trimmed_count else None
                if not trimmed_count:
                    merged.touch()
            else:
                trimmed_count = 0; trim_record = merge_record = None
                for path in (trim1, trim2, merged): path.touch()
            counts = _collapse(merged, collapsed, table, prefix, min_counts, scratch/'counts.sqlite')
            if counts['merged_reads'] > trimmed_count:
                raise ValueError('Merged read count exceeds trimmed input pair count')
            prov.set_counts(raw_read_pairs=raw_count, trimmed_read_pairs=trimmed_count,
                            trimmed_away_pairs=raw_count-trimmed_count,
                            unmerged_pairs=trimmed_count-counts['merged_reads'], **counts)
            prov.note(json.dumps({'bbduk': trim_record, 'bbmerge': merge_record}, sort_keys=True))
            if _tool_identity() != prov.parameters['tool_implementation']:
                raise ValueError('BBTools/Java implementation changed while processing')
            publish_tables([], prov, files=[
                (trim1, paths.trimmed_fastq/f'{prefix}_trimmed_R1.fastq', trimmed_count),
                (trim2, paths.trimmed_fastq/f'{prefix}_trimmed_R2.fastq', trimmed_count),
                (merged, paths.merged_fastq/f'{prefix}_trimmed_merged.fastq', counts['merged_reads']),
                (collapsed, paths.collapsed_fastq/f'{prefix}_collapsed.fastq', counts['retained_unique_sequences']),
                (table, final, counts['retained_unique_sequences'])])
