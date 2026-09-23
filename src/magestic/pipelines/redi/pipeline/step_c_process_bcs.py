"""Step C: registered molecular extraction with exhaustive read accounting.

Keys join by sample and library/barcoder identity. No missing sample, ambiguous
primer, malformed count or unregistered barcode becomes a physical clone.
Inputs are processed one file at a time to avoid a concatenated raw-read table.
"""
from __future__ import annotations
from collections import Counter
from dataclasses import asdict
from pathlib import Path
import re
import pandas as pd
from ..core.config import REDIConfig
from ..core.identity import keyed_frame
from ..core.io import provenance,publish_tables
from ..readers.fastq import extract_bc1_from_seq,extract_redi_bc_vectorized,_integer,_anchor
from ..readers.keyfiles import load_redi_bc_keys,_read_key

BC_ONLY='detect all bc1 on the plate'


def _key_paths(config):
    paths=[Path(p) if Path(p).is_absolute() else config.paths.keyfiles/p for p in config.redi_bc_keys]
    if len({p.resolve() for p in paths})!=len(paths):raise ValueError('duplicate REDI barcode key paths')
    return paths


def _load_redi_bc_keys_from_config(config: REDIConfig) -> pd.DataFrame:
    parts=[]
    for path in _key_paths(config):
        columns=set(_read_key(path).columns)
        if 'REDI_bc_seq' in columns and 'REDI_bc1_seq' not in columns:
            parts.append(load_redi_bc_keys(v5_key_path=path))
        elif {'REDI_bc1_seq','REDI_library'}<=columns and 'REDI_bc_seq' not in columns:
            parts.append(load_redi_bc_keys(yt182_yt183_key_path=path))
        else:raise ValueError(f'unsupported or ambiguous REDI barcode key schema: {path}')
    if not parts:return pd.DataFrame()
    return keyed_frame(pd.concat(parts,ignore_index=True),['REDI_library','REDI_bc'],label='all REDI barcode keys')


def _file_sample(path):
    if path.name.endswith('_R1_R2_counts.tsv'):return path.name[:-len('_R1_R2_counts.tsv')]
    if path.name.startswith('R1_R2_counts_') and path.suffix=='.tsv':return path.name[len('R1_R2_counts_'):-4]
    return None


def _validate_key_prefixes(library, size, positions):
    """A fixed-length prefix cannot identify two different physical wells.

    Keep every registered key. Same-position aliases are compatible, and
    another library is a distinct namespace. Compare only unequal lengths
    that could be sliced/padded into this assay's configured key length.
    """
    shorter_lengths={len(barcode) for barcode in positions if len(barcode)<size}
    for barcode,position in positions.items():
        if len(barcode)==size:
            for length in shorter_lengths:
                other=positions.get(barcode[:length])
                if other is not None and other!=position:
                    raise ValueError(f'REDI library {library!r}: prefix-ambiguous physical keys at barcode length {size}')
        elif len(barcode)>size:
            other=positions.get(barcode[:size])
            if other is not None and other!=position:
                raise ValueError(f'REDI library {library!r}: prefix-ambiguous physical keys at barcode length {size}')


def run_step_c(config: REDIConfig,combined_key_path: Path) -> tuple[Path,Path]:
    combined=config.paths.combined_data
    out=combined/'combined_REDI_bc_MAGESTIC_bc_merged_df.tsv'
    only_out=combined/'MAGESTIC_bc_only_grouped_df.tsv'
    rejected_out=combined/'unassigned_REDI_read_counts.tsv'
    prov=provenance('c_process_bcs',out)
    prov.add_input(combined_key_path,role='sample_assay_registration')
    meta=_read_key(combined_key_path)
    reserved=set(meta)&{'R1_seq','R2_seq','bc1','counts'}
    if reserved:raise ValueError(f'sample metadata has reserved measurement columns: {sorted(reserved)}')
    for required in ('sample_number','outer_primers','inner_primers','purpose'):
        meta=keyed_frame(meta,[required],label='sample metadata',unique=False)
    configured_date=str(config.sequencing_date)
    if 'sequencing_date' in meta and not meta.sequencing_date.eq(configured_date).all():
        raise ValueError('declared sequencing_date conflicts with configured dated assay registration')
    meta['sequencing_date']=configured_date
    generated=meta.sequencing_date+'_'+meta.outer_primers+'_'+meta.sample_number
    if 'sample_name' in meta and not meta.sample_name.equals(generated):
        raise ValueError('declared sample_name conflicts with dated assay/sample registration')
    meta['sample_name']=generated
    meta=keyed_frame(meta,['sample_name'],label='combined sample registration')
    if type(config.allow_partial_samples) is not bool:
        raise ValueError('allow_partial_samples must be boolean')
    if not config.primer_rules or any(side not in ('R1','R2') for side in config.primer_rules.values()):
        raise ValueError('primer rules must explicitly select R1 or R2')
    metadata={};sides={}
    for row in meta.to_dict(orient='records'):
        name=row['sample_name']
        matches=[side for primer,side in config.primer_rules.items()
                 if primer and re.search(r'(?<![A-Za-z0-9])'+re.escape(primer)+r'(?![A-Za-z0-9])',row['inner_primers'])]
        if len(matches)!=1:raise ValueError(f'{name}: missing or ambiguous primer rule')
        sides[name]=matches[0]
        status=row.get('registration_status') or ('declared_empty' if str(row.get('REDI_plate_barcode','')).casefold()=='empty' else 'registered')
        if status not in ('registered','declared_empty'):raise ValueError(f'{name}: invalid registration status')
        declared_empty=str(row.get('REDI_plate_barcode','')).casefold()=='empty'
        if (status=='declared_empty')!=declared_empty:
            raise ValueError(f'{name}: contradictory empty-control/physical plate registration')
        row['registration_status']=status
        if status=='registered' and row['purpose']!=BC_ONLY and not row.get('MAT_alpha_array_1536_barcode'):
            raise ValueError(f'{name}: missing physical array registration')
        if row['purpose']!=BC_ONLY:
            row['REDI_library']=row.get('REDI_library') or config.library_by_purpose.get(row['purpose'])
            if not row['REDI_library']:raise ValueError(f'{name}: explicit REDI library registration required')
            if 'MAT_a_barcode_prefix' not in row or 'MAT_a_barcode_length' not in row:
                raise ValueError(f'{name}: missing REDI molecular extraction schema')
        metadata[name]=row
    for path in _key_paths(config):prov.add_input(path,role='barcoder_position_key')
    keys=_load_redi_bc_keys_from_config(config)
    overlaps=(set(meta)&set(keys))-{'REDI_library'}
    if overlaps:raise ValueError(f'barcoder/sample metadata overlap/collision: {sorted(overlaps)}')
    lookup={} if keys.empty else {(r['REDI_library'],r['REDI_bc']):r for r in keys.to_dict(orient='records')}
    if any(r['purpose']!=BC_ONLY for r in metadata.values()) and not lookup:
        raise ValueError('REDI assays require a nonempty barcoder position key')
    key_lengths={};key_positions={}
    for (library,barcode),key in lookup.items():
        key_lengths.setdefault(library,set()).add(len(barcode))
        key_positions.setdefault(library,{})[barcode]=tuple(key[column] for column in
            ('REDI_plate_color_384','REDI_row_384','REDI_col_384'))
    checked_schemas=set()
    for name,row in metadata.items():
        if row['registration_status']!='registered' or row['purpose']==BC_ONLY:continue
        library=row['REDI_library']
        if library not in key_lengths:
            raise ValueError(f'{name}: REDI library {library!r} is absent from the barcode key')
        size=_integer(row['MAT_a_barcode_length'],'enabled REDI barcode length',minimum=1)
        _anchor(row['MAT_a_barcode_prefix'],'enabled REDI barcode anchor')
        if size not in key_lengths[library]:
            raise ValueError(f'{name}: REDI barcode length {size} is absent from library {library!r} key')
        if (library,size) not in checked_schemas:
            _validate_key_prefixes(library,size,key_positions[library])
            checked_schemas.add((library,size))
    all_paths=sorted(config.paths.counts_tables.glob('*.tsv'))
    paths=[p for p in all_paths if _file_sample(p) is not None]
    if not paths:raise FileNotFoundError(f'no declared read-count tables in {config.paths.counts_tables}')
    redi=Counter();only=Counter();rejected=Counter();seen=set();input_reads=input_rows=0
    for path in paths:
        file_sample=_file_sample(path)
        if file_sample not in metadata:raise ValueError(f'unregistered input sample: {file_sample}')
        if file_sample in seen:raise ValueError(f'duplicate sample count files: {file_sample}')
        seen.add(file_sample);prov.add_input(path,role='collapsed_read_counts')
        frame=_read_key(path)
        if not {'sample_name','R1_seq','R2_seq','counts'}<=set(frame):raise ValueError(f'{path}: invalid count-table schema')
        if set(frame.sample_name)-{file_sample}:raise ValueError(f'{path}: sample identities disagree with count filename')
        counts=[_integer(value,'read counts',minimum=1) for value in frame.counts]
        input_rows+=len(frame);input_reads+=sum(counts)
        row=metadata[file_sample];side=sides[file_sample];other='R2' if side=='R1' else 'R1'
        barcodes=frame[f'{side}_seq'].map(lambda seq:extract_bc1_from_seq(seq,**{
            'fwd_bc1_prefix':config.barcode.fwd_bc1_prefix,'fwd_bc1_suffix':config.barcode.fwd_bc1_suffix,
            'bc1_length':config.barcode.bc1_length,'length_tolerance':config.barcode.bc1_length_tolerance,
            'fallback_prefix_length':config.barcode.fallback_prefix_length,
            'allow_positional_fallback':config.barcode.allow_positional_fallback}))
        redi_barcodes=[None]*len(frame) if row['purpose']==BC_ONLY else extract_redi_bc_vectorized(
            frame[f'{other}_seq'],pd.Series(row['MAT_a_barcode_prefix'],index=frame.index),
            pd.Series(row['MAT_a_barcode_length'],index=frame.index))
        for bc1,rb,count in zip(barcodes,redi_barcodes,counts):
            reason=None
            if row['registration_status']=='declared_empty':reason='declared_empty_control'
            elif bc1 is None:reason='uncallable_clone_barcode'
            elif row['purpose']==BC_ONLY:only[file_sample,bc1]+=count;continue
            elif rb is None:reason='uncallable_REDI_barcode'
            elif (row['REDI_library'],rb) not in lookup:reason='unregistered_REDI_barcode'
            else:redi[file_sample,rb,bc1]+=count;continue
            rejected[file_sample,reason,bc1 or '',rb or '']+=count
    missing=set(metadata)-seen
    if missing and not config.allow_partial_samples:raise ValueError(f'missing registered sample count files: {sorted(missing)}')
    records=[]
    for (sample,rb,bc1),count in sorted(redi.items()):
        row=metadata[sample];key=lookup[row['REDI_library'],rb]
        overlaps=(set(row)&set(key))-{'REDI_library'}
        if overlaps:raise ValueError(f'barcoder/sample metadata collision: {sorted(overlaps)}')
        records.append({**row,**key,'bc1':bc1,'counts':count})
    fields=list(dict.fromkeys([*meta.columns,'REDI_library','REDI_bc','REDI_plate_color_384','REDI_row_384','REDI_col_384','bc1','counts']))
    annotated=pd.DataFrame(records) if records else pd.DataFrame(columns=fields)
    only_frame=pd.DataFrame([{'sample_name':s,'bc1':b,'counts':n} for (s,b),n in sorted(only.items())],columns=['sample_name','bc1','counts'])
    rejected_frame=pd.DataFrame([{'sample_name':s,'reason':why,'bc1':b,'REDI_bc':rb,'counts':n} for (s,why,b,rb),n in sorted(rejected.items())],columns=['sample_name','reason','bc1','REDI_bc','counts'])
    assigned_redi=sum(redi.values());assigned_only=sum(only.values());unassigned=sum(rejected.values())
    if input_reads!=assigned_redi+assigned_only+unassigned:raise AssertionError('read accounting did not reconcile')
    prov.set_parameters(barcode=asdict(config.barcode),primer_rules=config.primer_rules,
        library_by_purpose=config.library_by_purpose,allow_partial_samples=config.allow_partial_samples,
        missing_samples=sorted(missing),ignored_non_count_files=[str(p) for p in all_paths if p not in paths])
    prov.set_counts(input_rows=input_rows,input_reads=input_reads,registered_samples=len(metadata),
                    input_samples=len(seen),assigned_redi_reads=assigned_redi,
                    assigned_bc1_only_reads=assigned_only,unassigned_reads=unassigned,
                    physical_output_rows=len(annotated))
    publish_tables([(out,annotated),(only_out,only_frame),(rejected_out,rejected_frame)],prov)
    return out,only_out
