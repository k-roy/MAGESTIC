"""Atomic REDI artifact publication with source and input attestations."""
import fcntl
import hashlib
import os
import shutil
from pathlib import Path
import tempfile
from magestic.utils.provenance import StepProvenance, md5_file


# Modules OUTSIDE pipelines/redi that REDI code imports, relative to the magestic
# package. Their bytes change REDI behavior or receipts, so they belong to the REDI
# source identity. Only these files are hashed, not their import closure: following
# wgs/readers/__init__ would invalidate every REDI cache on each WGS-core edit.
# test_redi_source_identity_scope.py fails when REDI gains an unlisted cross-package import.
EXTERNAL_SOURCES = ('utils/provenance.py', 'pipelines/wgs/readers/wgs_calls.py')


def source_files():
    """Every source file in the REDI identity: the REDI package plus EXTERNAL_SOURCES.

    The Snakefile's COMMON input list uses this too, so scheduling and the digest
    cannot drift apart.
    """
    redi=Path(__file__).resolve().parents[1]
    package=redi.parents[1]
    return [*sorted(redi.rglob('*.py')), *(package/relative for relative in EXTERNAL_SOURCES)]


def source_digest():
    package=Path(__file__).resolve().parents[3]
    digest=hashlib.sha256()
    for path in source_files():
        digest.update(str(path.relative_to(package)).encode()+b'\0'+path.read_bytes()+b'\0')
    return digest.hexdigest()


def provenance(step,output):
    result=StepProvenance('redi',step,output)
    result.code['redi_source_sha256']=source_digest()
    return result


def publish_tables(tables,prov,*,files=()):
    """Stage every TSV before publication; write completion sidecar last.

    Cooperating writers serialize publication. Readers must use the sidecar
    as completion evidence; individual replacements are atomic, a file set is
    not a filesystem transaction. Input/source mutation prevents publication.
    """
    main=Path(prov.output);main.parent.mkdir(parents=True,exist_ok=True)
    temporary=[]
    try:
        for path,frame in tables:
            path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
            fd,name=tempfile.mkstemp(prefix=path.name+'.',suffix='.tmp',dir=path.parent)
            temporary.append((Path(name),path,len(frame)))
            with os.fdopen(fd,'w') as handle:frame.to_csv(handle,sep='\t',index=False)
        for staged,path,rows in files:
            path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
            fd,name=tempfile.mkstemp(prefix=path.name+'.',suffix='.tmp',dir=path.parent)
            os.close(fd)
            temporary.append((Path(name),path,rows))
            shutil.copyfile(staged,name)
        with Path(str(main)+'.lock').open('a') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX)
            if source_digest()!=prov.code['redi_source_sha256']:
                raise ValueError('REDI source changed while stage was running')
            for item in prov.inputs:
                if md5_file(item['path'])!=item['md5']:
                    raise ValueError(f"REDI input changed while stage was running: {item['path']}")
            Path(prov.sidecar_path()).unlink(missing_ok=True)
            for staged,path,rows in temporary:
                os.replace(staged,path)
                prov.add_output(path,**({'rows':rows} if rows is not None else {}))
            prov.write()
    finally:
        for staged,_,_ in temporary:staged.unlink(missing_ok=True)
