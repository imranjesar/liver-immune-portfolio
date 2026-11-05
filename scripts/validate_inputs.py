"""Validate downloaded provenance and tag-to-cell assignments before analysis."""
import csv,gzip,hashlib,json,pathlib,re
ROOT=pathlib.Path(__file__).resolve().parents[1]
def readrows(path):
    with gzip.open(path,'rt') as f:
        return list(csv.DictReader(l for l in f if l.strip() and not l.startswith('#')))
def main():
    records=json.loads((ROOT/'metadata/downloads-all.json').read_text())
    for r in records:
        p=ROOT/'data/raw'/r['filename'];h=hashlib.sha256()
        with p.open('rb') as f:
            while chunk:=f.read(1024*1024):h.update(chunk)
        assert p.stat().st_size==r['bytes'] and h.hexdigest()==r['sha256'],p.name
    with (ROOT/'config/sc_samples.tsv').open() as f:samples=list(csv.DictReader(f,delimiter='\t'))
    report=[];seen={}
    for s in samples:
        tag=re.search('SampleTag(\\d+)',s['matrix_dir'])[1]
        callpaths=list((ROOT/'data/raw').glob(s['gsm']+'*Sample_Tag_Calls.csv.gz'))
        assert len(callpaths)==1
        calls=readrows(callpaths[0]); assignments={r['Cell_Index']:r['Sample_Tag'] for r in calls}
        assert len(assignments)==len(calls),'Duplicate tag-call IDs'
        with gzip.open(ROOT/'data/processed'/s['matrix_dir']/'barcodes.tsv.gz','rt') as f:cells=f.read().splitlines()
        assert all(assignments.get(c)==f'SampleTag{tag}_mm' for c in cells),s
        overlap=set(cells)&seen.get(s['gsm'],set());assert not overlap,'Cell repeated across mouse tags'
        seen.setdefault(s['gsm'],set()).update(cells)
        report.append(dict(s,cells=len(cells),tag_calls_verified=True))
    (ROOT/'metadata/input-validation.json').write_text(json.dumps(report,indent=2)+'\n')
    print('Checksums and cell-level tag assignments verified for',len(report),'samples')
if __name__=='__main__':main()
