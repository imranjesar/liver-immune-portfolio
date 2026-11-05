"""Download only allowlisted public GEO files, recording SHA256 provenance."""
import argparse, csv, hashlib, json, pathlib, urllib.request
ROOT = pathlib.Path(__file__).resolve().parents[1]
def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--study', choices=['GSE338157','GSE247373','all'], default='all'); args=ap.parse_args()
    records=[]
    manifest=ROOT/"metadata"/f"downloads-{args.study}.json"
    previous={r["filename"]:r for r in json.loads(manifest.read_text())} if manifest.exists() else {}
    (ROOT/"data/raw").mkdir(parents=True,exist_ok=True)
    with (ROOT/'config/files.tsv').open() as f: rows=list(csv.DictReader(f,delimiter='\t'))
    for row in rows:
        if args.study != 'all' and row['accession'] != args.study: continue
        if not row['url'].startswith('https://ftp.ncbi.nlm.nih.gov/geo/'): raise ValueError('Non-GEO URL')
        dest=ROOT/'data/raw'/row['filename']; tmp=dest.with_suffix(dest.suffix+'.part')
        if not dest.exists():
            print('Downloading',dest.name,flush=True)
            with urllib.request.urlopen(row['url'],timeout=180) as src, tmp.open('wb') as out:
                while chunk:=src.read(1024*1024): out.write(chunk)
            tmp.replace(dest)
        h=hashlib.sha256()
        with dest.open('rb') as f:
            while chunk:=f.read(1024*1024): h.update(chunk)
        record=dict(row,bytes=dest.stat().st_size,sha256=h.hexdigest())
        if dest.name in previous and record['sha256'] != previous[dest.name]['sha256']:
            raise ValueError(f'Checksum changed: {dest.name}; investigate before updating provenance')
        records.append(record)
    (ROOT/'metadata'/f'downloads-{args.study}.json').write_text(json.dumps(records,indent=2)+'\n')
if __name__=='__main__': main()
