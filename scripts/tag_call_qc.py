"""Summarize GEO's sample-tag calls, including barcodes absent from per-tag matrices."""
import csv,gzip,collections,json,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[1]
records=[]
for path in sorted((ROOT/'data/raw').glob('*Sample_Tag_Calls.csv.gz')):
    with gzip.open(path,'rt') as f:
        rows=list(csv.DictReader(line for line in f if line.strip() and not line.startswith('#')))
    records.append(dict(file=path.name,total_called_barcodes=len(rows),calls=dict(collections.Counter(r['Sample_Tag'] for r in rows))))
(ROOT/'results/sc/review/tag_call_summary.json').write_text(json.dumps(records,indent=2)+'\n')
