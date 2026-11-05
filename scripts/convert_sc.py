"""Stream BD cell-by-gene CSV into sparse Matrix Market; never densify all cells."""
import csv,gzip,pathlib,tempfile,shutil,argparse
ROOT=pathlib.Path(__file__).resolve().parents[1]
def convert(src,out):
    out.mkdir(parents=True,exist_ok=True)
    with gzip.open(src,'rt') as f:
        lines=(l for l in f if l.strip() and not l.startswith('#'))
        reader=csv.reader(lines); header=next(reader)
        if header[0]!='Cell_Index' or len(set(header[1:]))!=len(header)-1: raise ValueError('Unexpected gene header')
        genes=header[1:]; cells=[]; nnz=0
        with tempfile.TemporaryFile(mode='w+t') as body:
            for j,row in enumerate(reader,1):
                if len(row)!=len(header): raise ValueError('Ragged counts')
                cells.append(row[0])
                for i,v in enumerate(row[1:],1):
                    n=int(v)
                    if n<0: raise ValueError('Negative counts')
                    if n: body.write(f'{i} {j} {n}\n'); nnz+=1
            if len(set(cells))!=len(cells) or not cells: raise ValueError('Invalid cell IDs')
            body.seek(0)
            with gzip.open(out/'matrix.mtx.gz','wt') as dest:
                dest.write(f'%%MatrixMarket matrix coordinate integer general\n{len(genes)} {len(cells)} {nnz}\n');shutil.copyfileobj(body,dest)
    for name,values in [('features.tsv.gz',genes),('barcodes.tsv.gz',cells)]:
        with gzip.open(out/name,'wt') as f:f.write('\n'.join(values)+'\n')
    return len(genes),len(cells),nnz
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--pattern',default='*RSEC_MolsPerCell.csv.gz');args=ap.parse_args()
    for src in sorted((ROOT/'data/raw').glob(args.pattern)):
        print(src.name,convert(src,ROOT/'data/processed'/src.name.split('_mm_')[0]),flush=True)
