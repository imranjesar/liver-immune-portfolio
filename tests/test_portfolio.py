import csv,gzip,pathlib,tempfile,unittest,sys
ROOT=pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from convert_sc import convert
class PortfolioTests(unittest.TestCase):
 def test_sparse_orientation_and_counts(self):
  with tempfile.TemporaryDirectory() as d:
   p=pathlib.Path(d);src=p/'input.csv.gz'
   with gzip.open(src,'wt') as f:f.write('# BD header\nCell_Index,GeneA,mt-GeneB\nc1,0,2\nc2,3,0\n')
   self.assertEqual(convert(src,p/'out'),(2,2,2))
   with gzip.open(p/'out/matrix.mtx.gz','rt') as f:self.assertEqual(f.read().splitlines()[2:],['2 1 2','1 2 3'])
 def test_negative_rejected(self):
  with tempfile.TemporaryDirectory() as d:
   p=pathlib.Path(d);src=p/'bad.gz'
   with gzip.open(src,'wt') as f:f.write('Cell_Index,G\nc,-1\n')
   with self.assertRaises(ValueError):convert(src,p/'out')
 def test_manifest_and_replicates(self):
  with (ROOT/'config/files.tsv').open() as f: rows=list(csv.DictReader(f,delimiter='\t'))
  self.assertEqual(len(rows),13)
  self.assertTrue(all(r['url'].startswith('https://ftp.ncbi.nlm.nih.gov/geo/') for r in rows))
  with (ROOT/'config/sc_samples.tsv').open() as f: samples=list(csv.DictReader(f,delimiter='\t'))
  self.assertEqual(len(samples),9);self.assertEqual(len({s['sample_id'] for s in samples}),9)
  self.assertEqual(sorted([sum(s['condition']==c for s in samples) for c in {s['condition'] for s in samples}]),[3,3,3])
if __name__=='__main__':unittest.main()
