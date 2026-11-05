"""Audit dimensions, embedded PDF fonts, vector output and PNG resolution with system tools."""
import pathlib,subprocess,csv,re,struct,json,hashlib,argparse
ap=argparse.ArgumentParser();ap.add_argument("--check",action="store_true");args=ap.parse_args()
ROOT=pathlib.Path(__file__).resolve().parents[1]
pub=ROOT/'figures/publication'
records=[]
for row in csv.DictReader((pub/'export_specifications.csv').open()):
 name=row['figure'];w=float(row['width_mm']);h=float(row['height_mm']);dpi=int(row['png_dpi'])
 pdf=pub/(name+'.pdf');png=pub/(name+'.png');svg=pub/(name+'.svg')
 for f in [pdf,png,svg]:assert f.exists() and f.stat().st_size>1000,f
 info=subprocess.check_output(['pdfinfo',str(pdf)],text=True)
 assert re.search(r'Pages:\s+1\b',info),info
 m=re.search(r'Page size:\s+([\d.]+) x ([\d.]+) pts',info)
 assert m and abs(float(m[1])*25.4/72-w)<=.353 and abs(float(m[2])*25.4/72-h)<=.353
 fonts=subprocess.check_output(['pdffonts',str(pdf)],text=True).splitlines()[2:]
 assert fonts and all(re.search(r'\byes\s+yes\s+yes\b',f) for f in fonts),fonts
 images=subprocess.check_output(['pdfimages','-list',str(pdf)],text=True).splitlines()[2:]
 assert not [line for line in images if line.strip()], 'Raster image in vector PDF'
 b=png.read_bytes();assert b[:8]==b'\x89PNG\r\n\x1a\n'
 width,height=struct.unpack('>II',b[16:24]);assert abs(width-w/25.4*dpi)<2 and abs(height-h/25.4*dpi)<2
 pos=8;resolution=None
 while pos<len(b):
  n=struct.unpack('>I',b[pos:pos+4])[0];kind=b[pos+4:pos+8];data=b[pos+8:pos+8+n]
  if kind==b'pHYs':
   px,py,unit=struct.unpack('>IIB',data);resolution=(px*.0254,py*.0254,unit)
  pos+=12+n
 assert resolution and abs(resolution[0]-dpi)<1 and abs(resolution[1]-dpi)<1 and resolution[2]==1
 assert '<svg' in svg.read_text()[:1500]
 records.append(dict(figure=name,width_mm=w,height_mm=h,png_pixels=[width,height],png_dpi=dpi,embedded_fonts=True,vector_pdf=True,sha256={f.suffix[1:]:hashlib.sha256(f.read_bytes()).hexdigest()for f in [pdf,png,svg]}))
manifest=ROOT/'metadata/publication_figure_manifest.json'
if args.check:
 assert json.loads(manifest.read_text())==records, 'Figure manifest is stale; regenerate and review exports'
else:
 manifest.write_text(json.dumps(records,indent=2)+'\n')
print('Verified',len(records),'figure sets: dimensions, 600 dpi, embedded fonts and vector PDF/SVG')
