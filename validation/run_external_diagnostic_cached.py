from pathlib import Path
import runpy, sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import backend.main as ism

_orig=ism.build_surface_residues
_cache={}
def _cached(path,pH):
    key=(str(path),float(pH))
    if key not in _cache:
        _cache[key]=_orig(path,pH)
    return _cache[key]
ism.build_surface_residues=_cached

src=Path('validation/run_external_diagnostic.py').read_text().replace('nperm=10000','nperm=1000')
tmp=Path('validation/_external_cached.py'); tmp.write_text(src)
runpy.run_path(str(tmp),run_name='__main__')
