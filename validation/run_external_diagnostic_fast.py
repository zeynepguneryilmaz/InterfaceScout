from pathlib import Path
import runpy
src=Path('validation/run_external_diagnostic.py').read_text()
src=src.replace('nperm=10000','nperm=1000')
tmp=Path('validation/_external_fast.py'); tmp.write_text(src)
runpy.run_path(str(tmp),run_name='__main__')
