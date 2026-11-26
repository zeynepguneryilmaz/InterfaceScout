"""Consolidate publication research files without recomputing or losing evidence."""
from pathlib import Path, PurePosixPath
import csv, hashlib, json, os, posixpath, re, subprocess, zipfile
ROOT=Path.cwd()
def git(*args): return subprocess.check_output(['git',*args],text=True).strip()
old_head=git('rev-parse','HEAD')
paths=git('ls-files').splitlines()
M='research/method_verification'
P='research/protein_surface_analysis'
INDEX='research/README.md'
mapping={}
def panel(name):
    n=name.lower()
    return P if re.search(r'isv|protein.sol|external.method|freesasa|grid.convergence.8c50|alpha.lactalbumin|transferrin|curated|expanded.kinase|rnase|npbe|pep.patch|surface.assembly',n) else M
for p in paths:
    q=PurePosixPath(p); name=q.name
    if p.startswith('datasets/TV48/'):
        mapping[p]=M+('/protein_structures.csv' if name=='protein_manifest.csv' else '/README.md' if name=='README.md' else '/results/'+name)
    elif p.startswith('datasets/ISV76/'):
        sub='results' if name=='ISV76_results.zip' else 'method_comparisons' if re.search(r'FreeSASA|PEP_Patch|Protein_sol|external_aggregate',name) else 'electrostatics'
        mapping[p]=P+('/protein_structures.csv' if name=='protein_manifest.csv' else '/README.md' if name=='README.md' else '/'+sub+'/'+name)
    elif p.startswith('datasets/case_studies/'):
        mapping[p]=M+'/historical_literature/'+name
    elif p=='datasets/MANIFEST.json': mapping[p]='research/MANIFEST.json'
    elif p in ('datasets/README.md','validation/README.md','validation/PUBLICATION_EVIDENCE.md','validation/ACTIVE_DATASETS.md'): 
        mapping[p] = INDEX if p.endswith('README.md') else 'research/'+name
    elif p=='validation/strict_apo_48_selected_candidates_v2_final.csv': mapping[p]=M+'/protein_structures.csv'
    elif p=='validation/ISV_active_candidates_76_v3.csv': mapping[p]=P+'/protein_structures.csv'
    elif p.startswith('validation/'):
        group=panel(name)
        if p.startswith('validation/evidence/'): dest=group+'/source_records/evidence/'+name
        elif name in ('verify_publication.py','benchmark_case_results.csv','benchmark_manifest.json','parameter_selection.json'): dest=M+'/historical_literature/'+name
        elif name=='repository_consistency_check.py': dest='research/verify_datasets.py'
        elif name.endswith('.py'): dest=group+'/source_scripts/'+name
        else: dest=group+'/source_records/'+name
        mapping[p]=dest
    elif p.startswith('paper_validation/'): mapping[p]=P+'/source_scripts/'+name
    elif p.startswith('examples/'):
        mapping[p]=P+'/examples/'+p[len('examples/'):]
    elif p.startswith('.github/workflows/') and name not in ('publication-ci.yml','platform-release.yml','platform-release-packages.yml','consolidate-publication-release.yml','electrostatics.yml','npbe-production-regression.yml','real-backbone-recovery.yml','capture-real-readme-ui.yml','readme-svg-layout.yml','cancel-superseded-pr-checks.yml','cancel-obsolete-platform-release.yml','research-reorganization.yml'):
        body=Path(p).read_text()
        if re.search(r'validation/|datasets/|paper_validation/|examples/',body):
            mapping[p]=panel(name)+'/source_workflows/'+name
# Preserve all original source records, but remove exact duplicate manifests and replace redundant indexes.
assert Path('datasets/TV48/protein_manifest.csv').read_bytes()==Path('validation/strict_apo_48_selected_candidates_v2_final.csv').read_bytes()
assert Path('datasets/ISV76/protein_manifest.csv').read_bytes()==Path('validation/ISV_active_candidates_76_v3.csv').read_bytes()
original_manifest=json.loads(Path('datasets/MANIFEST.json').read_text())
old_text={p:Path(p).read_text(encoding='utf-8') for p in paths if p.endswith(('.md','.py','.yml','.yaml','.cff','.txt','.json')) and not p.startswith(('frontend/vendor/','licenses/'))}
for old,new in mapping.items():
    if old in ('validation/README.md','datasets/README.md','datasets/TV48/README.md','datasets/ISV76/README.md','validation/repository_consistency_check.py','validation/strict_apo_48_selected_candidates_v2_final.csv','validation/ISV_active_candidates_76_v3.csv'):
        Path(old).unlink();continue
    target=Path(new);target.parent.mkdir(parents=True,exist_ok=True)
    assert not target.exists(),(old,new)
    Path(old).rename(target)
# Rewrite markdown relative links using original location, preserving anchors and external URLs.
def doc_links(text,old,new):
    def fix(m):
        dest=m.group(1)
        if dest.startswith(('#','http:','https:','mailto:')): return m.group(0)
        pathname,sep,anchor=dest.partition('#')
        absolute=posixpath.normpath(posixpath.join(posixpath.dirname(old),pathname))
        mapped=mapping.get(absolute,absolute)
        if absolute.endswith('/') or (not Path(mapped).exists() and any(x.startswith(absolute+'/') for x in mapping)):
            descendants={posixpath.dirname(v) for k,v in mapping.items() if k.startswith(absolute.rstrip('/')+'/')}
            if absolute in ('validation','datasets'): mapped=INDEX
            elif absolute=='examples': mapped=P+'/examples/README.md'
        if not Path(mapped).exists() and mapped!=INDEX and not mapped.endswith('/README.md'):
            return m.group(0)
        rel=posixpath.relpath(mapped,posixpath.dirname(new) or '.')
        return ']('+rel+(sep+anchor if sep else '')+')'
    text=re.sub(r'\]\(([^)\s]+)\)',fix,text)
    return text
# New active paths in live scripts; historical source scripts retain exact bytes and use frozen checkout.
def root_paths(s):
    for old,new in sorted(mapping.items(),key=lambda x:-len(x[0])):
        s=s.replace(old,new)
    return s
for old,s in old_text.items():
    new=mapping.get(old,old)
    if not Path(new).is_file() or old in ('validation/README.md','datasets/README.md','datasets/TV48/README.md','datasets/ISV76/README.md','validation/repository_consistency_check.py'): continue
    if old.endswith('.md'):
        s=doc_links(s,old,new)
        if not new.startswith(('research/method_verification/source_records/','research/protein_surface_analysis/source_records/')):
            s=root_paths(s)
        Path(new).write_text(s,encoding='utf-8')
    elif new.startswith(('research/method_verification/source_scripts/','research/protein_surface_analysis/source_scripts/','research/method_verification/source_workflows/','research/protein_surface_analysis/source_workflows/')): continue
    elif old.startswith(('tests/','tools/','scripts/','.github/workflows/')) and old.endswith(('.py','.yml','.yaml')):
        if old=='.github/workflows/research-reorganization.yml': continue
        s=root_paths(s)
        if old in ('tools/build_release_archives.py','scripts/package_platform_release.py'):
            s=s.replace('"examples",','').replace(',"examples"','')
        if old=='.github/workflows/platform-release.yml':
            s=s.replace('"InterfaceScout/examples/README.md", ','')
        if old=='.github/workflows/consolidate-publication-release.yml':
            s=s.replace('validation/PUBLICATION_EVIDENCE.md','research/PUBLICATION_EVIDENCE.md')
        Path(new).write_text(s,encoding='utf-8')
# Explicit active indexes with no duplication; exact historical source is pinned.
manifest=[]
for m in original_manifest:
    old='datasets/'+m['group']+'/'+m['archive']
    new=mapping[old]
    assert hashlib.sha256(Path(new).read_bytes()).hexdigest()==m['sha256'],new
    with zipfile.ZipFile(new) as z: assert z.testzip() is None,new
    manifest.append(dict(m,repository_path=new))
Path('research/MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
url='https://github.com/zeynepguneryilmaz/InterfaceScout'
lines=['# Research data and reproducibility','','Two structural studies accompany InterfaceScout v1.1.0 (model 2.5.0-dual-accessibility).','',
       '| Study | Scope |','| --- | --- |','| [Method verification](method_verification/README.md) | 48 structures; numerical verification, sensitivity and robustness |','| [Protein surface analysis](protein_surface_analysis/README.md) | 76 structures; 11 maps, external-method comparisons, electrostatics and literature-based examples |','',
       'The archive index [MANIFEST.json](MANIFEST.json) records the exact source runs, commits, file counts and SHA-256 checksums. [Publication evidence](PUBLICATION_EVIDENCE.md) links the calculations. Historical screening and retrospective records are retained within the relevant study. The original pre-reorganization tree is available at [frozen repository snapshot]('+url+'/tree/'+old_head+').','',
       'The archived source scripts and workflow definitions preserve their original layout assumptions; reproduce a calculation from its recorded source commit (see the instructions in each study). The worked 1LYZ example uses the current local API.','',
       'Verify the research archive from the repository root: `python research/verify_datasets.py`.']
Path(INDEX).write_text('\n'.join(lines)+'\n')
for group,title,desc in [(M,'Method verification — 48 structures','Technical verification, numerical sensitivity, region definition robustness and earlier exploratory literature analyses.'),(P,'Protein surface analysis — 76 structures','Canonical maps, structural applicability, external method comparisons, optional electrostatics and literature-associated examples.')]:
    rows=['# '+title,'',desc,'','[Protein structure manifest](protein_structures.csv) · [Research index](../README.md) · [Source provenance](../PUBLICATION_EVIDENCE.md)','','## Archived numerical results','','| Archive | Source run |','| --- | --- |']
    for m in manifest:
        path=m['repository_path']
        if path.startswith(group+'/'): rows.append('| ['+m['archive']+']('+posixpath.relpath(path,group)+') | ['+str(m['run_id'])+']('+m['run_url']+') |')
    rows+=['','## Reproduction','','Use the source SHA in the [archive manifest](../MANIFEST.json) for the specific analysis. For example:','','```bash','git clone '+url+'.git','cd InterfaceScout','git checkout <source_sha>','python -m pip install -r backend/requirements.txt','```','','Run the command recorded in the source workflow at that commit. [Source scripts](source_scripts/) and [source records](source_records/) preserve the research history; their original path assumptions belong to the frozen checkout.']
    if group==P: rows+=['','[Worked examples](examples/README.md) use the current local interface. Literature-based examples illustrate scientific interpretation of the maps.']
    Path(group+'/README.md').write_text('\n'.join(rows)+'\n')
# Correct current example links and active documentation paths after relocation.
for f in [Path('README.md'),Path('docs/USER_GUIDE.md'),Path(P+'/examples/README.md')]:
    s=f.read_text()
    s=s.replace('research/method_verification/historical_literature/verify_publication.py','research/verify_datasets.py')
    s=s.replace('The frozen historical parameter and benchmark check is separate from the archived publication data. ','')
    s=s.replace('The frozen historical parameter and benchmark check','The research archive integrity check')
    f.write_text(s)
# New archive checker; retains meaningful assertions rather than exact README phrasing.
Path('research/verify_datasets.py').write_text('''"""Validate the two research panels and all permanent numerical archives."""
from pathlib import Path
import csv, hashlib, json, zipfile
ROOT=Path(__file__).resolve().parents[1]
R=ROOT/"research"
for group,n in [("method_verification",48),("protein_surface_analysis",76)]:
    p=R/group/"protein_structures.csv"
    with p.open(encoding="utf-8-sig",newline="") as f: rows=list(csv.DictReader(f))
    assert len(rows)==n,(group,len(rows))
    if n==76:
        ids=[r["proposed_pdb_id"].upper() for r in rows]
        assert len(set(ids))==76 and "1STP" in ids
        assert not {"3GHG","1QJS","4DN0","6LQH"}.intersection(ids)
manifest=json.loads((R/"MANIFEST.json").read_text())
assert len(manifest)==16 and len({m["repository_path"] for m in manifest})==16
for m in manifest:
    p=ROOT/m["repository_path"]
    assert hashlib.sha256(p.read_bytes()).hexdigest()==m["sha256"],p
    with zipfile.ZipFile(p) as z: assert z.testzip() is None,p
print("RESEARCH_DATASETS_PASS: 48 + 76 structures, 16 source-verified archives")
''')
# CI checks moved: historical verifier remains in its source records but not a current manuscript check.
for f in (Path('.github/workflows/publication-ci.yml'),Path('.github/workflows/platform-release-packages.yml')):
    s=f.read_text().replace('research/method_verification/historical_literature/verify_publication.py','research/verify_datasets.py')
    # Avoid duplicate archive verification in publication CI.
    s=re.sub(r'\n      - name: Frozen validation checks\n        run: python research/verify_datasets.py','',s)
    f.write_text(s)
# Ensure active docs reference actual existing files. Directory references are allowed.
for f in [Path('README.md'),Path('docs/USER_GUIDE.md'),Path(INDEX),Path(M+'/README.md'),Path(P+'/README.md'),Path(P+'/examples/README.md')]:
    s=f.read_text()
    for x in re.findall(r'\]\(([^)\s]+)\)',s):
        if x.startswith(('#','http:','https:','mailto:')): continue
        q=(f.parent/x.split('#')[0]).resolve()
        assert q.exists(),(str(f),x)
# All archived binary files are still byte-identical; original script source remains in git history.
for p in paths:
    if p not in mapping: continue
    if p.endswith(('.zip','.pdb','.png','.csv')) and p not in ('validation/strict_apo_48_selected_candidates_v2_final.csv','validation/ISV_active_candidates_76_v3.csv'):
        q=mapping[p]
        assert Path(q).exists(),q
        assert git('hash-object',q)==git('rev-parse',old_head+':'+p),(p,q)
subprocess.run(['python','research/verify_datasets.py'],check=True)
print('REORGANIZATION_VERIFIED',len(mapping),'paths')
