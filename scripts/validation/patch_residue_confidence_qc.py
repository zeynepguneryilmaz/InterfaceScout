#!/usr/bin/env python3
"""Independent post-hoc patch-confidence / precursor-maturation QC; no V1 reranking."""
import argparse,csv,json,time,urllib.request
from pathlib import Path
AA={'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
def fetch(url):
    for attempt in range(3):
        try:
            with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'InterfaceScout-QC/1.0'}),timeout=65) as r:return r.read().decode()
        except Exception:
            if attempt==2:raise
            time.sleep(2*(attempt+1))
def ca_plddt(text):
    output={}
    for line in text.splitlines():
        if line.startswith('ATOM') and line[12:16].strip()=='CA' and line[16] in (' ','A'):
            try:output[f"{line[21]}:{int(line[22:26])}:{line[26].strip()}"]=float(line[60:66])
            except (ValueError,IndexError):pass
    return output
def processing(uniprot,length):
    if len(uniprot['sequence']['value'])!=length:return None,'UNIPROT_LENGTH_MISMATCH'
    regions=[]
    for f in uniprot.get('features',[]):
        if f.get('type','').lower() in ('signal','signal peptide','propeptide','transit peptide','transit'):
            try:regions.append((int(f['location']['start']['value']),int(f['location']['end']['value']),f['type']))
            except (KeyError,ValueError,TypeError):continue
    return regions,'ANNOTATED'
def write_csv(path,rows):
    if not rows:return
    with path.open('w',newline='') as out:
        w=csv.DictWriter(out,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--input',required=True)
    parser.add_argument('--output',required=True)
    args=parser.parse_args()
    out=Path(args.output);out.mkdir(parents=True,exist_ok=True)
    cache=out/'cache';cache.mkdir(exist_ok=True)
    source=list(Path(args.input).rglob('LP*_*.json'))
    patch_rows=[];protein_rows=[];seen=set()
    for path in source:
        d=json.loads(path.read_text());m=d['large_panel_metadata'];q=m['plddt_qc']
        length=int(m['uniprot_length']);cov=float(m['af_coverage_fraction'])
        if not(q['flag']=='STANDARD' and float(q['median'])>=90 and float(q['frac_lt50'])<=.05 and .95<=cov<=1.05):continue
        pid=m['panel_id']
        if pid in seen:raise RuntimeError('Duplicate protein '+pid)
        seen.add(pid);acc=m['uniprot_accession'].split('-')[0]
        plddt={};excluded=None;annot='UNVERIFIED';err=[]
        try:
            dest=cache/(m['af_entry_id']+'.pdb')
            if not dest.exists():dest.write_text(fetch(m['af_pdb_url']))
            plddt=ca_plddt(dest.read_text())
        except Exception as e:err.append('PDB '+str(e)[:90])
        try:
            dest=cache/(acc+'.json')
            if not dest.exists():dest.write_text(fetch('https://rest.uniprot.org/uniprotkb/'+acc+'.json'))
            excluded,annot=processing(json.loads(dest.read_text()),length)
        except Exception as e:err.append('UNIPROT '+str(e)[:90])
        n_pass=0;n_mask=0
        for name,z in d['maps'].items():
            best=[x for x in z['patches'] if x['display_rank']==1]
            if len(best)!=1:raise RuntimeError('Rank1 absent or ambiguous '+pid+name)
            p=best[0];members=p['members'];vals=[plddt[x] for x in members if x in plddt]
            mask=[]
            if excluded is not None:
                for member in members:
                    n=int(member.split(':')[1])
                    if any(a<=n<=b for a,b,_ in excluded):mask.append(n)
            avg=sum(vals)/len(vals) if vals else None
            high=sum(v>=90 for v in vals)/len(vals) if vals else None
            if excluded is None:state='UNVERIFIED_PROCESSING'
            elif len(vals)!=len(members):state='MISSING_PDB_RESIDUES'
            elif mask:state='PROCESSED_REGION_OVERLAP'
            elif avg is None or avg<90 or min(vals)<70 or high<.8:state='LOW_PATCH_CONFIDENCE'
            else:state='LOCAL_QC_PASS_PAE_PENDING'
            n_pass+=state=='LOCAL_QC_PASS_PAE_PENDING';n_mask+=bool(mask)
            patch_rows.append({'panel_id':pid,'accession':acc,'chemistry':name,'rank1_center':p['center_key'],'members':';'.join(members),'mean_plddt':round(avg,3) if avg is not None else '','min_plddt':min(vals) if vals else '','fraction_ge90':round(high,3) if high is not None else '','processing_status':annot,'processed_overlap':';'.join(map(str,mask)),'qc_state':state,'io_errors':';'.join(err)})
        protein_rows.append({'panel_id':pid,'protein':m['protein_name'],'accession':acc,'n_channels':len(d['maps']),'local_pass':n_pass,'processed_overlap_channels':n_mask,'all11_local_pass':n_pass==11,'pae_checked':False,'oligomer_context_checked':False,'io_errors':';'.join(err)})
    if len(seen)!=97:raise RuntimeError('Expected 97 strict pre-screen models, found '+str(len(seen)))
    write_csv(out/'patch_residue_qc.csv',patch_rows);write_csv(out/'protein_residue_qc.csv',protein_rows)
    (out/'summary.json').write_text(json.dumps({'n_proteins':len(protein_rows),'n_rank1_patches':len(patch_rows),'local_pass_all11':sum(x['all11_local_pass'] for x in protein_rows),'note':'Post-hoc local structural QC only; PAE/oligomer and experimental localization not validated'},indent=2))
if __name__=='__main__':main()
