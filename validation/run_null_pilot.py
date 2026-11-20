import csv, json, sys, numpy as np
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"backend"))
import core
from prepare import prepare_pdb_text
from interface_engine import analyze_all_maps
from geometry import extract_ca_nodes,build_surface_geometry,ca_distance,same_face
from coarse_patch import build_coarse_patches
ROOT=Path(__file__).resolve().parents[1]
maps={"silica":"anionic","paa_fe3o4":"anionic","citrate_au":"anionic","mpa_au":"anionic","alumoh":"cationic","polystyrene":"hydrophobic"}
cases=json.loads((ROOT/"validation/benchmark_manifest.json").read_text())["cases"]
frozen={r["case_id"]:r for r in csv.DictReader((ROOT/"validation/benchmark_case_results.csv").open(encoding="utf-8-sig"))}
def recall(patches,nodes,truth,k,threshold,equiv):
    lookup={n["key"]:i for i,n in enumerate(nodes)}
    members={lookup[m] for p in patches[:k] for m in p["members"] if m in lookup}
    if not members:return 0.
    coords=np.array([n["coord"] for n in nodes]); truth_idx=np.array(truth)
    dist=np.linalg.norm(coords[truth_idx,None,:]-coords[list(members)][None,:,:],axis=-1)
    hit=np.any(dist<=threshold+1.e-9,axis=1)
    if equiv:
        chains=np.array([nodes[i]["chain"] for i in truth])
        return max(float(hit[chains==c].mean()) for c in set(chains))
    return float(hit.mean())
def truth_indices(case,nodes):
    out=[]
    for i,n in enumerate(nodes):
        if n["chain"] not in case["chain"].split(","):continue
        seq=n["res_seq"]
        if case["gt_kind"]=="exact":
            if seq in case["gt_exact"]:out.append(i)
        elif any(a<=seq<=b for a,b in case["gt_ranges"]):out.append(i)
    return out
def pvalue(arr,obs):
    return (1+sum(v>=obs-1e-12 for v in arr))/(len(arr)+1)
rng=np.random.default_rng(20261116)
N=int(sys.argv[1]) if len(sys.argv)>1 else 100
for case in cases:
    cid=case["id"]; pdb=case["pdb_id"]; ph=case["pH"]
    files=sorted((ROOT/"examples"/pdb).rglob("*_Bfactor.pdb"))
    if pdb=="2CBA":files=[p for p in files if "pH_"+str(ph) in str(p)]
    if len(files)!=1:print(cid,"INPUT_MISSING",len(files),flush=True);continue
    prepared,_=prepare_pdb_text(files[0].read_text(),chain=case["chain"])
    nodes=extract_ca_nodes(prepared)
    truth=truth_indices(case,nodes)
    if not truth:print(cid,"NO_TRUTH",flush=True);continue
    res=analyze_all_maps(pdb_text=prepared,chain=None,pH=ph)
    orig=res["maps"][maps[case["surface"]]]["patches"]
    equiv=case["gt_kind"]=="equivalent_chain_ranges"
    observed=recall(orig,nodes,truth,5,8.,equiv)
    expected=float(frozen[cid]["top5_near_8A_recall"])
    if abs(observed-expected)>0.015:
        print(cid,"BASELINE_MISMATCH",observed,expected,flush=True);continue
    v1=core.analyze(core.AnalyzeRequest(pdb_text=prepared,env=core.EnvParams(pH=ph,ionic=150,temp=298)))
    surf=v1["surface_residues"];dist=core.build_distances(surf)
    geom=build_surface_geometry(v1,nodes);keys=list(geom["coords"])
    library=[[k for k in keys if ca_distance(c,k,geom)<=8 and same_face(c,k,geom)] for c in keys]
    sizes=np.array([len(p) for p in library])
    null={"geometry":[],"chemistry":[],"ranking":[]}
    for j in range(N):
        used=set();random_patches=[]
        for pp in orig[:5]:
            dif=np.abs(sizes-len(pp["members"])).astype(float)
            if used:dif[list(used)]=np.inf
            possibilities=np.where(dif==dif.min())[0]
            idx=int(rng.choice(possibilities)); used.add(idx)
            random_patches.append({"members":library[idx]})
        null["geometry"].append(recall(random_patches,nodes,truth,5,8.,equiv))
        shuffled=[dict(x,res_name=surf[int(v)]["res_name"]) for x,v in zip(surf,rng.permutation(len(surf)))]
        channel=core.chemistry_map(shuffled,dist,maps[case["surface"]],ph)
        altered=dict(v1,chemistries={maps[case["surface"]]:channel})
        chemical=build_coarse_patches(v1_result=altered,chemistry=maps[case["surface"]],ca_nodes=nodes)
        null["chemistry"].append(recall(chemical,nodes,truth,5,8.,equiv))
        reordered=[orig[int(i)] for i in rng.permutation(len(orig))]
        null["ranking"].append(recall(reordered,nodes,truth,5,8.,equiv))
    print(cid,"OK","observed",round(observed,5),"null",[(name,round(float(np.mean(vals)),5),round(pvalue(vals,observed),5)) for name,vals in null.items()],flush=True)
