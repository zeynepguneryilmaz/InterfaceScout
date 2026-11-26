import csv, json, io, zipfile, hashlib, pathlib, urllib.request, urllib.parse, concurrent.futures, time, re
import numpy as np
from scipy.stats import spearmanr
out=pathlib.Path("audit");out.mkdir(exist_ok=True)
root=pathlib.Path("datasets/ISV76")
if not root.exists(): root=pathlib.Path("research/protein_surface_analysis")
def locate(name): 
    matches=list(root.rglob(name)); assert len(matches)==1,(name,matches); return matches[0]
def rows(z,n):return list(csv.DictReader(io.StringIO(z.read(n).decode("utf-8-sig"))))
p=locate("ISV76_paired_PB.zip")
records=[];allx=[];ally=[]
with zipfile.ZipFile(p) as z:
    for n in z.namelist():
        if n.endswith("/paired.csv"):
            r=rows(z,n); pid=n.split("/")[-2];x=np.array([float(a["LPBE_kT_e"]) for a in r]);y=np.array([float(a["NPBE_kT_e"]) for a in r])
            assert len(set(a["key"] for a in r))==len(r) and np.isfinite(x).all() and np.isfinite(y).all()
            a=np.abs(x-y);records.append(dict(pdb=pid,n=len(x),rho=float(spearmanr(x,y).statistic),median_lpbe=float(np.median(x)),median_npbe=float(np.median(y)),p95=float(np.quantile(a,.95)),gt1=int(sum(a>1)),gt2=int(sum(a>2)),lpbe_gt10=int(sum(abs(x)>10)),npbe_gt10=int(sum(abs(y)>10)),sign_match=int(sum(np.sign(x)==np.sign(y)))))
            allx.extend(x);ally.extend(y)
    print("PAIRED_CSV_COUNT",len(records),flush=True)
    if len(records)!=76: print("ARCHIVE_MEMBER_SAMPLE",z.namelist()[:20],flush=True)
assert len(records)==76 and len({r["pdb"] for r in records})==76
x=np.array(allx);y=np.array(ally)
res=dict(structures=len(records),pairs=len(x),median_rho=float(np.median([r["rho"] for r in records])),rho_ge095=sum(r["rho"]>=.95 for r in records),p95_gt1=sum(r["p95"]>1 for r in records),gt1=int(sum(abs(x-y)>1)),gt2=int(sum(abs(x-y)>2)),lpbe_gt10=int(sum(abs(x)>10)),npbe_gt10=int(sum(abs(y)>10)),sign_agreement=float(np.mean(np.sign(x)==np.sign(y))),case_8C50=[r for r in records if r["pdb"].upper()=="8C50"],archive_sha256=hashlib.sha256(p.read_bytes()).hexdigest())
print("PAIRED_RAW_RESULT",json.dumps(res),flush=True)
(out/"paired_PB_summary.json").write_text(json.dumps(res,indent=2))
with (out/"paired_PB_recalculated.csv").open("w") as f:w=csv.DictWriter(f,fieldnames=list(records[0]));w.writeheader();w.writerows(records)
for name in ("ISV76_PB_grid_8C50.zip","ISV76_PB_performance.zip"):
    with zipfile.ZipFile(locate(name)) as z:
        print(name,"CSV_JSON", [n for n in z.namelist() if n.endswith((".csv",".json"))][:30],flush=True)
        for n in z.namelist():
            if n.endswith(("summary.json","aggregate.json","comparison.csv")) and len(z.read(n))<200000:
                dest=out/(name.replace(".zip","")+"_"+n.replace("/","__"));dest.write_bytes(z.read(n))
dois="""10.1038/nnano.2012.207
10.1021/acs.langmuir.8b03411
10.1016/0022-2836(73)90011-9
10.1371/journal.pone.0080635
10.1093/bioinformatics/btp163
10.12688/f1000research.7931.1
10.1016/0022-2836(82)90515-0
10.1073/pnas.181342398
10.1038/s41598-018-36950-8
10.1021/acs.jcim.3c01490
10.1021/acs.jcim.1c01269
10.1016/j.colsurfb.2020.110992
10.1016/S0022-2836(05)80364-X
10.1021/acsomega.5c13030
10.1021/ar100017n
10.1016/j.jinorgbio.2008.05.006
10.1098/rsif.2016.0462
10.1038/s41467-022-34749-w
10.1021/acs.langmuir.5c06090
10.1093/nar/28.1.235
10.1016/0022-2836(91)90027-4
10.1110/ps.051840806
10.1186/s13062-016-0159-9
10.1109/4235.996017
10.1093/nar/gkm276
10.1021/ct100578z
10.1021/ja00124a002
10.1093/nar/gkae1087
10.1016/j.str.2017.10.009
10.1016/j.actbio.2015.08.033
10.1039/C3CS60075K
10.1021/acs.langmuir.3c03163
10.1021/nl2044524
10.1016/j.jbc.2026.113533
10.1021/la503854a
10.1074/jbc.M004752200
10.1074/jbc.M604592200
10.1038/srep00999
10.1073/pnas.94.23.12291
10.1073/pnas.1306390110
10.1021/acs.jpcc.6b08469
10.1126/science.1279805
10.1194/jlr.M063388
10.1093/protein/12.6.439
10.1107/S0907444912027047
10.1007/s13399-026-07231-6
10.1073/pnas.95.7.3413
10.1016/j.jmb.2008.12.057
10.2210/pdb3ipj/pdb
10.1107/S1399004714004854
10.1371/JOURNAL.PPAT.1004228
10.1107/S0909049513020803
10.2210/pdb6ndq/pdb
10.1371/JOURNAL.PONE.0098042
10.1111/febs.16360
10.2210/pdb3eoi/pdb
10.1074/jbc.M808305200
10.1016/j.ibmb.2022.103785
10.1016/j.str.2010.03.006
10.1101/gad.187815.112
10.1093/nar/gks177
10.1038/srep35552
10.1046/j.1365-2958.2002.02809.x
10.1073/PNAS.97.12.6306
10.1016/j.molcel.2006.02.019
10.1016/J.JMB.2003.09.060
10.1023/A:1021257701676
10.2210/pdb2fb5/pdb
10.1074/jbc.M302926200
10.2210/pdb6inx/pdb
10.1016/j.bbrc.2010.04.130
10.1126/science.279.5353.1048
10.1016/j.bbagen.2019.129418
10.1038/32210
10.1021/acs.jafc.9b01290
10.1016/j.bbrc.2013.07.051
10.2210/pdb2pcn/pdb
10.1074/jbc.RA119.009713
10.1074/jbc.M504642200
10.3390/ijms22063219
10.2210/pdb3uci/pdb
10.1172/JCI131842
10.1074/jbc.M110.205161
10.1107/S2053230X13034407
10.1074/jbc.RA118.004012
10.1038/s41589-019-0415-2
10.1021/bi049071r
10.2210/pdb2dc4/pdb
10.1021/bi047711k
10.1016/S0969-2126(97)00192-5
10.1038/NCOMMS6701
10.1016/j.jmb.2004.04.001
10.1038/s41467-023-40787-9
10.1186/s13567-021-00969-z
10.1016/S0969-2126(99)80023-9
10.1093/emboj/18.22.6228
10.1016/S0022-2836(02)01382-7
10.1006/jmbi.2000.4078
10.1038/86188
10.1021/bi962514+
10.1110/ps.29002
10.1016/0014-5793(95)01331-8
10.1007/s007750050018
10.1107/S0907444994013521
10.1126/science.2063194
10.1038/nsb971
10.1038/367338a0
10.1038/nature04005
10.1038/nsmb1210
10.1107/S0907444911015423
10.1038/nature10823
10.1016/j.celrep.2012.02.012
10.1073/pnas.2200102119
10.1016/S0969-2126(99)80051-3
10.1126/sciadv.aax5068
10.1107/S0907444900006387
10.1021/bi973000t
10.1016/S0006-3495(99)77056-6
10.1006/jmbi.1993.1321
10.1016/S0969-2126(97)00205-0
10.1107/S0907444995008675
10.1016/0022-2836(92)90324-D
10.1038/nbt1096-1246
10.1016/0022-2836(89)90276-3
10.1016/S0021-9258(19)50095-3
10.1016/0022-2836(90)90200-6
10.1016/0022-2836(91)80185-W
10.1126/science.2911722
10.1016/0022-2836(87)90679-6
10.1016/S0969-2126(97)00167-6
10.1016/0022-2836(92)90531-N
10.1006/jmbi.1997.0970
10.1021/bi060773l
10.1021/ja0350684
10.1016/0022-2836(90)90313-B
10.2210/pdb7kwi/pdb
10.1021/bi00408a010
10.1016/0022-2836(91)90724-K
10.1126/science.8036495
10.1107/S0907444999003431
10.1006/jmbi.1999.3458
10.1021/bi00443a038
10.1006/jmbi.1994.0061
10.1074/JBC.M204571200
10.1038/417463A
10.1021/bi011727i
10.1016/0022-2836(87)90293-2
10.2210/pdb1s0q/pdb
10.1107/S0907444995014119
10.1002/prot.340100203
10.1016/S0969-2126(00)00031-9
10.1007/s007750050231
10.1016/0022-2836(92)90260-Q
10.1016/0022-2836(82)90174-7
10.2210/pdb5b1e/pdb
10.1016/0022-2836(84)90467-4
10.1038/NSB865
10.1016/S0092-8674(03)00809-2
10.1128/JB.00374-13
10.1074/jbc.M115.674929
10.1073/pnas.1718102115
10.2210/pdb5tcb/pdb
10.1073/pnas.1702798114
10.1074/jbc.RA119.008335
10.1038/s41467-023-38204-2
10.1073/pnas.2534280123""".splitlines()
def resolve(doi):
    api="https://api.datacite.org/dois/" if doi.lower().startswith("10.2210/") else "https://api.crossref.org/works/"
    for i in range(3):
        try:
            req=urllib.request.Request(api+urllib.parse.quote(doi,safe=""),headers={"User-Agent":"InterfaceScout publication DOI audit (GitHub Actions)","Accept":"application/json"})
            with urllib.request.urlopen(req,timeout=20) as r: d=json.load(r)
            if "datacite" in api:
                m=d["data"]["attributes"];return dict(doi=doi,status="resolved",title=m["titles"][0]["title"],year=m.get("publicationYear"),first_author=m.get("creators",[{}])[0].get("name",""),source="DataCite",journal=m.get("publisher",""),volume="",page="")
            m=d["message"];date=m.get("published",m.get("issued",{})).get("date-parts",[[None]])[0]
            return dict(doi=doi,status="resolved",title=(m.get("title") or [""])[0],year=date[0],first_author=(m.get("author") or [{}])[0].get("family",""),source="Crossref",journal=(m.get("container-title") or [""])[0],volume=m.get("volume",""),page=m.get("page",m.get("article-number","")))
        except Exception as e: err=str(e);time.sleep(i+1)
    return dict(doi=doi,status="unresolved",title="",year="",first_author="",source=err,journal="",volume="",page="")
with concurrent.futures.ThreadPoolExecutor(max_workers=5) as ex: refs=list(ex.map(resolve,dois))
with (out/"doi_metadata.csv").open("w",encoding="utf-8",newline="") as f:w=csv.DictWriter(f,fieldnames=list(refs[0]));w.writeheader();w.writerows(refs)
print("DOI_RESOLUTION",len(refs),sum(r["status"]=="resolved" for r in refs),flush=True)
for r in refs: print("DOI_METADATA",json.dumps(r,ensure_ascii=False),flush=True)
