/* Additive multi-PDB/consensus/complete-export UI. All original controls remain. */
let ensembleInputs=[], ensembleAnalyzedInputs=[], ensembleResult=null, ensembleActiveIndex=0, ensembleBusy=false;
let ensembleConsensusKey=null, ensembleSavedSnapshot=null;
let ensembleDownloadUrl=null;
function ensembleReadyLink(){
  let link=$("isReadyDownload");
  if(!link){
    link=document.createElement("a");
    link.id="isReadyDownload";
    link.className="btn secondary";
    link.style.cssText="display:none;text-align:center;text-decoration:none;margin-top:10px;font-size:13px";
    link.textContent="Save completed ZIP";
    const button=$("isFullExport");
    button.insertAdjacentElement("afterend",link);
  }
  return link;
}
function ensembleClearDownload(){
  const link=ensembleReadyLink();
  link.style.display="none";
  link.removeAttribute("href");
  ensembleDownloadUrl=null;
}


function ensembleSetup(){
  let css=document.createElement("style");
  css.textContent=".is-set{border:1px solid #dfe4ec;border-radius:8px;background:#fff;padding:10px;margin:8px 0}.is-set p{margin:6px 0;font-size:12px;color:#667085}.is-inline{display:flex;gap:6px;flex-wrap:wrap}.is-inline button{width:auto}.is-meta{font-size:11px;color:#667085;line-height:1.4;margin:6px 0}.is-choice{border:1px solid #dfe4ec;border-radius:6px;background:#fff;padding:7px;min-width:0;flex:1}.is-choice:focus{outline:2px solid #6d28d9}.is-tabs{display:flex;gap:7px;flex-wrap:wrap;margin:8px 0}.is-tabs button{background:white;border:1px solid #dfe4ec;border-radius:6px;padding:6px 9px;cursor:pointer}.is-tabs button.on{border-color:#6d28d9;background:#f5f3ff;color:#6d28d9}.is-errors{color:#b42318;font-size:12px;white-space:pre-wrap}.is-progress{font-size:11px;color:#2563eb;margin:8px 0}";
  document.head.appendChild(css);
  const place=$("run").closest("aside")||document.querySelector("aside");
  place.insertAdjacentHTML("beforeend",
    '<div class="section"><h3>Multi-PDB analysis (no trajectory)</h3>'+
    '<div class="note">One or several PDB structures. Each structure retains all 11 original maps. For multiple matching structures, compare and report descriptive residue consensus.</div>'+
    '<div class="field" style="margin-top:9px"><label>Upload one or several PDB files</label>'+
    '<input id="isMultiFiles" class="input" type="file" accept=".pdb" multiple onchange="ensembleLoadFiles(this.files)"></div>'+
    '<div class="is-inline" style="margin:7px 0">'+
    '<button class="workAction" onclick="ensembleAddCurrent()">Add current PDB</button>'+
    '<button class="workAction" onclick="ensembleClearQueue()">Clear selection</button></div>'+
    '<div class="field"><label>Or add by PDB ID</label><div class="row" style="grid-template-columns:1fr 80px">'+
    '<input class="input" id="isAnotherId" placeholder="e.g. 1AKE">'+
    '<button class="btn secondary" onclick="ensembleFetchId()">Add</button></div></div>'+
    '<div id="isQueued" class="is-meta">No structures added yet.</div>'+
    '<button class="btn secondary" onclick="ensembleAnalyzeSet()" id="isAnalyzeSet">Analyze selected PDBs</button>'+
    '<div class="is-progress" id="isStatus"></div><div class="is-errors" id="isError"></div>'+
    '</div>'+
    '<div class="section"><h3>Complete export</h3>'+
    '<button class="btn secondary" id="isFullExport" onclick="ensembleExportAll()" disabled>Download all PDFs, images, CSV, Excel and B-factor PDBs (ZIP)</button>'+
    '<div class="note" style="margin-top:8px">Exports all 11 colored map screenshots for every structure; includes full report, all-residue data, and optional comparison/consensus. Original Excel and current B-factor export stay available.</div></div>'
  );
  $("results").insertAdjacentHTML("afterend",'<section id="isSetResults" class="results" style="padding-top:0"></section>');
  ensembleRefresh();
}
function ensembleMessage(text,error=false){(error?$("isError"):$("isStatus")).textContent=text;}
function ensembleRefresh(){
  let root=$("isQueued");
  if(!root)return;
  root.innerHTML=ensembleInputs.length?
    ensembleInputs.map((x,i)=>'<div class="is-set">'+esc(x.name)+' · '+esc(x.chain||"all chains")+
      ' <button class="workAction" onclick="ensembleRemove('+i+')">Remove</button></div>').join("")+
      '<p>'+ensembleInputs.length+' PDB structures queued (maximum 20).</p>':
    "No structures added yet. You may export a single PDB from the existing analysis.";
  $("isFullExport").disabled=!(data||ensembleResult)||ensembleBusy;
  $("isAnalyzeSet").disabled=ensembleBusy||!ensembleInputs.length;
  if(ensembleResult)ensembleShowSetResults();
}
function ensembleUniqueName(raw){
  let name=String(raw||"structure").replace(/\.pdb$/i,"").trim()||"structure";
  const existing=new Set(ensembleInputs.map(x=>x.name));
  if(!existing.has(name))return name;
  let i=2;while(existing.has(name+"_"+i))i++;
  return name+"_"+i;
}
function ensembleAddInput(spec){
  if(ensembleInputs.length>=20)throw Error("Maximum 20 PDB structures per analysis");
  if(!spec.pdb_text||!spec.pdb_text.includes("ATOM"))throw Error("The file has no protein ATOM records");
  ensembleInputs.push({name:ensembleUniqueName(spec.name),pdb_text:spec.pdb_text,chain:spec.chain||null});
}
async function ensembleLoadFiles(filelist){
  ensembleMessage("");ensembleMessage("",true);
  try{
    for(const f of Array.from(filelist||[])){
      if(!/\.pdb$/i.test(f.name))throw Error("Only PDB files are supported for this release: "+f.name);
      ensembleAddInput({name:f.name,pdb_text:await f.text(),chain:$("chain").value.trim()||null});
    }
    ensembleRefresh();
  }catch(e){ensembleMessage(e.message,true)}
}
function ensembleAddCurrent(){
  try{
    if(!pdbText)throw Error("Load or fetch the first PDB using the original Protein input");
    ensembleAddInput({name:structureName,pdb_text:pdbText,chain:$("chain").value.trim()||null});
    ensembleMessage("",true);ensembleRefresh();
  }catch(e){ensembleMessage(e.message,true)}
}
async function ensembleFetchId(){
  const id=$("isAnotherId").value.trim().toUpperCase();
  try{
    if(!/^[0-9][A-Z0-9]{3}$/.test(id))throw Error("Enter a four-character PDB ID");
    ensembleMessage("Retrieving "+id+"...");
    const res=await fetch("https://files.rcsb.org/download/"+id+".pdb");
    if(!res.ok)throw Error("Unable to retrieve PDB "+id);
    ensembleAddInput({name:id,pdb_text:await res.text(),chain:$("chain").value.trim()||null});
    ensembleRefresh();ensembleMessage("Added "+id);
    ensembleMessage("",true);
  }catch(e){ensembleMessage(e.message,true)}
}
function ensembleRemove(i){ensembleInputs.splice(i,1);ensembleRefresh()}
function ensembleClearQueue(){ensembleInputs=[];ensembleRefresh()}
async function ensembleAnalyzeSet(){
  if(ensembleBusy||!ensembleInputs.length)return;
  ensembleBusy=true;ensembleRefresh();ensembleMessage("",true);
  ensembleMessage("Calculating 11 maps per PDB and matched-residue consensus...");
  const body={structures:ensembleInputs,pH:+$("pH").value,ionic_mM:+$("ionic").value,
    temp_K:+$("temp").value,combined_maps:selected.size>=2?[...selected]:null};
  try{
    const res=await fetch("/structure-set",{method:"POST",headers:{"Content-Type":"application/json"},
      body:JSON.stringify(body)});
    const obj=await res.json();
    if(!res.ok)throw Error(obj.detail||"Unable to analyze the structure set");
    ensembleResult=obj;
    ensembleAnalyzedInputs=ensembleInputs.map(x=>({...x}));
    ensembleActiveIndex=0;
    ensembleActivate(0);
    ensembleMessage("Analyzed "+obj.n_structures+" PDB structures; all 11 maps retained for every structure.");
  }catch(e){ensembleMessage(e.message,true)}
  finally{ensembleBusy=false;ensembleRefresh()}
}
function ensembleActivate(index){
  if(!ensembleResult||!ensembleResult.structures[index])return;
  const x=ensembleResult.structures[index];
  ensembleActiveIndex=index;ensembleConsensusKey=null;
  pdbText=ensembleAnalyzedInputs[index].pdb_text;structureName=x.name;data=x.analysis;
  $("chain").value=ensembleAnalyzedInputs[index].chain||"";
  $("excelBtn").disabled=false;$("pdbBtn").disabled=false;
  activeType="chemistry";activeMap=data.initial_map;focusedPatchRank=null;
  decisionCompare=null;decisionCombo=x.combined_maps||null;decisionComboError="";
  if(!selected.size)selected=new Set([activeMap]);
  render();ensembleRefresh();
}
function ensembleShowSetResults(){
  const root=$("isSetResults");
  if(!root||!ensembleResult)return;
  const runs=ensembleResult.structures;
  const buttons=runs.map((r,i)=>'<button class="'+(ensembleActiveIndex===i&&!ensembleConsensusKey?'on':'')+
    '" onclick="ensembleActivate('+i+')">'+esc(r.name)+'</button>').join("");
  let html='<div class="is-set"><h3 style="font-size:14px">Structure results</h3>'+
    '<p>Each PDB has all 11 canonical maps and the original Excel / B-factor exports.</p>'+
    '<div class="is-tabs">'+buttons+'</div></div>';
  const c=ensembleResult.consensus;
  if(!c && ensembleResult.consensus_unavailable_reason){
    html+='<div class="is-set"><h3 style="font-size:14px">Consensus not available</h3>'+
      '<p>'+esc(ensembleResult.consensus_unavailable_reason)+
      '</p><p>Each input PDB remains independently analyzed; no artificial consensus was calculated.</p></div>';
  }
  if(c){
    html+='<div class="is-set"><h3 style="font-size:14px">Structure consensus</h3>'+
      '<p>Only matched residue identities are compared. 100% means present in all supplied structures; this is not binding probability.</p>'+
      '<span class="is-meta">'+c.n_matched_residues+' matched residues across '+c.n_structures+' structures.</span>'+
      '<div class="is-tabs">'+c.map_order.map(k=>'<button class="'+(ensembleConsensusKey===k?'on':'')+
      '" onclick="ensembleShowConsensus(\''+k+'\')">'+esc(c.maps[k].label)+'</button>').join("")+'</div>'+
      '<table><thead><tr><th>Map</th><th>Strict Top-1</th><th>Strict Top-3</th><th>Strict Top-5</th><th>Majority Top-5</th></tr></thead><tbody>'+
      c.map_order.map(k=>{let m=c.maps[k];return '<tr><td>'+esc(m.label)+'</td><td>'+
        m.strict_top1_members.length+'</td><td>'+m.strict_top3_members.length+'</td><td>'+
        m.strict_top5_members.length+'</td><td>'+m.majority_top5_members.length+'</td></tr>'}).join("")+
      '</tbody></table>'+
      '<div class="is-meta">Strict and majority are residue-set membership, not necessarily spatially contiguous consensus patches.</div></div>';
  }
  if(ensembleResult.comparison.length){
    html+='<div class="is-set"><h3 style="font-size:14px">Pairwise structure comparisons</h3>';
    for(const p of ensembleResult.comparison){
      html+='<details><summary>'+esc(p.structure_a)+' / '+esc(p.structure_b)+
        ' · '+p.accessibility.n_exposure_status_changes+' exposure-class changes</summary>'+
        '<table><thead><tr><th>Map</th><th>Top-1 Jaccard</th><th>Top-3</th><th>Top-5</th></tr></thead><tbody>'+
        p.map_comparisons.map(m=>'<tr><td>'+esc(m.label)+'</td>'+
          m.top_k.map(row=>'<td>'+decisionFormat(row.jaccard)+'</td>').join("")+'</tr>').join("")+
        '</tbody></table></details>';
    }
    html+='</div>';
  }
  if((ensembleResult.comparison_skipped||[]).length){
    html+='<div class="is-set"><h3 style="font-size:14px">Unmatched structure comparisons</h3>'+
      ensembleResult.comparison_skipped.map(function(x){
        return '<p>'+esc(x.structure_a)+' / '+esc(x.structure_b)+': '+esc(x.reason)+'</p>';
      }).join("")+'</div>';
  }
  root.innerHTML=html;
}
function ensembleShowConsensus(k){
  if(!ensembleResult||!ensembleResult.consensus)return;
  ensembleConsensusKey=k;
  const c=ensembleResult.consensus.maps[k];
  const byKey=new Map(c.residues.map(r=>[r.key,r.top5_membership_fraction]));
  rebuildViewer(false);
  for(const r of (data.residue_profiles||[])){
    const frequency=byKey.get(r.key)||0;
    if(frequency<=0)continue;
    const color="hsl("+(228-207*frequency)+",82%,"+(49+8*frequency)+"%)";
    viewer.addStyle({chain:r.chain,resi:r.res_seq},{stick:{color,radius:.24},
      sphere:{color,radius:.39}});
  }
  $("legend").innerHTML="<b>Top-5 consensus / "+esc(c.label)+
    "</b><br>Blue: low structure frequency<br>Orange/red: high structure frequency"+
    "<br>Input-PDB frequency, NOT adsorption probability";
  $("legend").style.display="block";
  viewer.render();ensembleShowSetResults();
}
async function ensembleCaptureScreenshot(){
  if(!viewer||typeof viewer.pngURI!=="function")
    throw Error("3D viewer cannot capture PNG; browser/export not supported");
  viewer.render();
  await new Promise(r=>requestAnimationFrame(()=>requestAnimationFrame(r)));
  await new Promise(r=>setTimeout(r,110));
  const uri=viewer.pngURI();
  if(!uri||typeof uri!=="string"||!uri.startsWith("data:image/png;base64,"))
    throw Error("3D screenshot unavailable; refusing to create an incomplete visual report");
  return uri;
}
function ensembleCenterAtoms(){
  if(!viewer||typeof viewer.selectedAtoms!=="function")return [];
  const atoms=viewer.selectedAtoms({atom:"CA"})||[];
  return atoms.filter(a=>Number.isFinite(a.x)&&Number.isFinite(a.y)&&Number.isFinite(a.z));
}
function ensemblePatchCenters(key){
  if(!data||!viewer||typeof viewer.selectedAtoms!=="function")return [];
  const result=[];
  for(const patch of (data.maps[key]?.patches||[]).slice(0,3)){
    const residue=patch.center_residue||
      (patch.member_residues||[]).find(r=>r.key===patch.center_key);
    if(!residue)continue;
    let atoms=viewer.selectedAtoms({chain:residue.chain,resi:Number(residue.res_seq),atom:"CA"});
    if(!atoms?.length)atoms=viewer.selectedAtoms({chain:residue.chain,resi:Number(residue.res_seq)});
    if(atoms?.length && Number.isFinite(atoms[0].x))
      result.push({patch,atom:atoms[0]});
  }
  return result;
}
function ensembleNeedsOppositeView(centers,allAtoms){
  // Only provide another orientation when top candidate centers occupy
  // opposing hemispheres. Both images represent the SAME input structure.
  if(centers.length<2||allAtoms.length<3)return false;
  const c={x:0,y:0,z:0};
  allAtoms.forEach(a=>{c.x+=a.x;c.y+=a.y;c.z+=a.z});
  c.x/=allAtoms.length;c.y/=allAtoms.length;c.z/=allAtoms.length;
  const vectors=centers.map(({atom})=>({x:atom.x-c.x,y:atom.y-c.y,z:atom.z-c.z}));
  for(let i=0;i<vectors.length;i++)for(let j=i+1;j<vectors.length;j++){
    const a=vectors[i],b=vectors[j];
    const dot=a.x*b.x+a.y*b.y+a.z*b.z;
    const na=Math.hypot(a.x,a.y,a.z),nb=Math.hypot(b.x,b.y,b.z);
    if(na>0&&nb>0&&dot/(na*nb)<-0.35)return true;
  }
  return false;
}
function ensembleExportPatchLabels(centers,allAtoms){
  const objects=[];
  if(!viewer||typeof viewer.addLabel!=="function")return objects;
  const c={x:0,y:0,z:0};
  allAtoms.forEach(a=>{c.x+=a.x;c.y+=a.y;c.z+=a.z});
  const n=Math.max(allAtoms.length,1);
  c.x/=n;c.y/=n;c.z/=n;
  for(const {patch,atom} of centers){
    const vec={x:atom.x-c.x,y:atom.y-c.y,z:atom.z-c.z};
    const norm=Math.hypot(vec.x,vec.y,vec.z)||1;
    const factor=4.5/norm;
    const target={x:atom.x+vec.x*factor,y:atom.y+vec.y*factor,z:atom.z+vec.z*factor};
    try{
      // The leader originates at the actual candidate center atom.
      if(typeof viewer.addLine==="function"){
        const shape=viewer.addLine({
          start:{x:atom.x,y:atom.y,z:atom.z},end:target,
          color:"#6B4355",linewidth:2
        });
        if(shape)objects.push({kind:"shape",value:shape});
      }
      const label=viewer.addLabel("P"+patch.display_rank,{
        position:target,fontColor:"white",fontSize:12,
        backgroundColor:"#6B4355",backgroundOpacity:.93,
        showBackground:true,borderThickness:0
      });
      if(label)objects.push({kind:"label",value:label});
    }catch(e){/* Preserve scientific screenshot even if annotations unsupported. */}
  }
  return objects;
}
function ensembleRemovePatchLabels(objects){
  objects.forEach(o=>{
    try{
      if(o.kind==="shape")viewer.removeShape(o.value);
      else viewer.removeLabel(o.value);
    }catch(e){}
  });
}
async function ensembleCaptureEleven(){
  const shots={},alternates={};
  for(let i=0;i<data.map_order.length;i++){
    const key=data.map_order[i];
    activeType="chemistry";activeMap=key;focusedPatchRank=null;
    paintCurrent(false);
    const centers=ensemblePatchCenters(key),atoms=ensembleCenterAtoms();
    const labels=ensembleExportPatchLabels(centers,atoms);
    const view=typeof viewer.getView==="function"?viewer.getView():null;
    try{
      ensembleMessage("Capturing actual 3D map image "+(i+1)+"/11 for "+structureName+"...");
      shots[key]=await ensembleCaptureScreenshot();
      if(ensembleNeedsOppositeView(centers,atoms) && view &&
         typeof viewer.rotate==="function" && typeof viewer.setView==="function"){
        try{
          viewer.rotate(180,"y");
          alternates[key]=await ensembleCaptureScreenshot();
        }catch(error){
          // An optional secondary angle must never prevent the 11 canonical
          // screenshots / PDFs / XLSX / PDB files from downloading.
          console.warn("Optional reverse view skipped for",key,error);
        }finally{
          viewer.setView(view);
          viewer.render();
        }
      }
    }finally{
      ensembleRemovePatchLabels(labels);
    }
  }
  return {screenshots:shots,alternate_screenshots:alternates};
}

async function ensembleExportAll(){
  ensembleClearDownload();
  if(ensembleBusy)return;
  if(!data&&!ensembleResult){ensembleMessage("Analyze a PDB first",true);return}
  const backup={
    pdbText,structureName,data,activeType,activeMap,
    selected:new Set(selected),focusedPatchRank,ensembleActiveIndex,ensembleConsensusKey,
    decisionCombo,decisionCompare,decisionComboError,chain:$("chain").value,
  };
  ensembleBusy=true;ensembleRefresh();ensembleMessage("",true);
  try{
    const source=ensembleResult?ensembleResult.structures:
      [{name:structureName,analysis:data,combined_maps:decisionCombo||null}];
    const exported=[];
    for(let i=0;i<source.length;i++){
      if(ensembleResult){ensembleActivate(i)}else{activeType="chemistry"}
      const {screenshots,alternate_screenshots}=await ensembleCaptureEleven();
      // The interactive panel may show only the first 50 regions. Export
      // always fetches the entire co-location set for reproducibility.
      let combined=null;
      let combined_screenshot=null;
      if(selected.size>=2){
        const resp=await fetch("/combine",{method:"POST",
          headers:{"Content-Type":"application/json"},
          body:JSON.stringify({analysis:data,selected_maps:[...selected],max_regions:null})});
        const value=await resp.json();
        if(!resp.ok)throw Error(value.detail||"Cannot combine maps");
        combined=value;
        activeType="composite";activeMap=null;paintCurrent(false);
        combined_screenshot=await ensembleCaptureScreenshot();
      }
      exported.push({
        name:source[i].name, analysis:data,
        pdb_text:ensembleResult?ensembleAnalyzedInputs[i].pdb_text:pdbText,
        screenshots, alternate_screenshots, combined_maps:combined,combined_screenshot,
      });
    }
    const consensus=ensembleResult?.consensus||null;
    const comparison=ensembleResult?.comparison||[];
    const consensus_screenshots={};
    if(consensus){
      ensembleActivate(0);
      for(let i=0;i<consensus.map_order.length;i++){
        const k=consensus.map_order[i];
        ensembleShowConsensus(k);
        ensembleMessage("Capturing consensus image "+(i+1)+"/11...");
        consensus_screenshots[k]=await ensembleCaptureScreenshot();
      }
    }
    ensembleMessage("Assembling PDF reports, every CSV, Excel, images and B-factor PDBs...");
    // Generate on the local Python server. A regular GET attachment works
    // when Chrome/Edge block delayed javascript Blob downloads.
    const result=await fetch("/export-package/prepare",{method:"POST",
      headers:{"Content-Type":"application/json"},
      body:JSON.stringify({structures:exported,consensus,comparison,consensus_screenshots,
        comparison_skipped:ensembleResult?.comparison_skipped||[],
        consensus_unavailable_reason:ensembleResult?.consensus_unavailable_reason||null})});
    if(!result.ok){
      let details;
      try{details=(await result.json()).detail}catch(e){details="Export failed"}
      throw Error(typeof details==="string"?details:JSON.stringify(details));
    }
    const prepared=await result.json();
    if(!prepared.download_url || !prepared.download_url.startsWith("/export-package/download/") ||
       !/^[0-9a-f]{32}$/.test(prepared.download_url.split("/").pop()))
      throw Error("The server returned an invalid download link.");
    const link=ensembleReadyLink();
    link.href=prepared.download_url;
    link.download=prepared.filename||"InterfaceScout_Complete_Results.zip";
    link.style.display="block";
    link.textContent="Download ZIP now · "+(prepared.bytes/1048576).toFixed(1)+" MB";
    link.setAttribute("aria-label","Download complete research ZIP");
    // A server-served attachment does not rely on a transient blob URL and
    // can be downloaded again with an explicit click for up to one hour.
    window.location.assign(prepared.download_url);
    link.scrollIntoView({behavior:"smooth",block:"nearest"});
    ensembleMessage("Complete ZIP prepared on your computer ("+
      (prepared.bytes/1048576).toFixed(1)+
      " MB). Your browser should start the download. Otherwise click ‘Download ZIP now’ above.");
  }catch(e){ensembleMessage(e.message,true)}
  finally{
    ensembleBusy=false;
    pdbText=backup.pdbText;structureName=backup.structureName;data=backup.data;
    activeType=backup.activeType;activeMap=backup.activeMap;
    selected=backup.selected;focusedPatchRank=backup.focusedPatchRank;
    ensembleActiveIndex=backup.ensembleActiveIndex;
    ensembleConsensusKey=backup.ensembleConsensusKey;
    decisionCombo=backup.decisionCombo;decisionCompare=backup.decisionCompare;
    decisionComboError=backup.decisionComboError;$("chain").value=backup.chain;
    if(data){
      render();
      if(ensembleConsensusKey)ensembleShowConsensus(ensembleConsensusKey);
    }
    ensembleRefresh();
  }
}
ensembleSetup();

// Changing the original single-PDB input must invalidate the old multi-PDB result;
// otherwise the 'complete export' button can silently export stale structures.
const ensembleOriginalResetAnalysis=resetAnalysisForNewProtein;
resetAnalysisForNewProtein=function(){
  ensembleResult=null;ensembleAnalyzedInputs=[];ensembleConsensusKey=null;
  ensembleOriginalResetAnalysis();
  const old=$("isSetResults");if(old)old.innerHTML="";
  ensembleRefresh();
};

// The original Analyze button is asynchronous. Enable complete ZIP export
// only after the canonical single-structure analysis has completed.
const ensembleOriginalAnalyze=analyze;
analyze=async function(){
  await ensembleOriginalAnalyze();
  ensembleRefresh();
};

// Changing selected chemistry channels after initial multi-PDB analysis must
// also refresh exact/spatial intersection consensus without rerunning the core.
let ensembleConsensusRefreshId=0;
async function ensembleRefreshSelectedConsensus(){
  if(!ensembleResult||!ensembleResult.consensus||ensembleBusy)return;
  const requestId=++ensembleConsensusRefreshId;
  const maps=selected.size>=2?[...selected]:[];
  try{
    const response=await fetch("/derive-consensus",{
      method:"POST",headers:{"Content-Type":"application/json"},
      body:JSON.stringify({
        structures:ensembleResult.structures.map(x=>({name:x.name,analysis:x.analysis})),
        selected_maps:maps,
      }),
    });
    const result=await response.json();
    if(!response.ok)throw Error(result.detail||"Cannot refresh consensus");
    if(requestId!==ensembleConsensusRefreshId||!ensembleResult)return;
    ensembleResult.consensus=result;
    ensembleMessage("",true);
    if(ensembleConsensusKey)ensembleShowConsensus(ensembleConsensusKey);
    else ensembleShowSetResults();
  }catch(error){
    if(requestId===ensembleConsensusRefreshId)
      ensembleMessage("Map consensus refresh failed: "+error.message,true);
  }
}
const ensemblePreviousToggleChem=toggleChem;
toggleChem=function(key){
  ensemblePreviousToggleChem(key);
  ensembleRefreshSelectedConsensus();
};
const ensemblePreviousApplyGoal=decisionApplyGoal;
decisionApplyGoal=function(){
  ensemblePreviousApplyGoal();
  ensembleRefreshSelectedConsensus();
};
