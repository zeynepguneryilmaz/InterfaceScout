/* InterfaceScout research workflow extensions. No ANM, GNM, or network analysis. */
let decisionGoals=[],decisionGoal="explore",decisionCombo=null,decisionCompare=null;
let decisionSecondPdb=null,decisionSecondName="",decisionComboError="",decisionRequestId=0;
function decisionInstall(){
  const style=document.createElement("style");
  style.textContent=".workPanel{background:#fff;border:1px solid #dfe4ec;border-radius:8px;margin:10px 0;padding:12px}.workPanel h3{font-size:14px;margin:3px 0 9px}.workPanel p,.workNote{font-size:12px;color:#667085;line-height:1.48}.workMetric{display:inline-block;margin:0 17px 8px 0;font-size:12px}.workAction{font-size:11px;border:1px solid #dfe4ec;border-radius:6px;background:#fff;padding:5px 7px;cursor:pointer;color:#2563eb}.workError{color:#b42318;font-size:12px}";
  document.head.appendChild(style);
  $("run").insertAdjacentHTML("beforebegin",
    '<div class="section"><h3>Research question</h3><select class="input" id="researchGoal" onchange="decisionChooseGoal(this.value)"><option value="explore">Explore protein surface</option></select>'+
    '<div class="note" id="goalExplanation" style="margin:8px 0">All 11 maps remain available.</div>'+
    '<button class="btn secondary" type="button" onclick="decisionApplyGoal()">Suggest maps</button>'+
    '<div class="workNote">Suggestions do not modify chemistry scores.</div></div>');
  const exportSection=$("excelBtn").closest(".section");
  exportSection.insertAdjacentHTML("afterend",
    '<div class="section"><h3>Compare conformations</h3>'+
    '<div class="field"><label>Second PDB ID (or choose a file)</label><input class="input" id="secondPdbId" placeholder="e.g. 1AKE"></div>'+
    '<div class="field"><label>Second local PDB</label><input class="input" type="file" accept=".pdb" id="secondFile" onchange="decisionReadPdb(this)"></div>'+
    '<div class="field"><label>Second chain (optional)</label><input class="input" id="secondChain" placeholder="e.g. A"></div>'+
    '<button class="btn secondary" id="compareDecisionButton" type="button" onclick="decisionCompareStructures()" disabled>Compare 11 maps</button>'+
    '<div class="workNote">Uses the same pH; compares structures, not dynamic predictions.</div>'+
    '<div id="decisionCompareError" class="workError"></div></div>');
  $("results").insertAdjacentHTML("afterend",
    '<div id="decisionExport" style="padding:0 16px 15px;display:none">'+
    '<button class="btn secondary" onclick="decisionExportReport()">Export research report (JSON)</button></div>');
  fetch("/research-goals").then(function(r){return r.json()}).then(function(goals){
    decisionGoals=goals;
    $("researchGoal").innerHTML=goals.map(function(g){
      return '<option value="'+esc(g.id)+'">'+esc(g.label)+'</option>';
    }).join("");
    decisionChooseGoal("explore");
  }).catch(function(){ $("goalExplanation").textContent="Research-goal descriptions unavailable."; });
}
function decisionChooseGoal(id){
  decisionGoal=id;
  const g=decisionGoals.find(function(x){return x.id===id});
  if(g)$("goalExplanation").textContent=g.question+" "+g.caution;
}
function decisionApplyGoal(){
  if(!data){$("error").textContent="Analyze your protein first.";return}
  const g=decisionGoals.find(function(x){return x.id===decisionGoal});
  if(!g)return;
  const keys=(g.recommended_maps||[]).filter(function(k){return Boolean(data.maps[k])});
  decisionCombo=null;
  if(keys.length){
    selected=new Set(keys);
    activeType=keys.length>1?"composite":"chemistry";
    if(keys.length===1)activeMap=keys[0];
  }else{
    activeType="chemistry";activeMap=data.initial_map;
  }
  $("error").textContent=keys.length?"":"This goal does not require preset maps. Choose any relevant channels.";
  render();
  if(selected.size>1)decisionRefreshCombo();
}
function decisionDirectProperties(){
  if(!data||!data.protein_properties)return "";
  const s=data.protein_properties;
  const labels={nonpolar:"Nonpolar",polar_uncharged:"Polar uncharged",
    acidic_sidechain:"Acidic side chains",basic_sidechain:"Basic side chains",
    aromatic:"Aromatic",sulfur_containing:"Sulfur-containing",cysteine:"Cysteine",
    lysine:"Lysine",carboxylate_bearing:"Carboxylate-bearing"};
  const rows=Object.entries(s.groups).map(function(entry){
    const k=entry[0],g=entry[1];
    return '<tr><td>'+esc(labels[k]||k)+'</td><td>'+g.n_exposed_residues+'</td><td>'+
      esc((g.exposed_residues||[]).join(", "))+'</td></tr>';
  }).join("");
  return '<div class="workPanel"><h3>Direct protein surface properties</h3>'+
    '<p>Protein residue identity and side-chain solvent accessibility, distinct from material-compatibility channels. Group membership can overlap.</p>'+
    '<span class="workMetric"><b>'+s.n_residues+'</b> residues</span>'+
    '<span class="workMetric"><b>'+s.n_surface_residues+'</b> exposed</span>'+
    '<span class="workMetric"><b>'+Number(s.total_sidechain_sasa_A2).toFixed(1)+' Å²</b> side-chain SASA</span>'+
    '<details><summary>Exposed residue groups</summary><table><thead><tr><th>Protein property</th><th>Count</th><th>Residues</th></tr></thead><tbody>'+
    rows+'</tbody></table></details></div>';
}
function decisionFormat(x){return x===null||x===undefined?"—":Number(x).toFixed(3)}
function decisionComboPanel(){
  if(selected.size<2)return "";
  if(decisionComboError)return '<div class="workPanel"><p class="workError">'+esc(decisionComboError)+'</p></div>';
  if(!decisionCombo)return '<div class="workPanel"><p>Calculating same-residue intersection and spatial co-localization…</p></div>';
  const c=decisionCombo,all=c.regions.filter(function(x){return x.all_maps_colocated});
  const body=c.regions.slice(0,20).map(function(r,i){
    let tag=r.all_maps_colocated?(r.distinct_residue_colocation_only?"Different residues, same patch":"Includes exact intersection"):"Partial overlap";
    return '<tr><td>'+esc(r.source_center_key)+'</td><td>'+r.n_represented_maps+'/'+selected.size+'</td>'+
      '<td>'+esc(tag)+'</td><td><button class="workAction" onclick="decisionHighlightRegion('+i+')">Highlight</button></td></tr>';
  }).join("");
  return '<div class="workPanel"><h3>Combined maps — spatial co-localization</h3>'+
    '<p><b>Exact intersection</b>: same residue supports every selected channel. <b>Co-localization</b>: different residues supporting selected channels within one existing 8 Å patch. Not a combined binding score.</p>'+
    '<span class="workMetric"><b>'+c.n_exact_overlap+'</b> exactly shared residues</span>'+
    '<span class="workMetric"><b>'+all.length+'</b> fully co-localized candidate patches (shown)</span>'+
    '<details><summary>Exact shared residue list</summary><p>'+esc(c.exact_same_residue_overlap.join(", ")||"None")+'</p></details>'+
    '<table><thead><tr><th>Center</th><th>Map coverage</th><th>Type</th><th>Inspect</th></tr></thead><tbody>'+body+'</tbody></table>'+
    '<p>Candidate neighborhoods use the original single-channel patch geometry; the list does not constitute a separately validated ranking.</p></div>';
}
function decisionComparePanel(){
  if(!decisionCompare)return "";
  const c=decisionCompare,m=c.mapping,s=c.accessibility;
  const rows=c.map_comparisons.map(function(x){
    return '<tr><td>'+esc(x.label)+'</td><td>'+esc(x.top1_center_a||"—")+' → '+esc(x.top1_center_b||"—")+
      '</td><td>'+decisionFormat(x.top_k[0].jaccard)+'</td><td>'+decisionFormat(x.top_k[1].jaccard)+'</td><td>'+
      decisionFormat(x.top_k[2].jaccard)+'</td></tr>';
  }).join("");
  return '<div class="workPanel"><h3>Conformation comparison: '+esc(decisionSecondName||"second structure")+'</h3>'+
    '<p>Same protein, matched residue identities and identical solution pH. A comparison of candidate regions, not an MD or equilibrium ensemble.</p>'+
    '<span class="workMetric"><b>'+m.n_identity_matches+'</b> matched residues</span>'+
    '<span class="workMetric"><b>'+s.n_exposure_status_changes+'</b> exposure-class changes</span>'+
    '<span class="workMetric"><b>'+s.n_exposed_matched_a+' → '+s.n_exposed_matched_b+'</b> exposed residues</span>'+
    '<table><thead><tr><th>Channel</th><th>Top-1 centers</th><th>Top-1 J</th><th>Top-3 J</th><th>Top-5 J</th></tr></thead><tbody>'+
    rows+'</tbody></table><details><summary>Changed exposure residue IDs</summary><p>'+
    esc(c.residue_changes.filter(function(x){return x.changed_exposure_class}).map(function(x){return x.key}).join(", ")||"None")+
    '</p></details></div>';
}
const decisionOriginalRender=render;
render=function(){
  decisionOriginalRender();
  if(data){
    $("results").insertAdjacentHTML("beforeend",decisionDirectProperties()+decisionComboPanel()+decisionComparePanel());
    $("decisionExport").style.display="block";
    $("compareDecisionButton").disabled=false;
  }
};
const decisionOriginalAnalyze=analyze;
analyze=async function(){
  decisionCombo=null;decisionCompare=null;decisionComboError="";
  await decisionOriginalAnalyze();
  if(data){render();if(selected.size>1)decisionRefreshCombo()}
};
const decisionOriginalToggle=toggleChem;
toggleChem=function(k){
  decisionCombo=null;decisionComboError="";
  decisionOriginalToggle(k);
  if(selected.size>1)decisionRefreshCombo();
};
const decisionOriginalComposite=openComposite;
openComposite=function(){
  decisionOriginalComposite();
  if(data&&selected.size>1&&!decisionCombo)decisionRefreshCombo();
};
async function decisionRefreshCombo(){
  if(!data||selected.size<2)return;
  const id=++decisionRequestId;
  try{
    const r=await fetch("/combine",{method:"POST",headers:{"Content-Type":"application/json"},
      body:JSON.stringify({analysis:data,selected_maps:Array.from(selected),max_regions:50})});
    const x=await r.json();if(!r.ok)throw Error(x.detail||"Combined maps failed");
    if(id!==decisionRequestId)return;
    decisionCombo=x;decisionComboError="";render();
  }catch(e){if(id!==decisionRequestId)return;decisionComboError=e.message;render()}
}
function decisionHighlightRegion(i){
  if(!decisionCombo||!viewer)return;
  const r=decisionCombo.regions[i];if(!r)return;
  rebuildViewer(true);
  const sel=r.member_keys.map(function(key){
    const parts=key.split(":");return {chain:parts[0],resi:Number(parts[1])};
  });
  sel.forEach(function(s){viewer.addStyle(s,{stick:{color:"#ea580c",radius:0.25},sphere:{color:"#ea580c",radius:0.44}})});
  try{viewer.zoomTo({or:sel},400)}catch(e){}
  viewer.render();
}
function decisionReadPdb(input){
  const file=input.files&&input.files[0];decisionSecondPdb=null;decisionSecondName="";
  if(!file)return;
  const reader=new FileReader();
  reader.onload=function(){decisionSecondPdb=reader.result;decisionSecondName=file.name;$("decisionCompareError").textContent=""};
  reader.onerror=function(){$("decisionCompareError").textContent="Unable to read PDB file"};
  reader.readAsText(file);
}
async function decisionCompareStructures(){
  if(!data||!pdbText)return;
  const id=$("secondPdbId").value.trim().toUpperCase();
  const button=$("compareDecisionButton");
  $("decisionCompareError").textContent="";button.disabled=true;
  try{
    let second=decisionSecondPdb;
    if(!second){
      if(!/^[0-9][A-Z0-9]{3}$/.test(id))throw Error("Supply a four-character PDB ID or a local PDB");
      const r=await fetch("https://files.rcsb.org/download/"+encodeURIComponent(id)+".pdb");
      if(!r.ok)throw Error("Second structure could not be fetched");
      second=await r.text();decisionSecondName=id;
    }
    const env={pH:Number($("pH").value),ionic_mM:Number($("ionic").value),temp_K:Number($("temp").value)};
    const request={structure_a:Object.assign({pdb_text:pdbText,chain:$("chain").value.trim()||null},env),
      structure_b:Object.assign({pdb_text:second,chain:$("secondChain").value.trim()||null},env)};
    const r=await fetch("/compare",{method:"POST",headers:{"Content-Type":"application/json"},body:JSON.stringify(request)});
    const x=await r.json();if(!r.ok)throw Error(x.detail||"Structure comparison failed");
    decisionCompare=x;render();
  }catch(e){$("decisionCompareError").textContent=e.message}
  finally{button.disabled=false}
}
function decisionExportReport(){
  if(!data)return;
  const report={scope:"Descriptive structural decision support, not adsorption prediction",
    research_goal:decisionGoal,inputs:data.input,engine_version:data.version,method:data.method,
    protein_properties:data.protein_properties,selected_maps:Array.from(selected),
    map_combinations:decisionCombo,conformation_comparison:decisionCompare,
    warnings:["No adsorption free energy, unique pose or unfolding prediction",
      "Ionization and covalent reactivity require independent chemical assessment",
      "Ionic strength and temperature are recorded but not explicitly modeled"]};
  downloadText("InterfaceScout_"+safe(structureName)+"_research_report.json",
    JSON.stringify(report,null,2),"application/json");
}
decisionInstall();
