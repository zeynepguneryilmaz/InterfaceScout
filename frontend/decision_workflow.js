/* InterfaceScout research workflow extensions. No ANM, GNM, or network analysis. */
/* User-chosen visual exploration presets, not automated chemistry predictions.
   All eleven maps are always calculated and exported. */
const chemistryPresets=[
  {id:"individual",label:"All 11 maps",maps:[],
   rationale:"All maps are available individually, without an imposed overlap."},
  {id:"hydrophobic_aromatic",label:"Hydrophobic + aromatic",maps:["hydrophobic","pi_carbon"],
   rationale:"Explore nearby hydrophobic and aromatic-compatible surface regions."},
  {id:"cationic_donor",label:"Positive + donor",maps:["cationic","hbond_donor"],
   rationale:"Explore neighborhoods with positive-compatible and hydrogen-bond donor features."},
  {id:"anionic_acceptor",label:"Negative + acceptor",maps:["anionic","hbond_acceptor"],
   rationale:"Explore neighborhoods with negative-compatible and hydrogen-bond acceptor features."}
];
let decisionGoal="individual",decisionCombo=null,decisionCompare=null;
let decisionSecondPdb=null,decisionSecondName="",decisionComboError="",decisionRequestId=0;
function decisionInstall(){
  const style=document.createElement("style");
  style.textContent=".workPanel{background:#fff;border:1px solid #dfe4ec;border-radius:8px;margin:10px 0;padding:12px}.workPanel h3{font-size:14px;margin:3px 0 9px}.workPanel p,.workNote{font-size:12px;color:#667085;line-height:1.48}.workMetric{display:inline-block;margin:0 17px 8px 0;font-size:12px}.workAction{font-size:11px;border:1px solid #dfe4ec;border-radius:6px;background:#fff;padding:5px 7px;cursor:pointer;color:#2563eb}.workError{color:#b42318;font-size:12px}";
  document.head.appendChild(style);
  $("run").insertAdjacentHTML("beforebegin",
    '<div class="section"><h3>Quick views</h3>'+
    '<div id="presetChips" class="preset-chips" role="group" aria-label="Quick chemistry views"></div>'+
    '<p id="goalExplanation" class="workspace-quiet" aria-live="polite"></p>'+
    '<p class="workspace-quiet">One click changes the view. All 11 maps remain in the download. Multi-map views show nearby chemistry, not predicted binding.</p></div>');
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
  $("presetChips").innerHTML=chemistryPresets.map(p=>
    '<button type="button" class="preset-chip" data-preset="'+esc(p.id)+
    '" aria-pressed="false" title="'+esc(p.rationale)+
    '" onclick="decisionChooseGoal(this.dataset.preset)">'+esc(p.label)+'</button>').join("");
  decisionChooseGoal("individual");
}
function decisionChooseGoal(id){
  const p=chemistryPresets.find(x=>x.id===id)||chemistryPresets[0];
  decisionGoal=p.id;
  $("goalExplanation").textContent=p.rationale;
  document.querySelectorAll("#presetChips .preset-chip").forEach(button=>{
    button.setAttribute("aria-pressed",String(button.dataset.preset===p.id));
  });
  // A quick view is applied immediately: no dropdown + Apply step.
  if(data)decisionApplyGoal();
}
function decisionApplyGoal(){
  if(!data)return;
  const p=chemistryPresets.find(x=>x.id===decisionGoal)||chemistryPresets[0];
  const keys=p.maps.filter(k=>Object.hasOwn(data.maps,k));
  decisionCombo=null;decisionComboError="";
  if(keys.length>=2){
    selected=new Set(keys);
    activeType="composite";activeMap=null;
  }else{
    selected=new Set([data.initial_map]);
    activeType="chemistry";activeMap=data.initial_map;
  }
  $("error").textContent="";
  render();
  if(keys.length>=2)decisionRefreshCombo();
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
  const donorAcceptorNote=selected.has("hbond_donor")&&selected.has("hbond_acceptor")
    ? '<p class="workNote"><b>Interpretation note:</b> H-bond donor and acceptor are complementary interactions. Their exact same-residue intersection is not itself a binding criterion; some residues can have both capabilities. Examine the channels separately or inspect distinct residues in one neighborhood.</p>'
    : "";
  if(decisionComboError)return '<div class="workPanel"><p class="workError">'+esc(decisionComboError)+'</p></div>';
  if(!decisionCombo)return '<div class="workPanel"><p>Calculating same-residue intersection and spatial co-localization…</p></div>';
  const c=decisionCombo,all=c.regions.filter(function(x){return x.all_maps_colocated});
  const body=c.regions.slice(0,20).map(function(r,i){
    let tag=r.all_maps_colocated?(r.distinct_residue_colocation_only?"Different residues, same patch":"Includes exact intersection"):"Partial overlap";
    return '<tr><td>'+esc(r.source_center_key)+'</td><td>'+r.n_represented_maps+'/'+selected.size+'</td>'+
      '<td>'+esc(tag)+'</td><td><button class="workAction" onclick="decisionHighlightRegion('+i+')">Highlight</button></td></tr>';
  }).join("");
  return '<div class="workPanel"><h3>Combined maps — spatial co-localization</h3>'+donorAcceptorNote+
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
    // Compare mode uses the same full 11-map structure-set analysis as multi-PDB.
    // This enables complete export for BOTH structures, not just the active one.
    const inputs=[
      {name:String(structureName||"first"),pdb_text:pdbText,chain:$("chain").value.trim()||null},
      {name:String(decisionSecondName||id||"second").replace(/\.pdb$/i,""),pdb_text:second,chain:$("secondChain").value.trim()||null},
    ];
    if(inputs[0].name===inputs[1].name)inputs[1].name+="_2";
    const r=await fetch("/structure-set",{method:"POST",headers:{"Content-Type":"application/json"},
      body:JSON.stringify({structures:inputs,...env,combined_maps:selected.size>=2?[...selected]:null})});
    const x=await r.json();if(!r.ok)throw Error(x.detail||"Structure comparison failed");
    if(!x.comparison?.length)throw Error(x.comparison_skipped?.[0]?.reason||"Structures cannot be compared by matched residues");
    decisionCompare=x.comparison[0];
    ensembleResult=x;ensembleAnalyzedInputs=inputs;
    ensembleActiveIndex=0;ensembleConsensusKey=null;
    render();ensembleRefresh();
  }catch(e){$("decisionCompareError").textContent=e.message}
  finally{button.disabled=false}
}
function decisionExportReport(){
  if(!data)return;
  const report={scope:"Descriptive structural decision support, not adsorption prediction",
    chemistry_preset:decisionGoal,inputs:data.input,engine_version:data.version,method:data.method,
    protein_properties:data.protein_properties,selected_maps:Array.from(selected),
    map_combinations:decisionCombo,conformation_comparison:decisionCompare,
    warnings:["No adsorption free energy, unique pose or unfolding prediction",
      "Ionization and covalent reactivity require independent chemical assessment",
      "Ionic strength and temperature are recorded but not explicitly modeled"]};
  downloadText("InterfaceScout_"+safe(structureName)+"_research_report.json",
    JSON.stringify(report,null,2),"application/json");
}
decisionInstall();
