/* Optional external electrostatics: does not change canonical 11-map scoring. */
(function(){
 const side=document.querySelector("aside");if(!side)return;
 const section=document.createElement("section");section.className="section";
 section.style.cssText="padding:12px;border:1px solid #dfe4ec;border-radius:10px;background:#fff;margin-top:14px";
 section.innerHTML=`
 <h3>Optional electrostatics · PROPKA / APBS</h3>
 <p class="workspace-quiet">Structure-dependent protonation and continuum electrostatics.
  Calculate an experimental electrostatics-informed 11-map view without changing the standard model.</p>
 <div class="field"><label>Electrolyte</label><select class="input" id="esSalt">
  <option value="NaCl">NaCl (1:1)</option><option value="KCl">KCl (1:1)</option></select></div>
 <div class="row"><div class="field"><label>Forcefield</label>
  <select class="input" id="esFF"><option>PARSE</option><option>AMBER</option><option>CHARMM</option></select></div>
 <div class="field"><label>Protein dielectric</label><input class="input" id="esPD" type="number" value="2.0" min="1" max="30" step="0.5"></div></div>
 <div class="field"><label>Solvent dielectric</label><input class="input" id="esSD" type="number" value="78.5" min="20" max="110" step="0.5"></div>
 <label style="display:block;margin:6px 0;font-size:12px"><input id="esPropka" type="checkbox" checked> PROPKA protonation at the selected pH</label>
 <label style="display:block;margin:6px 0;font-size:12px"><input id="esAPBS" type="checkbox" checked> APBS electrostatic potential grid</label>
 <button type="button" class="btn secondary" id="esRun" disabled>Calculate electrostatics + updated maps</button>
 <div id="esOutputs" style="display:none;margin-top:10px">
 <div class="row"><button type="button" class="btn secondary" id="esShowPhysical">Show updated 11 maps</button>
 <button type="button" class="btn secondary" id="esShowBaseline">Show original 11 maps</button></div>
 <a id="esDownload" class="btn secondary" style="display:block;text-align:center;text-decoration:none;margin-top:7px"
 download="InterfaceScout_PhysicsAware_Electrostatics.zip">Download complete electrostatics + 11-map ZIP</a>
 </div>
 <p id="esMsg" class="workspace-quiet" role="status">Checking optional dependencies…</p>
 <p class="workspace-quiet">pH, ionic strength (mM) and temperature (K) are read from the solution fields above.
 Protein-only input: removed ligands, glycans, ions and other HETATM species are flagged in the report.
 Calculations use a linearized Poisson–Boltzmann model and are not adsorption predictions.</p>`;
 side.appendChild(section);
 const byId=id=>document.getElementById(id);
 let ready=false,busy=false,latest=null;
 const btn=byId("esRun"),msg=byId("esMsg");
 function refresh(){btn.disabled=!ready||busy||typeof pdbText==='undefined'||!pdbText;}
 async function check(){
   try{
      const r=await fetch("/electrostatics/status");
      if(!r.ok)throw new Error("Status service unavailable");
      const s=await r.json();
      ready=Boolean(s.pdb2pqr);
      msg.textContent=!s.pdb2pqr?"PDB2PQR not installed. Standard analysis remains available.":
         (!s.apbs?"PDB2PQR available; APBS missing. Disable APBS for PQR-only output.":"Optional PDB2PQR and APBS are available.");
   }catch(e){ready=false;msg.textContent="Optional electrostatics service unavailable: "+e.message;}
   refresh();
 }
 function showModel(analysis){
   if(!latest||!analysis||pdbText!==latest.pdb||
         (byId("chain").value.trim()||null)!==latest.chain){
      msg.textContent="Loaded PDB/chain changed. Recalculate electrostatics before viewing results.";
      return;
   }
   data=analysis;
   activeType="chemistry";
   activeMap=analysis.initial_map||analysis.map_order[0];
   selected=new Set([activeMap]);
   focusedPatchRank=null;
   byId("excelBtn").disabled=false;
   byId("pdbBtn").disabled=false;
   render();
   const experimental=analysis===latest.physical;
   msg.textContent=experimental?
       "Showing electrostatics-informed 11 maps (experimental model). These are chemical support regions, not adsorption predictions.":
       "Showing original validated 11-map model. Physics-aware results are preserved separately.";
 }
 byId("esShowPhysical").addEventListener("click",()=>showModel(latest?.physical));
 byId("esShowBaseline").addEventListener("click",()=>showModel(latest?.source));
 btn.addEventListener("click",async()=>{
   if(!pdbText||busy)return;
   busy=true;refresh();msg.textContent="Running PDB2PQR/PROPKA and APBS; preparing updated maps…";
   byId("esOutputs").style.display="none";
   try{
     const num=id=>Number(byId(id).value);
     const payload={
       pdb_text:pdbText,chain:byId("chain").value.trim()||null,
       pH:num("pH"),ionic_mM:num("ionic"),temperature_K:num("temp"),
       electrolyte:byId("esSalt").value,forcefield:byId("esFF").value,
       protein_dielectric:num("esPD"),solvent_dielectric:num("esSD"),
       propka_enabled:byId("esPropka").checked,apbs_enabled:byId("esAPBS").checked
     };
     const snap=pdbText;
     const r=await fetch("/electrostatics/prepare",{method:"POST",headers:{"Content-Type":"application/json"},body:JSON.stringify(payload)});
     if(!r.ok){const e=await r.json().catch(()=>({}));throw new Error(e.detail||("HTTP "+r.status));}
     const result=await r.json();
     if(pdbText!==snap){throw new Error("PDB changed during calculation. Rerun electrostatics.")}
     latest={pdb:snap,chain:payload.chain,source:result.canonical_analysis,
             physical:result.analysis,result:result.result,url:result.download_url};
     if(!result.analysis){
        msg.textContent="PQR charge preparation complete. Enable APBS to generate physical 11-map view.";
     }else{
        msg.textContent="Completed: 11 electrostatics-informed maps available. Compared against unchanged standard model.";
     }
     byId("esDownload").href=result.download_url;
     byId("esOutputs").style.display="block";
     byId("esShowPhysical").disabled=!result.analysis;
     byId("esShowBaseline").disabled=!result.canonical_analysis;
     if(result.analysis)showModel(result.analysis);
   }catch(e){msg.textContent="Electrostatics not completed: "+String(e.message||e);}
   finally{busy=false;refresh();}
 });
 setInterval(refresh,1500);
 check();
})();