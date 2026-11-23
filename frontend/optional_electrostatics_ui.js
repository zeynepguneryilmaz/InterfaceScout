/* Optional external electrostatics: does not change canonical 11-map scoring. */
(function(){
 const side=document.querySelector("aside");if(!side)return;
 const section=document.createElement("section");section.className="section";
 section.style.cssText="padding:12px;border:1px solid #dfe4ec;border-radius:10px;background:#fff;margin-top:14px";
 section.innerHTML=`
 <h3>Optional electrostatics · PROPKA / APBS</h3>
 <p class="workspace-quiet">Structure-dependent protonation and continuum electrostatics.
  The existing 11 chemistry maps and their ranks remain unchanged.</p>
 <div class="field"><label>Electrolyte</label><select class="input" id="esSalt">
  <option value="NaCl">NaCl (1:1)</option><option value="KCl">KCl (1:1)</option></select></div>
 <div class="row"><div class="field"><label>Forcefield</label>
  <select class="input" id="esFF"><option>PARSE</option><option>AMBER</option><option>CHARMM</option></select></div>
 <div class="field"><label>Protein dielectric</label><input class="input" id="esPD" type="number" value="2.0" min="1" max="30" step="0.5"></div></div>
 <div class="field"><label>Solvent dielectric</label><input class="input" id="esSD" type="number" value="78.5" min="20" max="110" step="0.5"></div>
 <label style="display:block;margin:6px 0;font-size:12px"><input id="esPropka" type="checkbox" checked> PROPKA protonation at the selected pH</label>
 <label style="display:block;margin:6px 0;font-size:12px"><input id="esAPBS" type="checkbox" checked> APBS electrostatic potential grid</label>
 <button type="button" class="btn secondary" id="esRun" disabled>Calculate &amp; download electrostatics ZIP</button>
 <p id="esMsg" class="workspace-quiet" role="status">Checking optional dependencies…</p>
 <p class="workspace-quiet">pH, ionic strength (mM) and temperature (K) are read from the solution fields above.
 Protein-only input: removed ligands, glycans, ions and other HETATM species are flagged in the report.
 Calculations use a linearized Poisson–Boltzmann model and are not adsorption predictions.</p>`;
 side.appendChild(section);
 const byId=id=>document.getElementById(id);
 let ready=false,busy=false;
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
 btn.addEventListener("click",async()=>{
   if(!pdbText||busy)return;
   busy=true;refresh();msg.textContent="Running external electrostatics (may take time)…";
   try{
     const num=id=>Number(byId(id).value);
     const payload={
       pdb_text:pdbText,chain:byId("chain").value.trim()||null,
       pH:num("pH"),ionic_mM:num("ionic"),temperature_K:num("temp"),
       electrolyte:byId("esSalt").value,forcefield:byId("esFF").value,
       protein_dielectric:num("esPD"),solvent_dielectric:num("esSD"),
       propka_enabled:byId("esPropka").checked,apbs_enabled:byId("esAPBS").checked
     };
     const r=await fetch("/electrostatics/export",{method:"POST",headers:{"Content-Type":"application/json"},body:JSON.stringify(payload)});
     if(!r.ok){const e=await r.json().catch(()=>({}));throw new Error(e.detail||("HTTP "+r.status));}
     const blob=await r.blob();
     const href=URL.createObjectURL(blob),a=document.createElement("a");
     a.href=href;a.download="InterfaceScout_Electrostatics_Optional.zip";a.click();
     setTimeout(()=>URL.revokeObjectURL(href),3000);
     msg.textContent="Electrostatics evidence ZIP generated. Read metadata/results.json for limitations and provenance.";
   }catch(e){msg.textContent="Electrostatics not completed: "+String(e.message||e);}
   finally{busy=false;refresh();}
 });
 setInterval(refresh,1500);
 check();
})();