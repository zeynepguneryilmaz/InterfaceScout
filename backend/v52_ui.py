"""Small, deterministic UI augmentation for the v5.2 candidate.

The publication frontend is kept intact on disk.  The v5.2 app injects only the
controls/output needed for structural context, auxiliary CX protrusion, material
profiles, and per-run applicability notes.
"""

from __future__ import annotations


def inject_v52_ui(html: str) -> str:
    # Add compact v5.2 controls beneath the chain selector.
    needle = '<div class="field"><label>CHAIN (optional; blank = all chains)</label><input class="input" id="chain" placeholder="e.g. A"></div>'
    controls = needle + '''
      <div class="field"><label>STRUCTURAL CONTEXT</label>
        <select class="input" id="structureContext">
          <option value="auto" selected>Auto · RCSB biological assembly 1 when available</option>
          <option value="biological_assembly_1">Biological assembly 1 · require RCSB assembly</option>
          <option value="deposited_structure">Deposited / uploaded coordinates</option>
          <option value="selected_chain_legacy">Selected chain only · legacy / diagnostic</option>
        </select>
      </div>
      <div class="field"><label>MATERIAL PROFILE (optional)</label>
        <select class="input" id="materialProfile"><option value="">None · choose chemistry maps manually</option></select>
      </div>
      <div class="field" style="display:flex;align-items:center;gap:7px">
        <input id="protrusion" type="checkbox" checked style="margin:0">
        <label for="protrusion" style="margin:0;font-family:inherit;font-size:11px">Report CX protrusion (auxiliary; does not change scores)</label>
      </div>'''
    if needle in html:
        html = html.replace(needle, controls, 1)

    # Add an applicability/result-scope box below the results grid.
    results_close = '    </section>\n  </main>'
    scope_box = '''    </section>
    <section id="applicabilityPanel" style="display:none;background:#fff;border-top:1px solid var(--border);padding:10px 14px;max-height:155px;overflow:auto">
      <div class="resultTitle">Applicability of this run</div>
      <div id="applicabilityBody" class="footerline"></div>
    </section>
  </main>'''
    if results_close in html:
        html = html.replace(results_close, scope_box, 1)

    # Candidate wording in theory modal without changing the frozen equations.
    html = html.replace('InterfaceScout publication-freeze model', 'InterfaceScout v5.2 structural candidate · frozen v5.1 scoring core')
    html = html.replace(
        '<p><strong>Not predicted:</strong> adsorption capacity, absolute free energy, unique adsorption orientation, adsorption-induced conformational change, or material-side porosity/transport/hydration.</p>',
        '<p><strong>v5.2 structural refinements:</strong> structural context may use RCSB biological assembly 1; Pintar CX protrusion is reported only as an auxiliary geometry descriptor; material profiles expose predeclared chemistry channels separately without a weighted combination.</p><p><strong>Not predicted:</strong> adsorption capacity, absolute free energy, unique adsorption orientation, adsorption-induced conformational change, explicit hydration/desolvation, or multi-protein corona organization.</p>'
    )

    # Inject v5.2 Javascript immediately after the existing state declaration.
    js_needle = "let pdbText=null, viewer=null, data=null, activeChem=null, colorMode='propensity', surfaces=[];"
    js = js_needle + r'''
let loadedPdbId=null;

async function loadMaterialProfiles(){
  const sel=document.getElementById('materialProfile'); if(!sel)return;
  try{
    const r=await fetch(API()+'/material_profiles'); if(!r.ok)return;
    const d=await r.json();
    (d.profiles||[]).forEach(p=>{
      const o=document.createElement('option');o.value=p.key;o.textContent=p.label;sel.appendChild(o);
    });
  }catch(e){}
}

function renderApplicability(){
  const panel=document.getElementById('applicabilityPanel');
  const body=document.getElementById('applicabilityBody');
  if(!panel||!body||!data||!data.applicability){return;}
  panel.style.display='block';
  const inc=data.applicability.included_in_this_run||[];
  const lim=data.applicability.not_included_or_interpretation_limits||[];
  const context=(data.settings&&data.settings.structure_context)||'—';
  const profile=data.material_profile?`${data.material_profile.label} · ${(data.material_profile.channels||[]).join(', ')}`:'none';
  body.innerHTML=`<strong>Structural context:</strong> ${escapeHtml(context)} &nbsp; · &nbsp; <strong>Material profile:</strong> ${escapeHtml(profile)}<br>`+
    `<strong>Included:</strong> ${inc.map(escapeHtml).join(' · ')}<br>`+
    `<strong>Interpretation limits:</strong> ${lim.map(escapeHtml).join(' · ')}`;
}
'''
    if js_needle in html:
        html = html.replace(js_needle, js, 1)

    # Replace fetchPDB/loadFile so the PDB ID is preserved for assembly lookup.
    old_fetch = "async function fetchPDB(){const id=document.getElementById('pdbId').value.trim().toUpperCase();if(!id)return;showLoading(`${id} downloading…`);showError('');try{const r=await fetch(`https://files.rcsb.org/download/${id}.pdb`);if(!r.ok)throw new Error('PDB not found');pdbText=await r.text();loadViewer(pdbText);document.getElementById('proteinInfo').innerHTML=`<strong>${id}</strong> loaded from RCSB.`;}catch(e){showError(e.message)}finally{hideLoading()}}"
    new_fetch = "async function fetchPDB(){const id=document.getElementById('pdbId').value.trim().toUpperCase();if(!id)return;showLoading(`${id} downloading…`);showError('');try{const r=await fetch(`https://files.rcsb.org/download/${id}.pdb`);if(!r.ok)throw new Error('PDB not found');pdbText=await r.text();loadedPdbId=id;loadViewer(pdbText);document.getElementById('proteinInfo').innerHTML=`<strong>${id}</strong> loaded from RCSB.`;}catch(e){showError(e.message)}finally{hideLoading()}}"
    html = html.replace(old_fetch, new_fetch, 1)
    old_file = "function loadFile(ev){const f=ev.target.files[0];if(!f)return;const rd=new FileReader();rd.onload=e=>{pdbText=e.target.result;loadViewer(pdbText);document.getElementById('proteinInfo').innerHTML=`<strong>${escapeHtml(f.name)}</strong> loaded.`;};rd.readAsText(f);}"
    new_file = "function loadFile(ev){const f=ev.target.files[0];if(!f)return;const rd=new FileReader();rd.onload=e=>{pdbText=e.target.result;loadedPdbId=null;loadViewer(pdbText);document.getElementById('proteinInfo').innerHTML=`<strong>${escapeHtml(f.name)}</strong> loaded.`;};rd.readAsText(f);}"
    html = html.replace(old_file, new_file, 1)

    # Replace the analysis payload/result summary, preserving all legacy display logic.
    old_analyze = "async function analyze(){if(!pdbText)return;showLoading('scRSA → chemistry mapping → 5/8 Å persistence → auxiliary descriptors…');showError('');try{const r=await fetch(API()+'/analyze_surface',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({pdb_text:pdbText,chain:document.getElementById('chain').value.trim()||null,env:env()}),signal:AbortSignal.timeout(600000)});const d=await r.json();if(!r.ok)throw new Error(d.detail||'Backend error');data=d;document.getElementById('proteinInfo').innerHTML=`<strong>${d.settings.chain==='ALL'?'All chains':'Chain '+escapeHtml(d.settings.chain)}</strong> · ${d.stats.n_residues} residues · ${d.stats.n_surface_res} surface residues<br>Electrostatics: ${escapeHtml(d.stats.electrostatics)}`;document.getElementById('csvBtn').disabled=false;document.getElementById('pdbBtn').disabled=false;activeChem=d.chemistry_list[0];colorMode='propensity';renderAll();if(viewer){viewer.zoomTo(analysisSelector(),500);viewer.render();}}catch(e){showError(e.message)}finally{hideLoading()}}"
    new_analyze = "async function analyze(){if(!pdbText)return;showLoading('structure context → scRSA → chemistry mapping → 5/8 Å persistence → auxiliary descriptors…');showError('');try{const payload={pdb_text:pdbText,pdb_id:loadedPdbId,chain:document.getElementById('chain').value.trim()||null,env:env(),structure_context:document.getElementById('structureContext').value,protrusion:document.getElementById('protrusion').checked,material_profile:document.getElementById('materialProfile').value||null};const r=await fetch(API()+'/analyze_surface',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(payload),signal:AbortSignal.timeout(600000)});const d=await r.json();if(!r.ok)throw new Error(d.detail||'Backend error');data=d;document.getElementById('proteinInfo').innerHTML=`<strong>${d.settings.chain==='ALL'?'All chains':'Chain '+escapeHtml(d.settings.chain)}</strong> · ${d.stats.n_residues} residues · ${d.stats.n_surface_res} surface residues<br>Context: ${escapeHtml(d.settings.structure_context||'—')} · Electrostatics: ${escapeHtml(d.stats.electrostatics)}`;document.getElementById('csvBtn').disabled=false;document.getElementById('pdbBtn').disabled=false;activeChem=(d.material_profile&&d.material_profile.channels&&d.material_profile.channels[0])||d.chemistry_list[0];colorMode='propensity';renderAll();renderApplicability();if(viewer){viewer.zoomTo(analysisSelector(),500);viewer.render();}}catch(e){showError(e.message)}finally{hideLoading()}}"
    html = html.replace(old_analyze, new_analyze, 1)

    # Add CX/context/profile columns to CSV exports. These are descriptors only.
    html = html.replace(
        "'analysis_version','selected_chain','pH','ionic_strength_mM','temperature_K',",
        "'analysis_version','selected_chain','structure_context','material_profile','pH','ionic_strength_mM','temperature_K',"
    )
    html = html.replace(
        "'charge_fraction','charge_descriptor','reference_pKa','ionization_sensitive',",
        "'charge_fraction','charge_descriptor','reference_pKa','ionization_sensitive','cx_residue_mean','cx_sidechain_mean','cx_max','cx_ca',"
    )
    html = html.replace(
        "data.settings.chain,\n      data.settings.pH,",
        "data.settings.chain,\n      data.settings.structure_context||'',\n      data.settings.material_profile||'',\n      data.settings.pH,"
    )
    html = html.replace(
        "r.charge_fraction,r.charge_descriptor,r.pka??'',r.ionization_sensitive,\n      r.ss",
        "r.charge_fraction,r.charge_descriptor,r.pka??'',r.ionization_sensitive,\n      r.cx_residue_mean??'',r.cx_sidechain_mean??'',r.cx_max??'',r.cx_ca??'',\n      r.ss"
    )

    # Populate profile selector after existing initialization.
    html = html.replace('initViewer();checkBackend();setInterval(checkBackend,10000);', 'initViewer();checkBackend();loadMaterialProfiles();setInterval(checkBackend,10000);', 1)
    return html
