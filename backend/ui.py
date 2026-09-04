"""UI augmentation for InterfaceScout 2.0 protein-centered development.

The UI exposes structural context and optional protected functional residues,
then reports a protein-derived target interface profile. No material library or
named-material selector is shown.
"""
from __future__ import annotations


def inject_ui(html: str) -> str:
    needle = '<div class="field"><label>CHAIN (optional; blank = all chains)</label><input class="input" id="chain" placeholder="e.g. A"></div>'
    controls = needle + '''
      <div class="field"><label>STRUCTURAL CONTEXT</label>
        <select class="input" id="structureContext">
          <option value="auto" selected>Auto · biological assembly 1 when available</option>
          <option value="biological_assembly_1">Biological assembly 1 · require RCSB assembly</option>
          <option value="deposited_structure">Deposited / uploaded coordinates</option>
          <option value="selected_chain_legacy">Selected chain only · v1 regression mode</option>
        </select>
      </div>
      <div class="field"><label>PROTECTED FUNCTIONAL RESIDUES (optional)</label>
        <input class="input" id="protectedResidues" placeholder="e.g. A:42, A:57, B:103">
        <div style="font-size:10px;opacity:.72;margin-top:3px">Use known catalytic, ligand-binding, epitope or other function-critical residues. These annotations do not change compatibility scores.</div>
      </div>
      <div class="field" style="display:flex;align-items:center;gap:7px">
        <input id="protrusion" type="checkbox" checked style="margin:0">
        <label for="protrusion" style="margin:0;font-family:inherit;font-size:11px">Report CX protrusion (auxiliary)</label>
      </div>'''
    if needle in html:
        html = html.replace(needle, controls, 1)

    results_close = '    </section>\n  </main>'
    extra = '''    </section>
    <section id="targetProfilePanel" style="display:none;background:#fff;border-top:1px solid var(--border);padding:12px 14px;max-height:280px;overflow:auto">
      <div class="resultTitle">Protein-derived target interface profile</div>
      <div id="targetProfileBody" class="footerline"></div>
    </section>
    <section id="applicabilityPanel" style="display:none;background:#fff;border-top:1px solid var(--border);padding:10px 14px;max-height:155px;overflow:auto">
      <div class="resultTitle">Applicability of this run</div>
      <div id="applicabilityBody" class="footerline"></div>
    </section>
  </main>'''
    if results_close in html:
        html = html.replace(results_close, extra, 1)

    html = html.replace('InterfaceScout publication-freeze model', 'InterfaceScout 2.0 development · protein-centered target interface profiling')
    html = html.replace(
        '<p><strong>Not predicted:</strong> adsorption capacity, absolute free energy, unique adsorption orientation, adsorption-induced conformational change, or material-side porosity/transport/hydration.</p>',
        '<p><strong>InterfaceScout 2.0:</strong> derives generalized target interface properties directly from the protein and maps the protein patches that could engage each property. It does not use a named-material library.</p><p><strong>Functional-site caution:</strong> PDB SITE records are generic annotations, not automatically catalytic active sites. Add known function-critical residues in the Protected Functional Residues field when available.</p>'
    )

    js_needle = "let pdbText=null, viewer=null, data=null, activeChem=null, colorMode='propensity', surfaces=[];"
    js = js_needle + r'''
let loadedPdbId=null;

function parseProtectedResidues(){
  const el=document.getElementById('protectedResidues');
  if(!el)return [];
  return el.value.split(',').map(x=>x.trim()).filter(Boolean);
}

function renderTargetProfile(){
  const panel=document.getElementById('targetProfilePanel');
  const body=document.getElementById('targetProfileBody');
  if(!panel||!body||!data)return;
  const p=data.protein_derived_target_interface_profile;
  if(!p){panel.style.display='none';return;}
  panel.style.display='block';
  const channels=(p.interface_channels||[]).slice().sort((a,b)=>
    (b.max_patch_density_8A_raw||0)-(a.max_patch_density_8A_raw||0)
  );
  const rows=channels.map(ch=>{
    const patches=(ch.top_patches||[]).slice(0,3).map(x=>{
      const site=x.functional_site_relation||'—';
      return `${escapeHtml(x.center_key||'—')} · M=${escapeHtml(String(x.multiscale_persistence??'—'))} · ${escapeHtml(site)}`;
    }).join('<br>');
    return `<tr><td style="padding:4px 8px 4px 0"><strong>${escapeHtml(ch.target_surface_property||ch.label||ch.key)}</strong><br><span style="opacity:.75">${escapeHtml(ch.label||ch.key)}</span></td>`+
      `<td style="padding:4px 8px">${escapeHtml(String(ch.n_compatible_surface_residues||0))} residues<br>raw ΣL=${escapeHtml(String(ch.total_accessible_compatibility_raw??0))}</td>`+
      `<td style="padding:4px 0">${patches||'No compatible patch'}</td></tr>`;
  }).join('');
  const siteCount=(p.pdb_site_annotations||[]).length;
  const protectedCount=(p.user_protected_residue_keys||[]).length;
  body.innerHTML=`<strong>Basis:</strong> protein only · no material library · no cross-channel weighted score<br>`+
    `<strong>Functional annotations:</strong> ${siteCount} PDB SITE record(s) · ${protectedCount} user-protected residue(s)<br>`+
    `<table style="width:100%;border-collapse:collapse;margin-top:7px"><thead><tr><th style="text-align:left">Target interface property</th><th style="text-align:left">Protein evidence</th><th style="text-align:left">Top protein patches / function relation</th></tr></thead><tbody>${rows}</tbody></table>`;
}

function renderApplicability(){
  const panel=document.getElementById('applicabilityPanel');
  const body=document.getElementById('applicabilityBody');
  if(!panel||!body||!data||!data.applicability)return;
  panel.style.display='block';
  const inc=data.applicability.included_in_this_run||[];
  const lim=data.applicability.not_included_or_interpretation_limits||[];
  body.innerHTML=`<strong>Included:</strong> ${inc.map(escapeHtml).join(' · ')}<br>`+
    `<strong>Interpretation limits:</strong> ${lim.map(escapeHtml).join(' · ')}`;
}
'''
    if js_needle in html:
        html = html.replace(js_needle, js, 1)

    old_fetch = "async function fetchPDB(){const id=document.getElementById('pdbId').value.trim().toUpperCase();if(!id)return;showLoading(`${id} downloading…`);showError('');try{const r=await fetch(`https://files.rcsb.org/download/${id}.pdb`);if(!r.ok)throw new Error('PDB not found');pdbText=await r.text();loadViewer(pdbText);document.getElementById('proteinInfo').innerHTML=`<strong>${id}</strong> loaded from RCSB.`;}catch(e){showError(e.message)}finally{hideLoading()}}"
    new_fetch = "async function fetchPDB(){const id=document.getElementById('pdbId').value.trim().toUpperCase();if(!id)return;showLoading(`${id} downloading…`);showError('');try{const r=await fetch(`https://files.rcsb.org/download/${id}.pdb`);if(!r.ok)throw new Error('PDB not found');pdbText=await r.text();loadedPdbId=id;loadViewer(pdbText);document.getElementById('proteinInfo').innerHTML=`<strong>${id}</strong> loaded from RCSB.`;}catch(e){showError(e.message)}finally{hideLoading()}}"
    html = html.replace(old_fetch, new_fetch, 1)

    old_file = "function loadFile(ev){const f=ev.target.files[0];if(!f)return;const rd=new FileReader();rd.onload=e=>{pdbText=e.target.result;loadViewer(pdbText);document.getElementById('proteinInfo').innerHTML=`<strong>${escapeHtml(f.name)}</strong> loaded.`;};rd.readAsText(f);}"
    new_file = "function loadFile(ev){const f=ev.target.files[0];if(!f)return;const rd=new FileReader();rd.onload=e=>{pdbText=e.target.result;loadedPdbId=null;loadViewer(pdbText);document.getElementById('proteinInfo').innerHTML=`<strong>${escapeHtml(f.name)}</strong> loaded.`;};rd.readAsText(f);}"
    html = html.replace(old_file, new_file, 1)

    old_analyze = "async function analyze(){if(!pdbText)return;showLoading('scRSA → chemistry mapping → 5/8 Å persistence → auxiliary descriptors…');showError('');try{const r=await fetch(API()+'/analyze_surface',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({pdb_text:pdbText,chain:document.getElementById('chain').value.trim()||null,env:env()}),signal:AbortSignal.timeout(600000)});const d=await r.json();if(!r.ok)throw new Error(d.detail||'Backend error');data=d;document.getElementById('proteinInfo').innerHTML=`<strong>${d.settings.chain==='ALL'?'All chains':'Chain '+escapeHtml(d.settings.chain)}</strong> · ${d.stats.n_residues} residues · ${d.stats.n_surface_res} surface residues<br>Electrostatics: ${escapeHtml(d.stats.electrostatics)}`;document.getElementById('csvBtn').disabled=false;document.getElementById('pdbBtn').disabled=false;activeChem=d.chemistry_list[0];colorMode='propensity';renderAll();if(viewer){viewer.zoomTo(analysisSelector(),500);viewer.render();}}catch(e){showError(e.message)}finally{hideLoading()}}"
    new_analyze = "async function analyze(){if(!pdbText)return;showLoading('protein surface → interaction channels → patch mapping → functional-site relation…');showError('');try{const payload={pdb_text:pdbText,pdb_id:loadedPdbId,chain:document.getElementById('chain').value.trim()||null,env:env(),structure_context:document.getElementById('structureContext').value,protrusion:document.getElementById('protrusion').checked,protected_residue_keys:parseProtectedResidues()};const r=await fetch(API()+'/analyze_surface',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(payload),signal:AbortSignal.timeout(600000)});const d=await r.json();if(!r.ok)throw new Error(d.detail||'Backend error');data=d;document.getElementById('proteinInfo').innerHTML=`<strong>${d.settings.chain==='ALL'?'All chains':'Chain '+escapeHtml(d.settings.chain)}</strong> · ${d.stats.n_residues} residues · ${d.stats.n_surface_res} surface residues<br>InterfaceScout ${escapeHtml(d.version||'2.0.0-dev')} · Context: ${escapeHtml(d.settings.structure_context||'—')}`;document.getElementById('csvBtn').disabled=false;document.getElementById('pdbBtn').disabled=false;activeChem=d.chemistry_list[0];colorMode='propensity';renderAll();renderTargetProfile();renderApplicability();if(viewer){viewer.zoomTo(analysisSelector(),500);viewer.render();}}catch(e){showError(e.message)}finally{hideLoading()}}"
    html = html.replace(old_analyze, new_analyze, 1)

    html = html.replace('initViewer();checkBackend();loadMaterialProfiles();setInterval(checkBackend,10000);', 'initViewer();checkBackend();setInterval(checkBackend,10000);', 1)
    html = html.replace('initViewer();checkBackend();setInterval(checkBackend,10000);', 'initViewer();checkBackend();setInterval(checkBackend,10000);', 1)
    return html
