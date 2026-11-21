/* Progressive disclosure only. Canonical computations and exports are unchanged. */
(function(){
  const side=document.querySelector("aside");
  if(!side)return;
  const sections=Array.from(side.querySelectorAll(".section"));
  const get=(name)=>sections.find(section=>section.querySelector("h3")?.textContent.toLowerCase().includes(name));
  const panels={
    protein:get("protein"),solution:get("solution"),research:get("research question"),
    compare:get("compare conformations"),multi:get("multi-pdb"),legacy:get("export"),
    complete:get("complete export")
  };
  if(!panels.protein||!panels.solution||!panels.multi||!panels.complete)return;
  const head=document.createElement("div");
  head.className="workspace-top";
  head.innerHTML='<div class="workspace-eyebrow">INTERFACESCOUT · WORKSPACE</div>'+
    '<h2 class="workspace-title">Explore what matters.</h2>'+
    '<p class="workspace-subtitle">Choose a workflow. Your 11 chemistry maps are always calculated together.</p>'+
    '<nav class="workspace-modes" aria-label="Analysis workflow">'+
    '<button class="workspace-mode" type="button" data-workspace-target="single" aria-pressed="true">Single protein</button>'+
    '<button class="workspace-mode" type="button" data-workspace-target="compare" aria-pressed="false">Compare two</button>'+
    '<button class="workspace-mode" type="button" data-workspace-target="multi" aria-pressed="false">Multi-PDB</button></nav>';
  side.insertBefore(head,side.firstChild);
  const exportPanel=document.createElement("section");
  exportPanel.className="workspace-export";
  exportPanel.setAttribute("aria-label","Complete results download");
  exportPanel.innerHTML='<h3>Your complete research package</h3>'+
    '<p>All 11 maps, always. PDFs with actual 3D screenshots, PNG images, B-factor PDBs, Excel, CSV, and original inputs.</p>';
  const full=$("isFullExport");
  full.textContent="Download complete results · ZIP";
  exportPanel.appendChild(full);
  const status=$("isStatus"), error=$("isError");
  if(status)exportPanel.appendChild(status);
  if(error)exportPanel.appendChild(error);
  const more=document.createElement("details");
  more.className="workspace-more";
  more.innerHTML="<summary>Individual files and advanced downloads</summary>";
  if(panels.legacy)more.appendChild(panels.legacy);
  const researchReport=$("decisionExport");
  if(researchReport)more.appendChild(researchReport);
  exportPanel.appendChild(more);
  if(panels.complete)panels.complete.remove();
  side.insertBefore(exportPanel,side.querySelector("#error")||null);
  const condition=panels.solution;
  const row=condition.querySelector(".row");
  if(row){
    const advanced=document.createElement("details");
    advanced.className="workspace-settings";
    advanced.innerHTML="<summary>Advanced conditions · recorded as metadata</summary>";
    const ionic=row.querySelectorAll(".field")[1];
    const temp=$("temp")?.closest(".field");
    if(ionic)advanced.appendChild(ionic);
    if(temp)advanced.appendChild(temp);
    const note=condition.querySelector(".note");
    if(note)advanced.appendChild(note);
    row.insertAdjacentElement("afterend",advanced);
    row.style.gridTemplateColumns="1fr";
  }
  if(panels.protein){
    const h=panels.protein.querySelector("h3");if(h)h.textContent="Protein structure";
  }
  if(panels.solution){
    const h=panels.solution.querySelector("h3");if(h)h.textContent="Solution pH";
  }
  if(panels.compare){
    const h=panels.compare.querySelector("h3");if(h)h.textContent="Second structure";
    panels.compare.insertAdjacentHTML("afterbegin","<p class='workspace-quiet'>Load and analyze the first structure above, then choose the second structure. The complete export includes both.</p>");
  }
  if(panels.multi){
    const h=panels.multi.querySelector("h3");if(h)h.textContent="Structure collection";
  }
  const titles={
    single:["Explore a protein.","Begin with a PDB ID or a local structure. Every chemistry map is calculated in one analysis."],
    compare:["Compare two structures.","Inspect matched residues, surface changes and all 11 chemistry maps for both PDBs."],
    multi:["Explore a collection.","Add 2–20 structures, inspect each independently and compare matched positions where appropriate."]
  };
  let current="single";
  function mode(next){
    if(!titles[next])return;
    const prev=current;current=next;
    document.body.dataset.workspaceMode=next;
    head.querySelector(".workspace-title").textContent=titles[next][0];
    head.querySelector(".workspace-subtitle").textContent=titles[next][1];
    for(const button of head.querySelectorAll("[data-workspace-target]")){
      button.setAttribute("aria-pressed",String(button.dataset.workspaceTarget===next));
    }
    const visible={
      single:new Set(["protein","solution","research"]),
      compare:new Set(["protein","solution","compare"]),
      multi:new Set(["solution","multi"])
    }[next];
    for(const [key,panel] of Object.entries(panels)){
      if(!panel||["legacy","complete"].includes(key))continue;
      panel.classList.toggle("workspace-hidden",!visible.has(key));
    }
    $("run")?.classList.toggle("workspace-hidden",next==="multi");
    const results=$("isSetResults");
    if(results)results.classList.toggle("workspace-hidden",next==="single");
    if(prev==="multi"&&next!=="multi"&&typeof ensembleResult!=="undefined"){
      ensembleResult=null;ensembleAnalyzedInputs=[];ensembleConsensusKey=null;
      if(results)results.innerHTML="";
      if(typeof ensembleRefresh==="function")ensembleRefresh();
    }
    if(typeof viewer!=="undefined"&&viewer&&typeof viewer.resize==="function"){
      requestAnimationFrame(()=>viewer.resize());
    }
  }
  head.querySelectorAll("[data-workspace-target]").forEach(button=>
    button.addEventListener("click",()=>mode(button.dataset.workspaceTarget)));
  mode("single");
  window.workspaceSetMode=mode;
})();
