/* Optional real high-resolution viewer capture for publication graphics.
   Separate user action: no change to all-11-map research ZIP semantics. */
let publicationImageUrl=null;
let publicationBusy=false;
async function exportPublicationPng(){
  const button=document.getElementById("publicationImageBtn");
  const link=document.getElementById("publicationImageReady");
  const status=document.getElementById("publicationImageStatus");
  if(publicationBusy)return;
  if(!data||!viewer||activeType!=="chemistry"||!activeMap){
    status.textContent="Analyze a protein and open one of its 11 individual maps first.";
    return;
  }
  publicationBusy=true;
  button.disabled=true;
  link.style.display="none";
  status.textContent="Rendering actual 3D structure at publication-image resolution…";
  const node=document.getElementById("viewer");
  const styles={width:node.style.width,height:node.style.height};
  const oldView=typeof viewer.getView==="function"?viewer.getView():null;
  try{
    node.style.width="1920px";
    node.style.height="1080px";
    viewer.resize();
    if(oldView&&typeof viewer.setView==="function")viewer.setView(oldView);
    const uri=await ensembleCaptureScreenshot();
    const picture=new Image();
    picture.src=uri;
    await picture.decode();
    if(picture.naturalWidth<1600||picture.naturalHeight<900){
      throw Error("Browser rendered only "+picture.naturalWidth+" × "+picture.naturalHeight+
        " px. Export not labeled high-resolution; increase WebGL capability or use full-size PNG in the complete ZIP.");
    }
    const blob=await (await fetch(uri)).blob();
    if(publicationImageUrl)URL.revokeObjectURL(publicationImageUrl);
    publicationImageUrl=URL.createObjectURL(blob);
    link.href=publicationImageUrl;
    link.download="InterfaceScout_"+safe(structureName)+"_"+safe(activeMap)+
      "_Publication_"+picture.naturalWidth+"x"+picture.naturalHeight+".png";
    link.textContent="Save publication PNG · "+picture.naturalWidth+" × "+picture.naturalHeight;
    link.style.display="block";
    status.textContent="Real 3Dmol render, no interpolation. Use the separate figure caption to identify map, protein and scoring convention.";
  }catch(error){
    status.textContent=error.message;
  }finally{
    node.style.width=styles.width;
    node.style.height=styles.height;
    viewer.resize();
    if(oldView&&typeof viewer.setView==="function")viewer.setView(oldView);
    viewer.render();
    publicationBusy=false;
    button.disabled=false;
  }
}
