// End-to-end exercise of the shipped frontend, checked by a separate PDB parser.
const fs=require('node:fs');
const vm=require('node:vm');
const path=require('node:path');
const crypto=require('node:crypto');
const assert=require('node:assert/strict');
const result=JSON.parse(fs.readFileSync(process.argv[2],'utf8'));
const html=fs.readFileSync(path.join(__dirname,'../../frontend/index.html'),'utf8');
const source=html.match(/<script>([\s\S]*?)<\/script>/)[1];
const element={classList:{toggle(){}},style:{}};
const context=vm.createContext({console,document:{getElementById:()=>element,addEventListener(){}},
  window:{addEventListener(){}},fetch:async()=>({json:async()=>({engine:'verification'})})});
vm.runInContext(source,context);
context.input=result;
vm.runInContext('data=input;pdbText=input.prepared_pdb;',context);
const isAtom=line=>line.startsWith('ATOM  ')||line.startsWith('HETATM');
const prepared=result.prepared_pdb.replace(/\r/g,'').split('\n');
const atoms=prepared.filter(isAtom);
const report={maps:[],checks:0};
for(const key of result.map_order){
  context.key=key;
  const text=vm.runInContext("activeType='chemistry';activeMap=key;(()=>{const r=bfactorRows();return rewriteBfactorPdb(r.rows,r.label,r.kind)})()",context);
  const lines=text.split('\n');
  const exported=lines.filter(isAtom);
  assert.equal(exported.length,atoms.length,'atom count');
  assert.deepEqual(lines.filter(l=>l.startsWith('MODEL ')),prepared.filter(l=>l.startsWith('MODEL ')),'model identity');
  const expected=new Map(result.maps[key].residues.map(r=>[`${r.chain.trim()}:${r.res_seq}:${r.icode.trim()}`,r.propensity]));
  let maxError=0;
  for(let i=0;i<atoms.length;i++){
    const a=atoms[i],b=exported[i];
    assert.equal(a.slice(0,60),b.slice(0,60),'identity/coordinates/occupancy');
    assert.equal(a.slice(66),b.slice(66),'non-B-factor suffix');
    const id=`${a.slice(21,22).trim()}:${Number(a.slice(22,26))}:${a.slice(26,27).trim()}`;
    const value=expected.get(id)??0;
    const field=Number(b.slice(60,66));
    const error=Math.abs(field-value);
    assert.ok(Number.isFinite(field)&&field>=0&&field<=100&&error<=0.0050000001,`${key}/${id}: B-factor precision`);
    maxError=Math.max(maxError,error);report.checks+=3;
  }
  report.maps.push({key,atoms:atoms.length,max_abs_bfactor_error:maxError,
    sha256:crypto.createHash('sha256').update(text).digest('hex')});
}
assert.equal(report.maps.length,11);
process.stdout.write(JSON.stringify(report));
