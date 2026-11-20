// Execute the real browser functions in a DOM-free VM; no browser dependencies.
const {test} = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const html = fs.readFileSync(path.join(__dirname, '../frontend/index.html'), 'utf8');
const script = html.match(/<script>([\s\S]*?)<\/script>/)[1];
const atom = 'ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C  ';
function runExport(raw, prepared, serial, rows) {
  const element = {classList: {toggle(){}}, style:{}};
  const ctx = vm.createContext({console, document:{getElementById:()=>element, addEventListener(){}},
    window:{addEventListener(){}}, fetch:async()=>({json:async()=>({engine:'test'})})});
  vm.runInContext(script, ctx);
  const response = {prepared_pdb:prepared, structure_preparation:{model_policy:'first_model_only',
    selected_model:1, selected_model_index:0, selected_model_serial:serial}};
  vm.runInContext(`pdbText=${JSON.stringify(raw)};data=${JSON.stringify(response)};`,ctx);
  return vm.runInContext(`rewriteBfactorPdb(${JSON.stringify(rows)},'test','canonical propensity')`,ctx);
}
const score = [{chain:'A',res_seq:1,icode:'',value:88.12345}];
test('B-factor quantization obeys the declared half-cent tolerance at score boundaries',()=>{
  const pdb=`${atom}\nEND\n`;
  for(const value of [0,1e-14,0.0049,0.0051,49.9949,49.9951,99.9949,99.9951,100]){
    const rows=[{chain:'A',res_seq:1,icode:'',value}];
    const line=runExport(pdb,pdb,null,rows).split('\n').find(l=>l.startsWith('ATOM'));
    assert.ok(Math.abs(Number(line.slice(60,66))-value)<=0.0050000001);
    assert.equal(line.slice(0,60),atom.slice(0,60));
    assert.equal(line.slice(66),atom.slice(66));
    assert.equal(rows[0].value,value,'export must not mutate the unrounded score');
  }
});
test('scores are exported only on the analyzed first model, with its original serial',()=>{
  const first = `MODEL        7\n${atom}\nENDMDL\nEND\n`;
  const raw = first.replace('END\n','') + `MODEL       42\n${atom}\nENDMDL\nEND\n`;
  const output = runExport(raw,first,7,score);
  assert.equal(output.split('\n').filter(l=>l.startsWith('ATOM')).length,1);
  assert.equal(output.split('\n').filter(l=>l.startsWith('MODEL')).length,1);
  assert.match(output,/MODEL\s+7/);
  assert.doesNotMatch(output,/MODEL\s+42/);
  assert.match(output,/FIRST MODEL ONLY/);
  assert.equal(output.split('\n').find(l=>l.startsWith('ATOM')).slice(60,66),' 88.12');
});
test('an implicit single model preserves coordinates and rounds only the B-factor field',()=>{
  const pdb = `${atom}\nEND\n`;
  const output = runExport(pdb,pdb,null,score);
  const rewritten=output.split('\n').find(l=>l.startsWith('ATOM'));
  assert.equal(rewritten.slice(0,60),atom.slice(0,60));
  assert.equal(rewritten.slice(66),atom.slice(66));
  assert.equal(rewritten.slice(60,66),' 88.12');
});
test('chain and insertion code remain part of residue identity',()=>{
  const inserted=atom.slice(0,26)+'A'+atom.slice(27);
  const otherChain=atom.slice(0,21)+'B'+atom.slice(22);
  const pdb=`${atom}\n${inserted}\n${otherChain}\nEND\n`;
  const rows=[...score,{chain:'A',res_seq:1,icode:'A',value:33},{chain:'B',res_seq:1,icode:'',value:22}];
  assert.deepEqual(runExport(pdb,pdb,null,rows).split('\n').filter(l=>l.startsWith('ATOM')).map(l=>+l.slice(60,66)),[88.12,33,22]);
});
test('missing prepared structure cannot silently fall back to the original ensemble',()=>{
  assert.throws(()=>runExport(`${atom}\nEND\n`,undefined,null,score),/prepared|analyzed/i);
});
test('mismatched model identity cannot receive scores',()=>{
  const pdb=`MODEL       42\n${atom}\nENDMDL\nEND\n`;
  assert.throws(()=>runExport(pdb,pdb,7,score),/model/i);
});
test('even repeated MODEL serials cannot authorize multiple models',()=>{
  const pdb=`MODEL        7\n${atom}\nENDMDL\nMODEL        7\n${atom}\nENDMDL\nEND\n`;
  assert.throws(()=>runExport(pdb,pdb,7,score),/model/i);
});
