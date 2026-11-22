"""Analytic invariance and nonredundancy checks; no experimental adsorption claims."""
from __future__ import annotations
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"backend"))
from regional_chemistry import annotate_regions

def template(a,b,coords):
 keys=list(coords)
 maps={}
 for name,scores in [("hydrophobic",a),("hbond_acceptor",b)]:
  maps[name]={"residues":[{"key":k,"local_score":v} for k,v in scores.items() if v>0]}
 return {
 "maps":maps,"diagnostics":{"surface_ca_xyz_A":coords}
 }
def rows(*members):
 return [{"member_keys":list(group),"source_center_key":group[0]} for group in members]

class InfoTests(unittest.TestCase):
 def test_support_strength_is_not_coverage(self):
  # Two patches with same number of positive supports, but different exposures.
  coords={f"A:{i}:":[float(i),0,0] for i in range(1,5)}
  m=template({"A:1:":.1,"A:3:":.9},{"A:2:":.1,"A:4:":.9},coords)
  x=rows(["A:1:","A:2:"],["A:3:","A:4:"])
  annotate_regions(m,["hydrophobic","hbond_acceptor"],x)
  pa=x[0]["regional_profile"];pb=x[1]["regional_profile"]
  self.assertEqual(pa["supporting_residue_fraction_by_map"],pb["supporting_residue_fraction_by_map"])
  self.assertNotEqual(pa["mean_local_support_by_map"],pb["mean_local_support_by_map"])
 def test_proximity_is_not_region_membership_or_chemistry(self):
  # Identical residue sets, scores, count, coverage, mean; only CA geometry differs.
  keys=["A:1:","A:2:","A:3:"]
  scores1={"A:1:":1.0};scores2={"A:2:":1.0}
  a=template(scores1,scores2,{keys[0]:[0,0,0],keys[1]:[2,0,0],keys[2]:[4,0,0]})
  b=template(scores1,scores2,{keys[0]:[0,0,0],keys[1]:[7,0,0],keys[2]:[4,0,0]})
  p=rows(keys);q=rows(keys)
  annotate_regions(a,["hydrophobic","hbond_acceptor"],p)
  annotate_regions(b,["hydrophobic","hbond_acceptor"],q)
  u=p[0]["regional_profile"];v=q[0]["regional_profile"]
  self.assertEqual(u["mean_local_support_by_map"],v["mean_local_support_by_map"])
  self.assertEqual(u["supporting_residue_fraction_by_map"],v["supporting_residue_fraction_by_map"])
  self.assertEqual(u["min_distinct_intermap_CA_distance_A"],2)
  self.assertEqual(v["min_distinct_intermap_CA_distance_A"],7)
 def test_contrast_is_affine_copy_not_independent(self):
  coords={f"A:{i}:":[i,0,0] for i in range(1,5)}
  a=template({"A:1:":.2,"A:3:":.8},{"A:2:":.3,"A:4:":.7},coords)
  p=rows(["A:1:","A:2:"],["A:3:","A:4:"])
  annotate_regions(a,["hydrophobic","hbond_acceptor"],p)
  for c in ["hydrophobic","hbond_acceptor"]:
   delta_m=p[1]["regional_profile"]["mean_local_support_by_map"][c]-p[0]["regional_profile"]["mean_local_support_by_map"][c]
   delta_c=p[1]["regional_profile"]["relative_to_candidate_median_by_map"][c]-p[0]["regional_profile"]["relative_to_candidate_median_by_map"][c]
   self.assertAlmostEqual(delta_m,delta_c)
 def test_label_independence(self):
  coords={"A:1:":[0,0,0],"A:2:":[1,0,0]}
  m=template({"A:1:":.7},{"A:2:":.4},coords)
  p=rows(["A:1:","A:2:"])
  x=annotate_regions(m,["hydrophobic","hbond_acceptor"],p)
  m["experimental_contacts"]=["A:1:"]
  q=rows(["A:1:","A:2:"])
  y=annotate_regions(m,["hydrophobic","hbond_acceptor"],q)
  self.assertEqual(x,y);self.assertEqual(p,q)

if __name__=="__main__":unittest.main()
