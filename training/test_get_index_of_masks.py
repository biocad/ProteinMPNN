from itertools import groupby
from pathlib import Path
import random
import string
from model_utils import get_mask_random,get_mask_cdrs
from proteinlib.structure.antibody_antigen_complex import AntibodyAntigenComplex, NumberingScheme

def test_mask_random():
    random.seed(42)
    len_seq=100
    max_parts=10
    max_length=5
    mask=get_mask_random(len_seq,max_parts,max_length)
    assert sum([k for k,_ in groupby(mask)])<=max_parts
    assert max(len(list(g)) for k, g in groupby(mask) if k==1) <=max_length
    assert mask== [0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    s=''.join((1+len_seq//len(string.ascii_uppercase))*string.ascii_uppercase)[:len_seq]
    s_masked=''.join(['_' if mask[i] else s[i] for i in range(len_seq)])
    assert s_masked=='ABC___GHIJKLMNOPQRSTUVWXYZABCDE__HIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUV'
  
def test_mask_cdrs():
    max_parts=10
    max_length=5

    pdb=Path('data/7pa6.pdb')
    heavy_chain_id="K"
    light_chain_id="k"
    antigen_chain_ids=["A"]
    numbering=NumberingScheme.CHOTHIA
    mask_vh,mask_vl=get_mask_cdrs(pdb=pdb,heavy_chain_id=heavy_chain_id,light_chain_id=light_chain_id,antigen_chain_ids=antigen_chain_ids,max_parts=max_parts,max_length=max_length,numbering=numbering)
    assert sum([k for k,_ in groupby(mask_vh)])<=max_parts
    assert sum([k for k,_ in groupby(mask_vl)])<=max_parts
    assert max(len(list(g)) for k, g in groupby(mask_vh) if k==1)<=max_length
    assert max(len(list(g)) for k, g in groupby(mask_vl) if k==1)<=max_length
    ab_complex = AntibodyAntigenComplex.from_pdb(pdb=pdb,heavy_chain_id=heavy_chain_id,light_chain_id=light_chain_id,antigen_chain_ids=antigen_chain_ids,numbering=numbering)
    r_h=ab_complex.antibody.heavy_chain.region_boundaries
    r_l=ab_complex.antibody.light_chain.region_boundaries
    assert len(r_h)==len(r_l)
    for i in range(0,len(r_h)-1,2):
        assert not any(mask_vh[r_h[i]:r_h[i+1]])
        assert not any(mask_vl[r_l[i]:r_l[i+1]])


test_mask_random()
test_mask_cdrs()
