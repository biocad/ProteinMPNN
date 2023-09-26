from itertools import groupby
import random
import string
from model_utils import get_index_of_masks
import numpy as np

def test_mask_indexing():
    random.seed(42)
    len_seq=100
    max_parts=10
    max_length=5
    mask=get_index_of_masks(len_seq,max_parts,max_length)
    assert sum([k for k,_ in groupby(mask)])<=max_parts
    assert max(len(list(g)) for k, g in groupby(mask) if k==1) <=max_length
    assert mask== [0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    s=''.join((1+len_seq//len(string.ascii_uppercase))*string.ascii_uppercase)[:len_seq]
    s_masked=''.join(['_' if mask[i] else s[i] for i in range(len_seq)])
    assert s_masked=='ABC___GHIJKLMNOPQRSTUVWXYZABCDE__HIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUV'

test_mask_indexing()
