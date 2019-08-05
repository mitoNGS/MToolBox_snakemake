from modules.general import nuc

def test_nuc():
    assert nuc("AACCCGT") == {'A': 2, 'C': 3, 'G': 1, 'T': 1, 'N':0}

