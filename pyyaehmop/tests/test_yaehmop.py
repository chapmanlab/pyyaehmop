import sys
import pyyaehmop

def test_imported():
    assert 'pyyaehmop' in sys.modules

def test_benz(benzene_universe):
    assert len(benzene_universe.atoms) == 12
