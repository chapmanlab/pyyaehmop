import sys
import pyyaehmop

def test_imported():
    assert 'pyyaehmop' in sys.modules

def test_benz(benzene_universe):
    assert len(benzene_universe.atoms) == 12


def test_run_bind(benzene_universe):
    ag = benzene_universe.atoms

    pyyaehmop.run_bind(ag.positions,
                       ag.types,
                       0.0)
