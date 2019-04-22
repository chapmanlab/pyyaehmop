import pytest

import MDAnalysis as mda
from MDAnalysis.lib.util import NamedStream
import io

benz = """\
CRYST1    0.000    0.000    0.000   0.00   0.00   0.00 P 1           1
ATOM      1  C1  BEN     1      -0.704  -1.218   0.000  1.00  0.00           C
ATOM      2  H1  BEN     1      -1.244  -2.155   0.000  1.00  0.00           H
ATOM      3  C2  BEN     1       0.704  -1.218   0.000  1.00  0.00           C
ATOM      4  H2  BEN     1       1.244  -2.155   0.000  1.00  0.00           H
ATOM      5  C3  BEN     1       1.406   0.000   0.000  1.00  0.00           C
ATOM      6  H3  BEN     1       2.488   0.000   0.000  1.00  0.00           H
ATOM      7  C4  BEN     1       0.704   1.218   0.000  1.00  0.00           C
ATOM      8  H4  BEN     1       1.244   2.155   0.000  1.00  0.00           H
ATOM      9  C5  BEN     1      -0.704   1.218   0.000  1.00  0.00           C
ATOM     10  H5  BEN     1      -1.244   2.155   0.000  1.00  0.00           H
ATOM     11  C6  BEN     1      -1.406   0.000   0.000  1.00  0.00           C
ATOM     12  H6  BEN     1      -2.488   0.000   0.000  1.00  0.00           H
CONECT    1    2    3   11
CONECT    2    1
CONECT    3    1    4    5
CONECT    4    3
CONECT    5    3    6    7
CONECT    6    5
CONECT    7    5    8    9
CONECT    8    7
CONECT    9    7   10   11
CONECT   10    9
CONECT   11    1    9   12
CONECT   12   11
END"""

@pytest.fixture
def benzene_universe():
    return mda.Universe(NamedStream(io.StringIO(benz), 'benzene.pdb'))
