# pyyaehmop

Python bindings to the [Yaehmop](https://github.com/greglandrum/yaehmop) package.

Requires:
 - Python 3
 - Cython
 - Yaehmop installed (including latest lib & headers)
 
With these in place, installation should just be `python setup.py install`

## Usage

Usage looks like:

```python
import yaehmop

# Can get atomic data from anywhere, such as MDAnalysis
u = mda.Universe(..., ...)
atomic_positions = u.atoms.positions
atomic_names = u.atoms.names

# Currently only return Hamiltonian and overlap matrices
# Hacking required for other data, but not impossible...
H, S = yaehmop.run_bind(atomic_positions, atomic_names, charge=0.0)

```
