## FiniteElement 

Codes for finite element simulations of the Toner-Tu equations in oblong 
geometry. The code makes use of the Python interface to
[FEniCS](https://fenicsproject.org/) for the finite element method.

The authors of the article acknowledge the earlier contribution of
Yoann Poupart to the making of this code.

### Usage
Modify the parameters in `main_simu_oblong_parallel_v2.py` and
run `python3 main_simu_oblong_parallel_v2.py`.

### Files
- `simu_oblong_parallel_v2.py`:
Main file to launch finite element simulations of the Toner-Tu equations
in oblong geometry.
the honeycomb lattice.
- `sub_simu_oblong_parallel_v2.py`:
Module for finite element simulations of the Toner-Tu equations
in oblong geometry.
- `sub_sub_simu_oblong_parallel_v2.py`:
Module using the python interface to FEniCS for finite element simulations of
the Toner-Tu equations in oblong geometry.

### Licence
This code has been produced at [ENS de Lyon](https://www.ens-lyon.fr/).
It is released under the
[GNU General Public License](https://www.gnu.org/licenses/gpl-3.0.en.html).
