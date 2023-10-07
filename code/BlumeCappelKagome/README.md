## BlumeCappelKagome 

Monte-Carlo worm simulations of a custom Blume-Cappel-like model on the
honeycomb lattice. See [[Jorge et al. 2023](https://arxiv.org/abs/2305.06078)
for details about the model.

### Usage
#### Compilation
Run `make` to compile the program. You may need to modify the `Makefile`.
See [here](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-link-line-advisor.html)
for linking MKL properly.

#### Command-line options
Run `./BlumeCappelKagome --help` for information on the command-line options.
```
Usage: ./BlumeCappelKagome options
Options:
  --nx arg                   Number of sites in x direction
  --ny arg                   Number of sites in y direction
  -s [ --simuls ] arg (=1)   Number simulations
  --nbIters arg              Number of iterations
  --nbItersTh arg (=0)       Number of iterations of thermalization
  --skip arg (=-1)           Number of iterations before computation (default:
                             nx*ny)
  --K arg                    Coupling constant for flow
  --J arg (=0)               Coupling constant for orientation
  -o [ --output ] arg (=BCK) Output file
  -p [ --per ]               Use periodic boundary conditions
  --spins                    Output the spins
  --obs                      Output the observables
  -c [ --correl ]            Compute correlations in loop
  --dr arg (=0.25)           Step for correlations
  --test                     Test mode
  --lattice                  Output lattice
  -v [ --verbose ]           Verbose mode
  -h [ --help ]              Print help message and exit

```

#### Output files
Below, `output` is the string parameter given with `--output`.
- `output_nloops.dat` (if `--obs`): number of loops / number of winding loops
/ number of sites associated with three 0 / fraction of parallel contacts /
height of topographic map (for non-periodic simulations only)
(each line corresponds to one configuration)
- `output_len.dat` (if `--obs`): lengths of all the loops that were observed
- `output_wind.dat` (if `--obs`): windings of all the loops that were observed
(1 if the loop winds around the system, 0 otherwise)
- `output_Rg.dat` (if `--obs`): giration radii of all the loops that were
observed
- `output_correl.dat` (if `--obs` and `--correl`): loop correlations:
distance / number of edges in the same loop at distance / total number of
edges at distance (each line corresponds to one configuration=
- `output_spins.dat` (if `--spins`): spins at each iteration
(one configuration per line)

### Dependencies
- [Boost Program Options](https://www.boost.org/doc/libs/1_83_0/doc/html/program_options.html)
- [Intel Math Kernel Library (MKL)](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html)

### Licence
This code has been produced at [ENS de Lyon](https://www.ens-lyon.fr/).
It is released under the
[GNU General Public License](https://www.gnu.org/licenses/gpl-3.0.en.html).
