/*
BlumeCappelKagome, custom Blume-Cappel-like model on the honeycomb lattice.
Copyright (C) 2023 ENS de Lyon
Contributor: Alexis Poncet <alexis.poncet@ens-lyon.fr>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/
/*
 * BlumeCappelKagome
 * main.cpp
 *
 * Author: Alexis Poncet
 * Email: alexis.poncet@ens-lyon.fr
 *
 * Simulate a Blume-Cappel-like model on the edges of the honeycomb lattice.
*/

#include <iostream>
#include "parameters.h"
#include "lattice.h"
#include "simul.h"

int main(int argc, char **argv) {
	// Load parameters
    Parameters p;
    int status = p.fromCommandLine(argc, argv);
	if (status) {
		return status - 1;
	}
	status = p.check();
	if (status) {
		return status;
	}

	if (p.test) {
		status = testLattice();
		testLatticeNotPeriodic();
		return status;
	}
	if (p.outputLattice) {
		Lattice lattice(p.nx, p.ny, p.periodic, p.computeCorrel, p.dr);
		lattice.output(p.output);
		return 0;
	}

    if (p.verbose) {
		std::cout << "--- This is " << argv[0] << ", compiled on " << __DATE__
			<< " at " << __TIME__ << " ---" << std::endl;
        p.print();
		std::cout << std::endl;
    }

	status = runSimulations(p);

	return status;
}
