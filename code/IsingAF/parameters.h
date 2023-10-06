/*
IsingAF, Anti-ferromagnetic Ising model on the triangular lattice.
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
 * IsingAF
 * parameters.h
 *
 * Author: Alexis Poncet
 * Email: alexis.poncet@ens-lyon.fr
 *
 * Header for parameters.cpp.
*/

#ifndef PARAMETERS_H
#define PARAMETERS_H

#define DEFAULT_OUTPUT_FILE "ising_AF"

#include <iostream>
#include <string>
#include <vector>

struct Parameters {
	// Methods
	int check() const;
	void print(std::ostream &stream = std::cout) const;
	int fromCommandLine(int argc, char **argv);

	long nx;  // Number of hexagons in x direction
	long ny;  // Number of hexagons in x direction
	long nbSimuls;  // Number of simulations
	long nbIters; // Number of iterations
	long nbItersTh; // Number of iterations of thermalization
	double K; // Coupling constant

	std::string output;  // Name of the output file

	long ntot; // Total number of sites

	bool outputSpins; // Output the raw spins
	bool outputObs; // Output the observables
	bool computeCorrel;  // Compute loop correlations
	double dr; // Step for correlations
	bool test;  // Test mode
	bool verbose;  // Verbose mode
};

#endif
