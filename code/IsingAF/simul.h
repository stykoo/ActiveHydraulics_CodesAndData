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
 * simul.h
 *
 * Author: Alexis Poncet
 * Email: alexis.poncet@ens-lyon.fr
 *
 * Header for simul.cpp.
*/

#ifndef SIMUL_H
#define SIMUL_H

#include <vector>
#include <deque>
#include <mkl.h>
#include <mkl_vsl.h>
#include "parameters.h"
#include "lattice.h"

struct Observables {
	std::deque<long> numberOfLoops;
	std::deque<long> numberOfLoopsWind;

	std::deque<long> lenOfLoop;
	std::deque<bool> windOfLoop;
	std::deque<double> girationOfLoop;

	// For correlations
	std::vector<long> countAll;
	std::vector<long> countSameLoop;
};

int runSimulations(const Parameters &p);
void runSimulation(const Parameters &p, const Lattice &lattice,
		           Observables &obs, const unsigned int seed);
void outputObservables(const Parameters &p, const Observables &obs);

#endif
