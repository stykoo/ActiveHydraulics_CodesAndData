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
 * simul.cpp
 *
 * Author: Alexis Poncet
 * Email: alexis.poncet@ens-lyon.fr
 *
 * Run the simulation of the Ising model.
*/

#include <algorithm>
#include <fstream>
#include <sstream>
#include <iterator>
#include <thread>
#include <iomanip>
#include <random>
#include "simul.h"

// Run all the simulations and store the moments.
int runSimulations(const Parameters &p) {
	std::random_device rd;
	Observables obs;

	if (p.outputSpins) {
		std::ofstream file(p.output + ".dat");
		file << "# IsingAF (" << __DATE__ <<  ", " << __TIME__ << ")\n";
		p.print(file);
		file.close();
	}

	Lattice lattice(p.nx, p.ny, p.computeCorrel, p.dr);

	// Initialize stuff for correlations
	if (p.computeCorrel) {
		size_t n_step_correl = (size_t) (lattice.getMaxDist() / p.dr + 1.);
		obs.countAll.assign(n_step_correl, 0);
		obs.countSameLoop.assign(n_step_correl, 0);
		lattice.countAll(obs.countAll);
	}

	// Launch simulations
	for (int i = 0 ; i < p.nbSimuls ; ++i) {
		runSimulation(p, lattice, obs, rd()); 
	}

	if (p.outputObs) {
		if (p.verbose)
			std::cout << "Outputing observables" << std::endl;
		outputObservables(p, obs);
	}

	return 0;
}

void update(std::vector<int> &spins, const Lattice &lattice, const double K,
			int i, double u) {
	double dE = 0.;
	for (long j : lattice.getNbrsOfSite(i)) {
		dE += spins[j];
	}
	dE *= -2 * K * spins[i];
	if (dE < 0 || u < -dE) {
		spins[i] *= -1;
	}
}

void runSimulation(const Parameters &p, const Lattice &lattice,
		           Observables &obs, const unsigned int seed) {
	// Seed useful if simulations are ran on different threads
	// (useless in the present case)

	// Random numbers with MKL
	VSLStreamStatePtr stream;
	vslNewStream(&stream, VSL_BRNG_SFMT19937, seed);
	std::vector<int> indices(p.ntot);
	std::vector<double> us(p.ntot);

	std::vector<bool> edgeOnWall;
	std::vector<int> spins(p.ntot);
	viRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, p.ntot, spins.data(),
				 0, 2);
	for (int i = 0 ; i < p.ntot ; ++i) {
		spins[i] = 2 * spins[i] - 1;
	}

	// Thermalization
	for (long k = 0 ; k < p.nbItersTh ; k++) {
		if (p.verbose)
			std::cout << k << " / " << p.nbItersTh << "\r" << std::flush;

		viRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, p.ntot, indices.data(),
					 0, p.ntot);
		vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, p.ntot, us.data(),
				     0.0, 1.0);
		vdLn(p.ntot, us.data(), us.data());

		for (long l = 0 ; l < p.ntot ; l++) {
			update(spins, lattice, p.K, indices[l], us[l]);
		}
	}

	if (p.verbose)
		std::cout << "\n";

	std::ofstream file;
	if (p.outputSpins)
		file.open(p.output + ".dat", std::ios_base::app); // Append

	// For the computation of correlations
	std::vector<long> loopOfEdge;
	if (p.computeCorrel) {
		loopOfEdge.resize(lattice.getNEdges());
	}

	for (long k = 0 ; k < p.nbIters ; k++) {
		if (p.verbose)
			std::cout << k << " / " << p.nbIters << "\r" << std::flush;

		viRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, p.ntot, indices.data(),
					 0, p.ntot);
		vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, p.ntot, us.data(),
				     0.0, 1.0);
		vdLn(p.ntot, us.data(), us.data());

		for (long l = 0 ; l < p.ntot ; l++) {
			update(spins, lattice, p.K, indices[l], us[l]);
		}

		if (p.outputSpins) {
			for (long l = 0 ; l < p.ntot ; l++) {
				file << spins[l] << " ";
			}
			file << "\n";
		}

		// Loop stuff
		if (p.outputObs) {
			lattice.computeWalls(edgeOnWall, spins);
			auto edgesOfLoop = p.computeCorrel ?
				lattice.computeLoops(edgeOnWall, loopOfEdge) :
				lattice.computeLoops(edgeOnWall);

			long n_loops = edgesOfLoop.size();
			long n_loops_wind = 0;
			for (auto loop : edgesOfLoop) {
				long wx = lattice.computeWindingX(loop);
				long wy = lattice.computeWindingY(loop);
				bool wind = (wx != 0) || (wy != 0);
				double Rg = lattice.computeGirationRadius(loop);
				obs.lenOfLoop.push_back(loop.size());
				obs.windOfLoop.push_back(wind);
				obs.girationOfLoop.push_back(Rg);
				if (wind)
					n_loops_wind++;
			}
			obs.numberOfLoops.push_back(n_loops);
			obs.numberOfLoopsWind.push_back(n_loops_wind);

			if (p.computeCorrel) {
				lattice.countSameLoop(loopOfEdge, obs.countSameLoop);
			}
		}
	}

	if (p.outputSpins)
		file.close();
}

void outputObservables(const Parameters &p, const Observables &obs) {
	std::ofstream file;

	// Number of loops
	file.open(p.output + "_nloops.dat");
	file << "# IsingAF (" << __DATE__ <<  ", " << __TIME__ << ")\n";
	p.print(file);
	file << "# n_loops n_loops_wind\n";
	for (size_t i = 0 ; i < obs.numberOfLoops.size() ; ++i) {
		file << obs.numberOfLoops[i] << " " << obs.numberOfLoopsWind[i] << "\n";
	}
	file.close();

	// Length of loops
	file.open(p.output + "_len.dat");
	file << "# IsingAF (" << __DATE__ <<  ", " << __TIME__ << ")\n";
	p.print(file);
	for (auto l : obs.lenOfLoop) {
		file << l << "\n";
	}
	file.close();

	// Winding of loops
	file.open(p.output + "_wind.dat");
	file << "# IsingAF (" << __DATE__ <<  ", " << __TIME__ << ")\n";
	p.print(file);
	for (auto l : obs.windOfLoop) {
		file << l << "\n";
	}
	file.close();

	// Giration radii of loops
	file.open(p.output + "_Rg.dat");
	file << "# IsingAF (" << __DATE__ <<  ", " << __TIME__ << ")\n";
	p.print(file);
	for (auto l : obs.girationOfLoop) {
		file << l << "\n";
	}
	file.close();

	// Loop correlations
	if (p.computeCorrel) {
		file.open(p.output + "_correl.dat");
		file << "# IsingAF (" << __DATE__ <<  ", " << __TIME__ << ")\n";
		p.print(file);
		for (size_t i = 0 ; i < obs.countAll.size() ; ++i) {
			file << i * p.dr << " " << obs.countSameLoop[i]
				<< " " << obs.countAll[i] << "\n";
		}
		file.close();
	}
}
