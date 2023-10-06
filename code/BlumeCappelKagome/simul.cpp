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
 * simul.cpp
 *
 * Author: Alexis Poncet
 * Email: alexis.poncet@ens-lyon.fr
 *
 * Run the Monte-Carlo simulation of the model.
*/

#include <algorithm>
#include <fstream>
#include <sstream>
#include <iterator>
#include <thread>
#include <iomanip>
#include <random>
#include "simul.h"
#include "lattice.h"

// Run all the simulations and store the moments.
int runSimulations(const Parameters &p) {
	Observables obs;

	// Random numbers with MKL
	std::random_device rd;
	VSLStreamStatePtr stream;
	vslNewStream(&stream, VSL_BRNG_SFMT19937, rd());

	if (p.outputSpins) {
		std::ofstream file(p.output + ".dat");
		file << "# BlumeCappelKagome (" << __DATE__ <<  ", " << __TIME__
			<< ")\n";
		p.print(file);
		file.close();
	}

	Lattice lattice(p.nx, p.ny, p.periodic, p.computeCorrel, p.dr);

	// Initialize stuff for correlations
	if (p.computeCorrel) {
		size_t n_step_correl = (size_t) (lattice.getMaxDist() / p.dr + 1.);
		obs.countAll.assign(n_step_correl, 0);
		obs.countSameLoop.assign(n_step_correl, 0);
		lattice.countAll(obs.countAll);
	}

	// Launch simulations
	for (int i = 0 ; i < p.nbSimuls ; ++i) {
		if (p.verbose) {
			std::cout << "Simulation " << i << std::endl;
		}
		runSimulation(p, lattice, obs, stream); 
	}

	if (p.outputObs) {
		if (p.verbose)
			std::cout << "Outputing observables" << std::endl;
		outputObservables(p, obs);
	}

	return 0;
}

/*
void update(std::vector<int> &spins, const Lattice &lattice, const double K,
		    VSLStreamStatePtr stream) {
	std::deque<long> edgesVisited;
	std::deque<int> newSpins;

	// Make worm
	long n = lattice.makeWorm(spins, edgesVisited, newSpins, stream);

	// Explore worm starting from n
	double dE = 0.;
	for (size_t i = n ; i < edgesVisited.size() ; ++i) {
		double sOld = spins[edgesVisited[i]];
		double sNew = newSpins[i];
		dE += K * (sOld * sOld - sNew * sNew);
	}

	// Do worm with Monte-Carlo rule
	double u;
	vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, 1, &u, 0., 1.);
	if (dE < 0 || u < exp(-dE)) {
		for (size_t i = n ; i < edgesVisited.size() ; ++i)
			spins[edgesVisited[i]] = newSpins[i];
	}
}
*/

void assignNewSpins(
		std::vector<int> &spins, const std::deque<long> &edgesVisited,
		const std::deque<int> &newSpins, long n
	) {
	for (size_t i = n ; i < edgesVisited.size() ; ++i) {
		spins[edgesVisited[i]] = newSpins[i];
	}
}

// From spins2 to spins1
void copySpins(
		std::vector<int> &spins1, const std::vector<int> &spins2,
		const std::deque<long> &edgesVisited, long n
	) {
	for (size_t i = n ; i < edgesVisited.size() ; ++i) {
		spins1[edgesVisited[i]] = spins2[edgesVisited[i]];
	}
}

void update(std::vector<int> &spins, std::vector<int> &spins_cpy,
		    const Lattice &lattice, const Parameters &p,
			VSLStreamStatePtr stream) {
	// At the beginning of the function spins and spins_cpy should be identical
	
	std::deque<long> edgesVisited;
	std::deque<int> newSpins;

	// Make worm
	long n = lattice.makeWorm(spins, edgesVisited, newSpins, stream);
	double dE = 0.;

	// Flow energy
	for (size_t i = n ; i < edgesVisited.size() ; ++i) {
		double sOld = spins[edgesVisited[i]];
		double sNew = newSpins[i];
		dE += p.K * (sOld * sOld - sNew * sNew);
	}

	// Coupling energy (only if J != 0)
	if (p.withCouplings) {
		assignNewSpins(spins_cpy, edgesVisited, newSpins, n);
		dE += p.J * lattice.diffCouplingEnergy(spins, spins_cpy,
				edgesVisited, n);
	}

	// Do worm with Monte-Carlo rule
	double u;
	vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, 1, &u, 0., 1.);
	if (dE < 0 || u < exp(-dE)) {
		// If move is accepted, assign the new spins
		assignNewSpins(spins, edgesVisited, newSpins, n);
	} else if (p.withCouplings) {
		// If move is rejected, enforce that spins and spins_cpy remain identical
		copySpins(spins_cpy, spins, edgesVisited, n);
	}
}

void runSimulation(const Parameters &p, const Lattice &lattice,
		           Observables &obs, VSLStreamStatePtr stream) {
	std::vector<int> spins(lattice.getNEdges(), 0); // Start with 0 everywhere
	std::vector<int> spins_cpy(spins);

	// Thermalization
	for (long k = 0 ; k < p.nbItersTh ; k++) {
		if (p.verbose)
			std::cout << "Therm: " << k << " / " << p.nbItersTh
				<< "\r" << std::flush;

		for (long l = 0 ; l < p.skip ; l++) {
			update(spins, spins_cpy, lattice, p, stream);
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

	// Main loop
	for (long k = 0 ; k < p.nbIters ; k++) {
		if (p.verbose)
			std::cout << "Main: " << k << " / " << p.nbIters
				<< "\r" << std::flush;

		for (long l = 0 ; l < p.skip ; l++) {
			update(spins, spins_cpy, lattice, p, stream);
		}

		if (p.outputSpins) {
			for (long l = 0 ; l < lattice.getNEdges() ; l++) {
				file << spins[l] << " ";
			}
			file << "\n";
		}

		// Loop stuff
		if (p.outputObs) {
			auto edgesOfLoop = p.computeCorrel ?
				lattice.computeLoops(spins, loopOfEdge) :
				lattice.computeLoops(spins);
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
			obs.numberOfZeroSites.push_back(lattice.countZeroSites(spins));
			obs.fracPara.push_back(lattice.computeFracPara(spins));
			if (!p.periodic) {
				obs.height.push_back(lattice.computeHeight(spins));
			}

			if (p.computeCorrel) {
				lattice.countSameLoop(loopOfEdge, obs.countSameLoop);
			}
		}
	}

	if (p.verbose)
		std::cout << "\n";

	if (p.outputSpins)
		file.close();
}

void outputObservables(const Parameters &p, const Observables &obs) {
	std::ofstream file;
	std::string header = "# BlumeCappelKagome (";
	header += __DATE__;
	header +=  ", ";
	header += __TIME__;
	header += ")\n";

	// Number of loops
	file.open(p.output + "_nloops.dat");
	file << header;
	p.print(file);
	file << "# n_loops n_loops_wind n_zero_sites frac_para (height)\n";
	for (size_t i = 0 ; i < obs.numberOfLoops.size() ; ++i) {
		file << obs.numberOfLoops[i] << " " << obs.numberOfLoopsWind[i] << " "
			<< obs.numberOfZeroSites[i] << " " << obs.fracPara[i];
		if (!p.periodic) {
			file << " " << obs.height[i];
		}
		file << "\n";
	}
	file.close();

	// Length of loops
	file.open(p.output + "_len.dat");
	file << header;
	p.print(file);
	for (auto l : obs.lenOfLoop) {
		file << l << "\n";
	}
	file.close();

	// Winding of loops
	if (p.periodic) {
		file.open(p.output + "_wind.dat");
		file << header;
		p.print(file);
		for (auto l : obs.windOfLoop) {
			file << l << "\n";
		}
		file.close();
	}

	// Giration radii of loops
	file.open(p.output + "_Rg.dat");
	file << header;
	p.print(file);
	for (auto l : obs.girationOfLoop) {
		file << l << "\n";
	}
	file.close();

	// Loop correlations
	if (p.computeCorrel) {
		file.open(p.output + "_correl.dat");
		file << header;
		p.print(file);
		for (size_t i = 0 ; i < obs.countAll.size() ; ++i) {
			file << i * p.dr << " " << obs.countSameLoop[i]
				<< " " << obs.countAll[i] << "\n";
		}
		file.close();
	}
}
