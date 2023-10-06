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
 * lattice.h
 *
 * Author: Alexis Poncet
 * Email: alexis.poncet@ens-lyon.fr
 *
 * Header for lattice.cpp.
*/

#ifndef LATTICE_H
#define LATTICE_H

#include <vector>
#include <deque>
#include <array>
#include <cmath>
#include <mkl.h>
#include <mkl_vsl.h>

#define N_BITS_64 64
//#define USE_MAP

/*  This works only on Linux
#define DIV_ROUND_CLOSEST(x, divisor)(          \
{                           \
    typeof(x) __x = x;              \
    typeof(divisor) __d = divisor;          \
    (((typeof(x))-1) > 0 ||             \
     ((typeof(divisor))-1) > 0 || (__x) > 0) ?  \
        (((__x) + ((__d) / 2)) / (__d)) :   \
        (((__x) - ((__d) / 2)) / (__d));    \
}                           \
)
*/

#define DIV_ROUND_CLOSEST(n, d) ((((n) < 0) ^ ((d) < 0)) ? (((n) - (d)/2)/(d)) : (((n) + (d)/2)/(d)))

inline long modulo(long x, long N){
    return (x % N + N) % N;
}

inline long moduloSym(long x, long N){
    return x - N * DIV_ROUND_CLOSEST(x, N);
}

class Lattice {
	public:
		Lattice(long nx_, long ny_, bool periodic=true,
				bool correl=false, double dr_=1.);

		long getNVertices() const {
			return ntot;
		};
		long getNEdges() const {
			return n_edges;
		};
		const auto& getSitesOfEdge(long e) const {
			return sitesOfEdge[e];
		};
		double getMaxDist() const {
			return maxDist;
		};

		void print() const;
		void output(const std::string ofname) const;

		int checkIceRule(const std::vector<int> &spins) const;
		std::vector<std::vector<long> > 
			computeLoops(const std::vector<int> &spins) const;
		std::vector<std::vector<long> > 
			computeLoops(const std::vector<int> &spins,
	  				     std::vector<long> &loopOfEdge) const;
		long computeWindingX(const std::vector<long> &loop) const;
		long computeWindingY(const std::vector<long> &loop) const;
		double computeGirationRadius(const std::vector<long> &loop) const;
		long computeHeight(const std::vector<int>& spins) const;
		long makeWorm(const std::vector<int> &spins,
		       	      std::deque<long> &edgesVisited,
					  std::deque<int> &newSpins,
		              VSLStreamStatePtr stream) const;
		//long couplingEnergy(const std::vector<int> &spins) const;
		long couplingEnergyEdge(const long e, const std::vector<int> &spins)
		    	const;
		long diffCouplingEnergy(const std::vector<int> &spins1,
							    const std::vector<int> &spins2,
								const std::deque<long> &edgesVisited,
								const long start) const;
		long countZeroSites(const std::vector<int> &spins) const;
		double computeFracPara(const std::vector<int> &spins) const;
		void countAll(std::vector<long> &count) const;
		void countSameLoop(const std::vector<long> &loopOfEdge,
						   std::vector<long> &count) const;

	private:
		void makePeriodicLattice();
		void makeNonPeriodicLattice();

		const long nx, ny;
		const bool periodic; // Periodic boundary conditions or not
		const double sx, sy; // Scaling units in x and y directions
		const long Lx, Ly; // Lengths in scaled units
		const double dr; // Step for correlations
		long ntot, n_edges, n_faces; // Number of sites, number of edges

		std::vector<long> xOfSite;
		std::vector<long> yOfSite;
		std::vector<long> xOfEdge;
		std::vector<long> yOfEdge;
		std::vector<std::vector<long> > nbrsOfSite;
		std::vector<std::vector<long> > edgesOfSite;
		std::vector<std::array<long, 2> > sitesOfEdge;
		std::vector<std::vector<long> > nbrsOfEdge;
		std::vector<bool> edgeInBulk; // Only for non-periodic
		std::vector<std::vector<long> > edgesOfFace;
		std::vector<std::vector<long> > nbrsOfFace;
		std::vector<std::vector<int> > dirNbrOfFace;
		std::vector<long> edgesThroughFaces;
		std::vector<int> dirThroughFaces;

		double maxDist;
		std::vector<size_t> boxEdges;
};

int testLattice();
int testLatticeNotPeriodic();

#endif
