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
#include <array>
#include <cmath>

// Maximum size for which we store the data for the correlations (boxEdges)
// The size of box edges is (9/2) N^4
// For N=100 -> 5e8 
#define NMAX_CORREL 100

/*
 * Divide positive or negative dividend by positive divisor and round
 * to closest integer. Result is undefined for negative divisors and
 * for negative dividends if the divisor variable type is unsigned.
 */
/*
This works only on Linux
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

/*template<typename T>
inline T pbc(const T x, const T L){
	return x - L * std::floor(x / L);
}

template<typename T>
inline T pbcSym(const T x, const T L) {
	return x - L * std::round(x / L);
}*/

class Lattice {
	public:
		Lattice(long nx_, long ny_, bool correl=false, double dr_=1.);

		long getNEdges() const {
			return n_edges;
		};
		double getMaxDist() const {
			return maxDist;
		};
		const auto& getNbrsOfSite(long i) const {
			return nbrsOfSite[i];
		};
		const auto& getSitesOfEdge(long e) const {
			return sitesOfEdge[e];
		};

		void print() const;
		void computeWalls(std::vector<bool> &edgeOnWall,
						  const std::vector<int> &spins) const;
		std::vector<std::vector<long> > 
			computeLoops(const std::vector<bool> &edgeOnWall) const;
		std::vector<std::vector<long> > 
			computeLoops(const std::vector<bool> &edgeOnWall,
	  				     std::vector<long> &loopOfEdge) const;
		long computeWindingX(const std::vector<long> &loop) const;
		long computeWindingY(const std::vector<long> &loop) const;
		double computeGirationRadius(const std::vector<long> &loop) const;
		void countAll(std::vector<long> &count) const;
		void countSameLoop(const std::vector<long> &loopOfEdge,
						   std::vector<long> &count) const;

	private:
		const long nx, ny, ntot, n_edges;
		const double sx, sy; // Scaling units in x and y directions
		const long Lx, Ly; // Lengths in scaled units
		const double dr; // Step for correlations

		std::vector<long> xOfSite;
		std::vector<long> yOfSite;
		std::vector<long> xOfEdge;
		std::vector<long> yOfEdge;
		std::vector<std::vector<long> > nbrsOfSite;
		std::vector<std::vector<long> > edgesOfSite;
		std::vector<std::array<long, 2> > sitesOfEdge;
		std::vector<std::vector<long> > nbrsOfEdge;
		double maxDist;
		//std::vector<std::vector<double> > distsEdges;
		//std::vector<std::vector<size_t> > boxEdges;
		std::vector<size_t> boxEdges;
};

int testLattice();

#endif
