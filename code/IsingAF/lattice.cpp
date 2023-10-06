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
 * Implement the honeycomb lattice, and its dual the triangular lattice.
*/

#include <iostream>
#include <stack>
#include <cmath>
#include <cassert>
#include "lattice.h"

Lattice::Lattice(long nx_, long ny_, bool correl, double dr_) :
   		nx(nx_), ny(ny_), ntot(nx_ * ny_),
		n_edges(3 * nx_ * ny_), sx(sqrt(3.) / 4.), sy(0.25),
		Lx(2 * nx_), Ly(4 * ny_), dr(dr_) {
	assert(nx % 2 == 0);
	assert(ny % 2 == 0);

	xOfSite.resize(ntot);
	yOfSite.resize(ntot);
	nbrsOfSite.resize(ntot);
	edgesOfSite.resize(ntot);
	for (long i = 0 ; i < ntot ; ++i)
		edgesOfSite[i].resize(6);
	nbrsOfEdge.resize(n_edges);

	long x, y, xp, yp;

	for (long i = 0 ; i < ntot ; ++i) {
		x = i / ny;
		y = i % ny;

		// Positions
		xOfSite[i] = 2 * x;
		yOfSite[i] = 2 * (2 * y + (x % 2));

		// For the three first direction, we add edges
		// top
		yp = (y + 1) % ny;
		nbrsOfSite[i].push_back(x * ny + yp);
		sitesOfEdge.push_back({i, x * ny + yp});
		edgesOfSite[i][0] = sitesOfEdge.size() - 1;
		edgesOfSite[x * ny + yp][3] = sitesOfEdge.size() - 1;
		// top right
		xp = (x + 1) % nx;
		yp = (x % 2 == 0) ? y : ((y + 1) % ny);
		nbrsOfSite[i].push_back(xp * ny + yp);
		sitesOfEdge.push_back({i, xp * ny + yp});
		edgesOfSite[i][1] = sitesOfEdge.size() - 1;
		edgesOfSite[xp * ny + yp][4] = sitesOfEdge.size() - 1;
		// bottom right
		yp = (x % 2 == 0) ? ((y - 1 + ny) % ny) : y;
		nbrsOfSite[i].push_back(xp * ny + yp);
		sitesOfEdge.push_back({i, xp * ny + yp});
		edgesOfSite[i][2] = sitesOfEdge.size() - 1;
		edgesOfSite[xp * ny + yp][5] = sitesOfEdge.size() - 1;

		// bottom
		yp = (y - 1 + ny) % ny;
		nbrsOfSite[i].push_back(x * ny + yp);
		// bottom left
		xp = (x - 1 + nx) % nx;
		yp = (x % 2 == 0) ? ((y - 1 + ny) % ny) : y;
		nbrsOfSite[i].push_back(xp * ny + yp);
		// bottom left
		yp = (x % 2 == 0) ? y : ((y + 1) % ny);
		nbrsOfSite[i].push_back(xp * ny + yp);
	}
	assert(sitesOfEdge.size() == (size_t) n_edges);
	
	for (long i = 0 ; i < ntot ; ++i) {
		for (long k = 0 ; k < 6 ; ++k) {
			long e1 = edgesOfSite[i][k];
			long e2 = edgesOfSite[i][(k+1)%6];
			nbrsOfEdge[e1].push_back(e2);
			nbrsOfEdge[e2].push_back(e1);
		}
	}

	// Positions of edges are the middle of the segments
	long d, r;
	for (auto ss : sitesOfEdge) {
		d = moduloSym(xOfSite[ss[1]] - xOfSite[ss[0]], Lx);
		r = xOfSite[ss[0]] + d / 2;
		xOfEdge.push_back(modulo(r, Lx));
		d = moduloSym(yOfSite[ss[1]] - yOfSite[ss[0]], Ly);
		r = yOfSite[ss[0]] + d / 2;
		yOfEdge.push_back(modulo(r, Ly));
	}

	// Distances between edges
	if (correl) {
		// If the size of the array won't be too large
		if (nx <= NMAX_CORREL && ny <= NMAX_CORREL) {
			double dx, dy, d;
			maxDist = 0.;
			for (long i = 0 ; i < n_edges-1 ; ++i) {
				for (long j = 0 ; j < i ; ++j) {
					dx = sx * moduloSym(xOfEdge[i]-xOfEdge[j], Lx);
					dy = sy * moduloSym(yOfEdge[i]-yOfEdge[j], Ly);
					d = std::sqrt(dx*dx+dy*dy);
					boxEdges.push_back((size_t) (d / dr));
					if (d > maxDist)
						maxDist = d;
				}
			}
		} else {
			maxDist = 0.67 * std::max(nx, ny);
		}
		/*std::cout << "n_edges: " << n_edges << std::endl;
		std::cout << "Box edges: " << boxEdges.size() << std::endl;
		std::cout << "maxDist " << maxDist << std::endl;
		std::cout << "maxDist/dr " << maxDist/dr << std::endl;*/
	}
}

void Lattice::print() const {
	std::cout << "sx = " << sx << ", sy = " << sy << "\n";
	std::cout << "Lx = " << Lx << ", Ly = " << Ly << "\n\n";

	std::cout << "Sites (" << nbrsOfSite.size() << ")\n";
	for (long i = 0 ; i < ntot ; ++i) {
		std::cout << i << ": pos=(" << xOfSite[i] << ", " << yOfSite[i] << "), ";
		std::cout << "nbrs=";
		for (long j: nbrsOfSite[i])
			std::cout << j << " ";
		std::cout << " ; edges=";
		for (long j: edgesOfSite[i])
			std::cout << j << " ";
		std::cout << "\n";
	}
	std::cout << "\nEdges (" << sitesOfEdge.size() << ")\n";
	for (long e = 0 ; e < n_edges ; ++e) {
		std::cout << e << ": sites=[" << sitesOfEdge[e][0] << " " << sitesOfEdge[e][1];
		std::cout << "], pos=(" << xOfEdge[e] << ", " << yOfEdge[e] << "), "
			<< " ; nbrs=";
		for (long j: nbrsOfEdge[e])
			std::cout << j << " ";
		std::cout << "\n";
	}
}

void Lattice::computeWalls(std::vector<bool> &edgeOnWall,
			          	   const std::vector<int> &spins) const {
	edgeOnWall.resize(n_edges);
	for (long e = 0 ; e < n_edges ; ++e) {
		edgeOnWall[e] = (spins[sitesOfEdge[e][0]] != spins[sitesOfEdge[e][1]]);
	}
}

std::vector<std::vector<long> > 
		Lattice::computeLoops(const std::vector<bool> &edgeOnWall) const {
	std::vector<bool> visited(n_edges, false);
	std::vector<std::vector<long> > edgesOfLoop;

	// DFS
	for (long e = 0 ; e < n_edges ; ++e) {
		if (visited[e])
			continue;
		if (!edgeOnWall[e]) {
			visited[e] = true;
			continue;
		}

		// Otherwise create and explore new loop
		edgesOfLoop.push_back({});
		std::stack<long> S;
		S.push(e);
		while (!S.empty()) {
			long e1 = S.top();
			S.pop();
			if (visited[e1])
				continue;
			visited[e1] = true;
			if (!edgeOnWall[e1])
				continue;

			// Otherwise add edge to loop, and add neighbors to stack
			edgesOfLoop.back().push_back(e1);
			for (long e2 : nbrsOfEdge[e1]) 
				S.push(e2);
		}
	}

	return edgesOfLoop;
}

std::vector<std::vector<long> > 
		Lattice::computeLoops(const std::vector<bool> &edgeOnWall,
				              std::vector<long> &loopOfEdge) const {
	std::vector<bool> visited(n_edges, false);
	std::vector<std::vector<long> > edgesOfLoop;

	long iLoop = 0;
	for (long e = 0 ; e < n_edges ; ++e) {
		loopOfEdge[e] = -1;
	}

	// DFS
	for (long e = 0 ; e < n_edges ; ++e) {
		if (visited[e])
			continue;
		if (!edgeOnWall[e]) {
			visited[e] = true;
			continue;
		}

		// Otherwise create and explore new loop
		edgesOfLoop.push_back({});
		iLoop++;
		std::stack<long> S;
		S.push(e);
		while (!S.empty()) {
			long e1 = S.top();
			S.pop();
			if (visited[e1])
				continue;
			visited[e1] = true;
			if (!edgeOnWall[e1])
				continue;

			// Otherwise add edge to loop, and add neighbors to stack
			edgesOfLoop.back().push_back(e1);
			loopOfEdge[e1] = iLoop;
			for (long e2 : nbrsOfEdge[e1]) 
				S.push(e2);
		}
	}

	return edgesOfLoop;
}

long Lattice::computeWindingX(const std::vector<long> &loop) const {
	long x = 0;
	for (size_t e = 0 ; e < loop.size() - 1 ; ++e) {
		x += moduloSym(xOfEdge[loop[e+1]] - xOfEdge[loop[e]], Lx);
	}
	x += moduloSym(xOfEdge[loop[0]] - xOfEdge[loop.back()], Lx);
	return x / Lx;
}

long Lattice::computeWindingY(const std::vector<long> &loop) const {
	long y = 0;
	for (size_t e = 0 ; e < loop.size() - 1 ; ++e) {
		y += moduloSym(yOfEdge[loop[e+1]] - yOfEdge[loop[e]], Ly);
	}
	y += moduloSym(yOfEdge[loop[0]] - yOfEdge[loop.back()], Ly);
	return y / Ly;
}

double Lattice::computeGirationRadius(const std::vector<long> &loop) const {
	size_t N = loop.size();
	double Nf = (double) N;

	long x=0, y=0, sum_x=0, sum_y=0, sum_x2=0, sum_y2=0;
	for (size_t e = 0 ; e < N - 1 ; ++e) {
		x += moduloSym(xOfEdge[loop[e+1]] - xOfEdge[loop[e]], Lx);
		y += moduloSym(yOfEdge[loop[e+1]] - yOfEdge[loop[e]], Ly);
		sum_x += x;
		sum_y += y;
		sum_x2 += x * x;
		sum_y2 += y * y;
	}

	double xm = sum_x / Nf;
	double ym = sum_y / Nf;
	double x2v = (sum_x2 / Nf) - xm * xm;
	x2v *= sx * sx;
	double y2v = (sum_y2 / Nf) - ym * ym;
	y2v *= sy * sy;
	return sqrt(x2v + y2v);
}

void Lattice::countAll(std::vector<long> &count) const {
	// If we had enough memory to store the boxes
	if (!boxEdges.empty()) {
		for (size_t k = 0 ; k < boxEdges.size() ; ++k) {
			count[boxEdges[k]]++;
		}
	} else {
		double dx, dy, d;
		for (long i = 0 ; i < n_edges-1 ; ++i) {
			for (long j = 0 ; j < i ; ++j) {
				dx = sx * moduloSym(xOfEdge[i]-xOfEdge[j], Lx);
				dy = sy * moduloSym(yOfEdge[i]-yOfEdge[j], Ly);
				d = std::sqrt(dx*dx+dy*dy);
				count[(size_t) (d / dr)]++;
			}
		}
	}
}

void Lattice::countSameLoop(const std::vector<long> &loopOfEdge,
						    std::vector<long> &count) const {
	size_t k = 0;
	double dx, dy, d;

	for (long i = 0 ; i < n_edges-1 ; ++i) {
		for (long j = 0 ; j < i ; ++j) {
			++k;
			if (loopOfEdge[i] >= 0 && loopOfEdge[i] == loopOfEdge[j]) {
				// If we had enough memory to store the boxes
				if (!boxEdges.empty()) {
					count[boxEdges[k]]++;
				} else {
					dx = sx * moduloSym(xOfEdge[i]-xOfEdge[j], Lx);
					dy = sy * moduloSym(yOfEdge[i]-yOfEdge[j], Ly);
					d = std::sqrt(dx*dx+dy*dy);
					count[(size_t) (d / dr)]++;
				}
			}
		}
	}
}

int testLattice() {
	const long nx = 6, ny = 4;
	std::cout << "# nx = " << nx << ", ny = " << ny << "\n\n";

	Lattice lattice(nx, ny);
	lattice.print();

	std::vector<int> spins = {
		1, -1, -1, 1,
		1, 1, -1, -1,
		1, -1, 1, 1,
		1, 1, -1, -1,
		-1, -1, -1, 1,
		1, 1, 1, -1
	};

	std::vector<bool> edgeOnWall;
	lattice.computeWalls(edgeOnWall, spins);

	auto edgesOfLoop = lattice.computeLoops(edgeOnWall);
	std::cout << "\nLoops (" << edgesOfLoop.size() << ")\n";
	for (auto loop : edgesOfLoop) {
		std::cout << "Size " << loop.size() << ": ";
		for (long e : loop) 
			std::cout << "[" << lattice.getSitesOfEdge(e)[0] << " "
				<< lattice.getSitesOfEdge(e)[1] << "] ";
		long wx = lattice.computeWindingX(loop);
		long wy = lattice.computeWindingY(loop);
		std::cout << "; windings=" << wx << " " << wy;
		double Rg = lattice.computeGirationRadius(loop);
		std::cout << " ; Rg=" << Rg << "\n";
		/*double xG = lattice.computeCenterOfMassX(loop);
		double yG = lattice.computeCenterOfMassY(loop);
		std::cout << " ; CM=" << xG << " " << yG << "\n";*/
	}

	return 0;

	/*std::cout << "\nEdge on wall:\n";
	for (bool b : edgeOnWall) {
		std::cout << b << " ";
	}
	std::cout << "\n\n";*/
}
