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
 * lattice.cpp
 *
 * Author: Alexis Poncet
 * Email: alexis.poncet@ens-lyon.fr
 *
 * Implement the honeycomb lattice.
*/

#include <iostream>
#include <fstream>
#include <stack>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <random>
#include <cmath>
#include <cassert>
#include "lattice.h"

Lattice::Lattice(long nx_, long ny_, bool periodic_, bool correl, double dr_) :
		nx(nx_), ny(ny_), periodic(periodic_), sx(0.25), sy(sqrt(3.) / 4.),
		Lx(6 * nx_), Ly(4 * ny_), dr(dr_)  {
	assert(nx % 2 == 0);

	if (periodic) {
		makePeriodicLattice();
	} else {
		makeNonPeriodicLattice();
	}
	
	// Distances between edges
	if (correl) {
		double dx, dy, d;
		maxDist = 0.;
		for (long i = 0 ; i < n_edges-1 ; ++i) {
			for (long j = 0 ; j < i ; ++j) {
				if (periodic) {
					dx = sx * moduloSym(xOfEdge[i]-xOfEdge[j], Lx);
					dy = sy * moduloSym(yOfEdge[i]-yOfEdge[j], Ly);
				} else {
					dx = sx * (xOfEdge[i]-xOfEdge[j]);
					dy = sy * (yOfEdge[i]-yOfEdge[j]);
				}
				d = std::sqrt(dx*dx+dy*dy);
				boxEdges.push_back((size_t) (d / dr));
				if (d > maxDist)
					maxDist = d;
			}
		}
	}
}

void Lattice::makePeriodicLattice() {
	ntot = 2 * nx * ny;
	n_edges = 3 * nx * ny;

	xOfSite.resize(ntot);
	yOfSite.resize(ntot);
	nbrsOfSite.resize(ntot);
	edgesOfSite.resize(ntot);
	for (long i = 0 ; i < ntot ; ++i)
		edgesOfSite[i].resize(3);
	nbrsOfEdge.resize(n_edges);

	long i, xp, yp, ip;
	for (long x = 0 ; x < nx ; ++x) {
		for (long y = 0 ; y < ny ; ++y) {
			i = (2 * y) + (2 * ny * x);

			for (long a = 0 ; a < 2 ; ++a) {
				xOfSite[i+a] = 6 * x + 4 * a;
				yOfSite[i+a] = 4 * y + 2 * (1 - (x % 2));
			}
			
			// Internal edge
			nbrsOfSite[i].push_back(i+1);
			nbrsOfSite[i+1].push_back(i);
			sitesOfEdge.push_back({i, i+1}); // {a=0, a=1}
			edgesOfSite[i][0] = sitesOfEdge.size() - 1;
			edgesOfSite[i+1][0] = sitesOfEdge.size() - 1;
			// Bottom right edge
			xp = (x + 1) % nx;
			yp = (y - (x % 2) + ny) % ny;
			ip = (2 * yp) + (2 * ny * xp);
			nbrsOfSite[i+1].push_back(ip);
			nbrsOfSite[ip].push_back(i+1);
			sitesOfEdge.push_back({ip, i+1}); // {a=0, a=1}
			edgesOfSite[i+1][1] = sitesOfEdge.size() - 1;
			edgesOfSite[ip][1] = sitesOfEdge.size() - 1;
			// Top right edge
			yp = (y + 1 - (x % 2)) % ny;
			ip = (2 * yp) + (2 * ny * xp);
			nbrsOfSite[i+1].push_back(ip);
			nbrsOfSite[ip].push_back(i+1);
			sitesOfEdge.push_back({ip, i+1}); // {a=0, a=1}
			edgesOfSite[i+1][2] = sitesOfEdge.size() - 1;
			edgesOfSite[ip][2] = sitesOfEdge.size() - 1;
		}
	}
	assert(sitesOfEdge.size() == (size_t) n_edges);

	// We add the neighbors at sites (a=0) before those at (a=1)
	// The neighbors of an edge end up being listed in circular order
	nbrsOfEdge.resize(n_edges);
	for (int a = 0 ; a < 2 ; ++a) {
		for (long i = 0 ; i < ntot ; i+=2) {
			for (long k = 0 ; k < 3 ; ++k) {
				long e1 = edgesOfSite[i+a][k];
				long e2 = edgesOfSite[i+a][(k+1)%3];
				nbrsOfEdge[e1].push_back(e2);
				nbrsOfEdge[e2].push_back(e1);
			}
		}
	}

	edgeInBulk.assign(n_edges, true);

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
}

void Lattice::makeNonPeriodicLattice() {
	ntot = 2 * nx * (ny + 1) + 2 * ny;
	n_faces = nx * ny;

	xOfSite.resize(ntot);
	yOfSite.resize(ntot);
	nbrsOfSite.resize(ntot);
	edgesOfSite.resize(ntot);
	for (long i = 0 ; i < ntot ; ++i)
		edgesOfSite[i].assign(3, -1);
	edgesOfFace.resize(n_faces);
	nbrsOfFace.resize(n_faces);
	// Does a spin 1 on the edge goes in the ccc or cc direction
	dirNbrOfFace.resize(n_faces);

	long i, f, f2, f3, e, xp, yp, ip;
	bool f_exists, f2_exists, f3_exists;
	long ny2 = ny + 1;

	// Bulk
	for (long x = 0 ; x < nx ; ++x) {
		for (long y = 0 ; y < ny2 ; ++y) {
			i = (2 * y) + (2 * ny2 * x);
			// Neighboring faces
			f = y + ny * x;
			f_exists = (y < ny);
			f2 = (y - 1) + ny * x;
			f2_exists = (y > 0);
			f3 = y - (x % 2) + ny * (x + 1);
			f3_exists = (y - (x % 2) >= 0) & (y - (x % 2) < ny) & (x < nx - 1);

			for (long a = 0 ; a < 2 ; ++a) {
				xOfSite[i+a] = 6 * x + 4 * a;
				yOfSite[i+a] = 4 * y + 2 * (1 - (x % 2));
			}
			
			// Internal edge (always there)
			nbrsOfSite[i].push_back(i+1);
			nbrsOfSite[i+1].push_back(i);
			sitesOfEdge.push_back({i, i+1}); // {a=0, a=1}
			e = sitesOfEdge.size() - 1;
			edgesOfSite[i][0] = e;
			edgesOfSite[i+1][0] = e;
			// f-f2 edge
			if (f_exists && f2_exists) {
				edgesOfFace[f].push_back(e);
				edgesOfFace[f2].push_back(e);
				nbrsOfFace[f].push_back(f2);
				nbrsOfFace[f2].push_back(f);
				dirNbrOfFace[f].push_back(+1);
				dirNbrOfFace[f2].push_back(-1);
			}

			xp = x + 1;
			if (xp < nx) {
				// Bottom right edge
				yp = y - (x % 2);
				if (yp >= 0) {
					ip = (2 * yp) + (2 * ny2 * xp);
					nbrsOfSite[i+1].push_back(ip);
					nbrsOfSite[ip].push_back(i+1);
					sitesOfEdge.push_back({ip, i+1}); // {a=0, a=1}
					e = sitesOfEdge.size() - 1;
					edgesOfSite[i+1][1] = e;
					edgesOfSite[ip][1] = e;
					// f2-f3 edge
					if (f2_exists && f3_exists) {
						edgesOfFace[f2].push_back(e);
						edgesOfFace[f3].push_back(e);
						nbrsOfFace[f2].push_back(f3);
						nbrsOfFace[f3].push_back(f2);
						dirNbrOfFace[f2].push_back(+1);
						dirNbrOfFace[f3].push_back(-1);
					}
				}

				// Top right edge
				yp = y + 1 - (x % 2);
				if (yp < ny2) {
					ip = (2 * yp) + (2 * ny2 * xp);
					nbrsOfSite[i+1].push_back(ip);
					nbrsOfSite[ip].push_back(i+1);
					sitesOfEdge.push_back({ip, i+1}); // {a=0, a=1}
					e = sitesOfEdge.size() - 1;
					edgesOfSite[i+1][2] = e;
					edgesOfSite[ip][2] = e;

					// f-f3 edge
					if (f_exists && f3_exists) {
						edgesOfFace[f].push_back(e);
						edgesOfFace[f3].push_back(e);
						nbrsOfFace[f].push_back(f3);
						nbrsOfFace[f3].push_back(f);
						dirNbrOfFace[f].push_back(-1);
						dirNbrOfFace[f3].push_back(+1);
					}
				}
			}
		}
	}

	// Left and right boundaries
	i = 2 * nx * ny2;
	for (long y = 0 ; y < ny ; ++y) {
		// RIGHT
		//i = 2 * nx * ny2 + ny + y;
		xOfSite[i] = 6 * nx;
		yOfSite[i] = 4 * y + 2 * (1 - (nx % 2));

		// Bottom left edge
		ip = 2 * (nx - 1) * ny2 + 1 + 2 * y;
		nbrsOfSite[i].push_back(ip);
		nbrsOfSite[ip].push_back(i);
		sitesOfEdge.push_back({i, ip}); // {a=0, a=1}
		edgesOfSite[i][2] = sitesOfEdge.size() - 1;
		edgesOfSite[ip][2] = sitesOfEdge.size() - 1;

		// Top right edge
		ip = 2 * (nx - 1) * ny2 + 3 + 2 * y;
		nbrsOfSite[i].push_back(ip);
		nbrsOfSite[ip].push_back(i);
		sitesOfEdge.push_back({i, ip}); // {a=0, a=1}
		edgesOfSite[i][1] = sitesOfEdge.size() - 1;
		edgesOfSite[ip][1] = sitesOfEdge.size() - 1;
		++i;

		// LEFT
		//i = 2 * nx * ny2 + y;
		xOfSite[i] = -2;
		yOfSite[i] = 4 * (y + 1);

		// Bottom right edge
		ip = 2 * y;
		nbrsOfSite[i].push_back(ip);
		nbrsOfSite[ip].push_back(i);
		sitesOfEdge.push_back({ip, i}); // {a=0, a=1}
		edgesOfSite[i][1] = sitesOfEdge.size() - 1;
		edgesOfSite[ip][1] = sitesOfEdge.size() - 1;

		// Top right edge
		ip = 2 * (y + 1);
		nbrsOfSite[i].push_back(ip);
		nbrsOfSite[ip].push_back(i);
		sitesOfEdge.push_back({ip, i}); // {a=0, a=1}
		edgesOfSite[i][2] = sitesOfEdge.size() - 1;
		edgesOfSite[ip][2] = sitesOfEdge.size() - 1;
		++i;
	}

	n_edges = sitesOfEdge.size();

	// We add the neighbors at sites (a=0) before those at (a=1)
	// The neighbors of an edge end up being listed in circular order
	nbrsOfEdge.resize(n_edges);
	edgeInBulk.assign(n_edges, true);
	for (int a = 0 ; a < 2 ; ++a) {
		for (long i = 0 ; i < ntot ; i+=2) {
			for (long k = 0 ; k < 3 ; ++k) {
				long e1 = edgesOfSite[i+a][k];
				long e2 = edgesOfSite[i+a][(k+1)%3];
				if (e1 >= 0)
					nbrsOfEdge[e1].push_back(e2);
				else
					edgeInBulk[e2] = false;
				if (e2 >= 0)
					nbrsOfEdge[e2].push_back(e1);
				else
					edgeInBulk[e1] = false;
			}
		}
	}

	// Remove (-1)s in edgesOfSite
	for (long i = 0 ; i < ntot ; ++i) {
		edgesOfSite[i].erase(
				std::remove(edgesOfSite[i].begin(), edgesOfSite[i].end(), -1),
				edgesOfSite[i].end()
		);
	}

	// Positions of edges are the middle of the segments
	for (auto ss : sitesOfEdge) {
		xOfEdge.push_back((xOfSite[ss[1]] + xOfSite[ss[0]]) / 2);
		yOfEdge.push_back((yOfSite[ss[1]] + yOfSite[ss[0]]) / 2);
	}

	// Create loop through the faces
	std::vector<long> of; // ordered faces
	for (long x = 0 ; x < nx ; ++x) {
		if (x % 2 == 0)
			for (long y = 0 ; y < ny ; ++y)
				of.push_back(y + ny * x);
		else
			for (long y = ny-1 ; y >= 0 ; --y)
				of.push_back(y + ny * x);
	}
	for (long f = 0 ; f < n_faces - 1 ; ++f) {
		size_t k = std::distance(
			nbrsOfFace[of[f]].begin(),
			std::find(nbrsOfFace[of[f]].begin(), nbrsOfFace[of[f]].end(),
					  of[f+1])
			);
		edgesThroughFaces.push_back(edgesOfFace[of[f]][k]);
		dirThroughFaces.push_back(dirNbrOfFace[of[f]][k]);
	}
}

void Lattice::print() const {
	std::cout << "sx = " << sx << ", sy = " << sy << "\n";
	std::cout << "Lx = " << Lx << ", Ly = " << Ly << "\n\n";

	std::cout << "Sites (" << nbrsOfSite.size() << ")\n";
	for (long i = 0 ; i < ntot ; ++i) {
		std::cout << i << ": pos=(" << xOfSite[i] << ", "
			<< yOfSite[i] << "), ";
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
		std::cout << e << (edgeInBulk[e] ? " (bulk)" : " (bord)")
			<< ": sites=[" << sitesOfEdge[e][0] << " "
			<< sitesOfEdge[e][1];
		std::cout << "], pos=(" << xOfEdge[e] << ", " << yOfEdge[e] << "), "
			<< " nbrs=";
		for (long j: nbrsOfEdge[e]) {
			if (j >= 0) {
				std::cout << "[" << sitesOfEdge[j][0] << " "
					<< sitesOfEdge[j][1] << "] ";
			}
			//std::cout << j << " ";
		}
		std::cout << "\n";
	}
	
	if (!periodic) {
		std::cout << "\nFaces (" << edgesOfFace.size() << ")\n";
		for (long f = 0 ; f < n_faces ; ++f) {
			std::cout << f << ": ";
			std::cout << "nbrs=";
			for (long j: nbrsOfFace[f])
				std::cout << j << " ";
			std::cout << " ; edges=";
			for (long j: edgesOfFace[f])
				std::cout << j << " ";
			std::cout << " ; dirs=";
			for (int j: dirNbrOfFace[f])
				std::cout << j << " ";
			std::cout << "\n";
		}
	}

	std::cout << "\nLoop though faces:\nedges ";
	for (long i: edgesThroughFaces)
		std::cout << i << " ";
	std::cout << "\ndirs: ";
	for (int i: dirThroughFaces)
		std::cout << i << " ";
	std::cout << "\n";
}

void Lattice::output(const std::string ofname) const {
	std::ofstream file;
	file.open(ofname + "_pos.dat");
	for (long i = 0 ; i < ntot ; ++i) {
		file << xOfSite[i] << " " << yOfSite[i] << "\n";
	}
	file.close();

	file.open(ofname + "_edges.dat");
	std::cout << n_edges << " edges\n";
	for (long e = 0 ; e < n_edges ; ++e) {
		file << sitesOfEdge[e][0] << " " << sitesOfEdge[e][1] << "\n";
	}
	file.close();

	if (!periodic) {
		file.open(ofname + "_faces.dat");
		for (long f = 0 ; f < n_faces ; ++f) {
			for (long j: nbrsOfFace[f])
				file << j << " ";
			file << "\n";
		}
		file.close();
	}
}

int Lattice::checkIceRule(const std::vector<int> &spins) const {
	for (long i = 0 ; i < ntot ; ++i) {
		int s = 0;
		for (long e : edgesOfSite[i]) 
			s += spins[e];
		if (s != 0) {
			std::cout << i << "\n";
			return 0;
		}
	}
	return 1;
}

std::vector<std::vector<long> > 
		Lattice::computeLoops(const std::vector<int> &spins) const {
	std::vector<bool> visited(n_edges, false);
	std::vector<std::vector<long> > edgesOfLoop;

	// DFS
	for (long e = 0 ; e < n_edges ; ++e) {
		if (visited[e])
			continue;
		if (!spins[e]) {
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
			if (!spins[e1])
				continue;

			// Otherwise add edge to loop, and add neighbors to stack
			edgesOfLoop.back().push_back(e1);
			for (long e2 : nbrsOfEdge[e1]) 
				if (e2 >= 0)
					S.push(e2);
		}
	}

	return edgesOfLoop;
}

std::vector<std::vector<long> > 
		Lattice::computeLoops(const std::vector<int> &spins,
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
		if (!spins[e]) {
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
			if (!spins[e1])
				continue;

			// Otherwise add edge to loop, and add neighbors to stack
			edgesOfLoop.back().push_back(e1);
			loopOfEdge[e1] = iLoop;
			for (long e2 : nbrsOfEdge[e1]) 
				if (e2 >= 0)
					S.push(e2);
		}
	}

	return edgesOfLoop;
}

long Lattice::computeWindingX(const std::vector<long> &loop) const {
	// Only for periodic boundaries
	long x = 0;
	for (size_t e = 0 ; e < loop.size() - 1 ; ++e) {
		x += moduloSym(xOfEdge[loop[e+1]] - xOfEdge[loop[e]], Lx);
	}
	x += moduloSym(xOfEdge[loop[0]] - xOfEdge[loop.back()], Lx);
	return x / Lx;
}

long Lattice::computeWindingY(const std::vector<long> &loop) const {
	// Only for periodic boundaries
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
		if (periodic) {
			x += moduloSym(xOfEdge[loop[e+1]] - xOfEdge[loop[e]], Lx);
			y += moduloSym(yOfEdge[loop[e+1]] - yOfEdge[loop[e]], Ly);
		} else {
			x += (xOfEdge[loop[e+1]] - xOfEdge[loop[e]]);
			y += (yOfEdge[loop[e+1]] - yOfEdge[loop[e]]);
		}
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

long Lattice::computeHeight(const std::vector<int>& spins) const {
	// We previously found a loop through the faces
	// we follow it to assign the heights and record the min and max.
	long h = 0, min_h = 0, max_h = 0;

	for (long k = 0 ; k < n_faces - 1; ++k) {
		h += dirThroughFaces[k] * spins[edgesThroughFaces[k]];
		if (h < min_h)
			min_h = h;
		else if (h > max_h)
			max_h = h;
	}
	return max_h - min_h;
}

/*
long Lattice::makeWorm(const std::vector<int> &spins,
					   std::deque<long> &edgesVisited,
					   std::deque<int> &newSpins,
		               VSLStreamStatePtr stream) const {
	// We assume that edgesVisited and newSpins are empty
	int e, sNew, a;
	long i;

#ifndef USE_MAP
	// Visited sites implemented either with a vector (or deque)
	std::vector<long> sitesVisited;
#else
	// Or with a unordered_map (or map); key: site visited, value: time of visit
	std::unordered_map<long, long> sitesVisited;
#endif
	// It turns out that the vector solution is at least as fast as
	// the unordered_map probably because average loop size = constant
	// (wrt system size) and smaller overhead

	unsigned long long rndBits; // random bits
	unsigned bit = 0;
	viRngUniformBits64(VSL_RNG_METHOD_UNIFORMBITS64_STD, stream, 1, &rndBits);

	// Choose edge at random
	viRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, 1, &e, 0, n_edges);
	edgesVisited.push_back(e);

	// New state
	if (spins[e] != 0) {
		sNew = 0;
	} else {
		//viRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, 1, &sNew, 0, 2);
		//sNew = 2 * sNew - 1;
		sNew = 2 * (rndBits & 1) - 1;
		rndBits >>= 1;
		bit++;
	}
	newSpins.push_back(sNew);

	// Expand the worm once
	//viRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, 1, &a, 0, 2);
	a = rndBits & 1;
	rndBits >>= 1;
	bit++;
	i = sitesOfEdge[e][a]; // Chosen site
#ifndef USE_MAP
	sitesVisited.push_back(sitesOfEdge[e][1-a]); // Other site ; O(1)
#else
	long t = 0;
	sitesVisited[sitesOfEdge[e][1-a]] = t; // Other site ; O(log(N))
#endif

	// Expand the worm until comeback to visited site
	bool explore = true;
	auto search = sitesVisited.end();
	do {
#ifndef USE_MAP
		sitesVisited.push_back(i); // O(1)
#else
		t++;
		sitesVisited[i] = t; // O(log(N))
#endif

		// The two neighbors of the edge on side a in random order
		//vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, 1, &u, 0., 1.);
		long e1 = nbrsOfEdge[e][2*a]; 
		long e2 = nbrsOfEdge[e][2*a+1]; 
		if (rndBits & 1)
			std::swap(e1, e2);
		rndBits >>= 1;
		bit++;

		// New random bits if needed
		if (bit == N_BITS_64) {
			bit = 0;
			viRngUniformBits64(VSL_RNG_METHOD_UNIFORMBITS64_STD, stream,
					           1, &rndBits);
		}
		
		// If we make a defect with a zero edge, we can repair it either
		// by zeroing the only non-zero edge at the site
		// or by setting the other zero edge to a non-zero value
		if (sNew == 0) {
			e = e1; // New edge
			sNew = -spins[e2]; // New spin
		// If we previously had a 3-zero site, we repair it by assigning
		// the correct value to one of the two remaining 0s.
		} else if (spins[e1] == 0) {
			e = e1;
			sNew = -sNew;
		// The last case corresponds to a new +/-1 spin coming into
		// a vertex with a +1 and a -1 spin. Only one way to repair it.
		} else {
			e = (spins[e1] == sNew) ? e1 : e2;
			sNew = 0;
		}

		edgesVisited.push_back(e);
		newSpins.push_back(sNew);
		a = 1 - a; // Other side of edge
		i = sitesOfEdge[e][a]; // New site

#ifndef USE_MAP
		search = std::find(sitesVisited.begin(), sitesVisited.end(), i); // O(N)
		explore = (search == sitesVisited.end());
#else
		search = sitesVisited.find(i); // O(log(N))
		explore = (search == sitesVisited.end());
#endif
	} while(explore);

	 // Return time at which final site was first found
#ifndef USE_MAP
	return search - sitesVisited.begin();
#else
	return search->second;
#endif
}
*/

long Lattice::makeWorm(const std::vector<int> &spins,
					   std::deque<long> &edgesVisited,
					   std::deque<int> &newSpins,
		               VSLStreamStatePtr stream) const {
	// We assume that edgesVisited and newSpins are empty
	int e, sNew, a;
	long i;

	// Visited sites implemented with a vector
	std::vector<long> sitesVisited;

	unsigned long long rndBits; // random bits
	unsigned bit = 0;
	viRngUniformBits64(VSL_RNG_METHOD_UNIFORMBITS64_STD, stream, 1, &rndBits);

	// Choose edge at random
	viRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, 1, &e, 0, n_edges);
	edgesVisited.push_back(e);

	// New state
	if (spins[e] != 0) {
		sNew = 0;
	} else {
		sNew = 2 * (rndBits & 1) - 1;
		rndBits >>= 1;
		bit++;
	}
	newSpins.push_back(sNew);

	// Expand the worm once
	a = rndBits & 1;
	rndBits >>= 1;
	bit++;
	i = sitesOfEdge[e][a]; // Chosen site
	sitesVisited.push_back(sitesOfEdge[e][1-a]); // Other site ; O(1)

	// Expand the worm until comeback to visited site
	bool explore = true;
	auto search = sitesVisited.end();
	do {
		sitesVisited.push_back(i); // O(1)

		// The two neighbors of the edge on side a
		long e1 = nbrsOfEdge[e][2*a]; 
		long e2 = nbrsOfEdge[e][2*a+1]; 

		// Rules on border (only one choice)
		if (e2 < 0) {
			e = e1; // New edge
			sNew = -sNew; // To satisfy ice rule with 2 edges
		} else if (e1 < 0) {
			e = e2; // New edge
			sNew = -sNew; // To satisfy ice rule with 2 edges
		} else { // Rules in bulk
			// e1 and e2 in random order
			if (rndBits & 1)
				std::swap(e1, e2);
			rndBits >>= 1;
			bit++;

			// New random bits if needed
			if (bit == N_BITS_64) {
				bit = 0;
				viRngUniformBits64(VSL_RNG_METHOD_UNIFORMBITS64_STD, stream,
								   1, &rndBits);
			}
			
			// If we make a defect with a zero edge, we can repair it either
			// by zeroing the only non-zero edge at the site
			// or by setting the other zero edge to a non-zero value
			if (sNew == 0) {
				e = e1; // New edge
				sNew = -spins[e2]; // New spin
			// If we previously had a 3-zero site, we repair it by assigning
			// the correct value to one of the two remaining 0s.
			} else if (spins[e1] == 0) {
				e = e1;
				sNew = -sNew;
			// The last case corresponds to a new +/-1 spin coming into
			// a vertex with a +1 and a -1 spin. Only one way to repair it.
			} else {
				e = (spins[e1] == sNew) ? e1 : e2;
				sNew = 0;
			}
		}

		edgesVisited.push_back(e);
		newSpins.push_back(sNew);
		a = 1 - a; // Other side of edge
		i = sitesOfEdge[e][a]; // New site

		search = std::find(sitesVisited.begin(), sitesVisited.end(), i); // O(N)
		explore = (search == sitesVisited.end());
	} while(explore);

	 // Return time at which final site was first found
	return search - sitesVisited.begin();
}

/*long Lattice::couplingEnergy(const std::vector<int> &spins) const {
	long E = 0;
	for (long e = 0 ; e < n_edges ; ++e)
		if(spins[e] == 0)
			E += couplingEnergyEdge(e, spins);
	return E;
}*/

long Lattice::couplingEnergyEdge(const long e, const std::vector<int> &spins)
		const {
	if (edgeInBulk[e]) {
		return (spins[nbrsOfEdge[e][0]] - spins[nbrsOfEdge[e][1]])
			 * (spins[nbrsOfEdge[e][2]] - spins[nbrsOfEdge[e][3]]);
	} else {
		return 0;
	}
}

long Lattice::diffCouplingEnergy(const std::vector<int> &spins1,
		                         const std::vector<int> &spins2,
								 const std::deque<long> &edgesVisited,
								 const long start) const {
	std::vector<bool> visited(n_edges, false);
	long E = 0;

	for (size_t i = start ; i < edgesVisited.size() ; ++i) {
		long e = edgesVisited[i];
		visited[e] = true;
		if (spins2[e] == 0) {
			E += couplingEnergyEdge(e, spins2);
		} else if (spins1[e] == 0) {
			E -= couplingEnergyEdge(e, spins1);
		}
		
		for (long e2 : nbrsOfEdge[e]) {
			if (e2 >= 0 && !visited[e2]) {
				visited[e2] = true;
				if (spins2[e2] == 0) {
					E += couplingEnergyEdge(e2, spins2);
				}
				if (spins1[e2] == 0) {
					E -= couplingEnergyEdge(e2, spins1);
				}
			}
		}
	}

	return E;
}

long Lattice::countZeroSites(const std::vector<int> &spins) const {
	long c = 0, a;
	for (long i = 0 ; i < ntot ; ++i) {
		a = 1;
		for (long e : edgesOfSite[i]) 
			if (spins[e] != 0)
				a = 0;
		c += a;
	}
	return c;
}

double Lattice::computeFracPara(const std::vector<int> &spins) const {
	long para = 0, tot = 0, c;

	for (long e = 0 ; e < n_edges ; ++e) {
		if (spins[e] == 0) {
			c = couplingEnergyEdge(e, spins);
			if (c == 4) {
				para++;
				tot++;
			} else if (c == -4) {
				tot++;
			}
		}
	}
	return para / ((double) tot);
}

void Lattice::countAll(std::vector<long> &count) const {
	for (size_t k = 0 ; k < boxEdges.size() ; ++k) {
		count[boxEdges[k]]++;
	}
}

void Lattice::countSameLoop(const std::vector<long> &loopOfEdge,
						    std::vector<long> &count) const {
	size_t k = 0;
	for (long i = 0 ; i < n_edges-1 ; ++i) {
		for (long j = 0 ; j < i ; ++j) {
			++k;
			if (loopOfEdge[i] >= 0) {
				if (loopOfEdge[i] == loopOfEdge[j]) {
					count[boxEdges[k]]++;
				}
			}
		}
	}
}

int testLattice() {
	//const long nx = 6, ny = 4;
	const long nx = 4, ny = 4;
	std::cout << "# Periodic: nx = " << nx << ", ny = " << ny << "\n\n";

	Lattice lattice(nx, ny);
	lattice.print();

	std::vector<int> spins = {
		0, 0, 0, 0, -1, 1,
		0, 0, 0, 1, -1, 0,
		0, -1, 1, 1, 0, -1,
		-1, 1, 0, 1, -1, 0,
		-1, 0, 1, 0, -1, 1,
		1, 0, -1, 1, 0, -1,
		1, 0, -1, 0, 1, -1,
		-1, 1, 0, 1, 0, -1
	};

	if (lattice.checkIceRule(spins))
		std::cout << "\nIce rule is satisfied!\n";
	else
		std::cout << "\nIce rule is NOT satisfied!!??\n";

	{
	auto edgesOfLoop = lattice.computeLoops(spins);
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
		std::cout << " ; Rg=" << Rg;
		std::cout << "\n";
	}
	}

	// Worm
	VSLStreamStatePtr stream;
	std::random_device rd;
	auto seed = rd();
	vslNewStream(&stream, VSL_BRNG_SFMT19937, seed);

	std::deque<long> edgesVisited;
	std::deque<int> newSpins;
	long n = lattice.makeWorm(spins, edgesVisited, newSpins, stream);

	std::cout << "\n Worm (seed=" << seed << "):\n";

	for (size_t i = n ; i < edgesVisited.size() ; ++i) {
		long e = edgesVisited[i];
		double sOld = spins[e];
		double sNew = newSpins[i];
		std::cout << "[" << lattice.getSitesOfEdge(e)[0] << " "
			<< lattice.getSitesOfEdge(e)[1] << "]: "
			<< sOld << " -> " << sNew << "\n";
		spins[e] = sNew;
	}
	if (lattice.checkIceRule(spins))
		std::cout << "\nIce rule is satisfied!\n";
	else
		std::cout << "\nIce rule is NOT satisfied!!??\n";

	{
	auto edgesOfLoop = lattice.computeLoops(spins);
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
		std::cout << " ; Rg=" << Rg;
		std::cout << "\n";
	}
	}

	return 0;
}

int testLatticeNotPeriodic() {
	const long nx = 4, ny = 3;
	std::cout << "\n\n# Non periodic: nx = " << nx << ", ny = " << ny << "\n\n";

	Lattice lattice(nx, ny, false);
	lattice.print();

	std::vector<int> spins = {
		-1, 0, 1, 1, 0,
		-1, 1, 0, -1, -1,
		1, 0, 0, -1, 0,
		1, 1, -1, 0, 0,
		1, -1, 0, 1, -1,
		0, 1, -1, -1, 1,
		0, 1, -1, -1, 0,
		0, 1, 1, -1, 1,
		-1, 1, -1, 0, 0,
		1, -1, -1, 1
	};

	if (lattice.checkIceRule(spins))
		std::cout << "\nIce rule is satisfied!\n";
	else
		std::cout << "\nIce rule is NOT satisfied!!??\n";

	{
	auto edgesOfLoop = lattice.computeLoops(spins);
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
		std::cout << " ; Rg=" << Rg;
		std::cout << "\n";
	}
	}

	// Worm
	VSLStreamStatePtr stream;
	std::random_device rd;
	auto seed = rd();
	vslNewStream(&stream, VSL_BRNG_SFMT19937, seed);

	std::deque<long> edgesVisited;
	std::deque<int> newSpins;
	long n = lattice.makeWorm(spins, edgesVisited, newSpins, stream);

	std::cout << "\n Worm (seed=" << seed << "):\n";

	for (size_t i = n ; i < edgesVisited.size() ; ++i) {
		long e = edgesVisited[i];
		double sOld = spins[e];
		double sNew = newSpins[i];
		std::cout << "[" << lattice.getSitesOfEdge(e)[0] << " "
			<< lattice.getSitesOfEdge(e)[1] << "]: "
			<< sOld << " -> " << sNew << "\n";
		spins[e] = sNew;
	}
	if (lattice.checkIceRule(spins))
		std::cout << "\nIce rule is satisfied!\n";
	else
		std::cout << "\nIce rule is NOT satisfied!!??\n";

	{
	auto edgesOfLoop = lattice.computeLoops(spins);
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
		std::cout << " ; Rg=" << Rg;
		std::cout << "\n";
	}
	}

	return 0;
}
