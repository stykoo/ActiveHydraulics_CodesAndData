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
 * parameters.cpp
 *
 * Author: Alexis Poncet
 * Email: alexis.poncet@ens-lyon.fr
 *
 * Parse, check and print command-line arguments.
*/

#include <exception>
#include <boost/program_options.hpp>
#include "parameters.h"

namespace po = boost::program_options;
#include "parameters.h"

// Check if the parameters are valid. Return 0 if they are, 1 otherwise.
int Parameters::check() const {
	// TODO
	return 0;
}

// Print the parameters to stream.
void Parameters::print(std::ostream &stream) const {
	stream << "# " << (periodic ? "[periodic] " : "[non-periodic] ")
		<< "nx=" << nx << ", ny=" << ny << ", K=" << K << ", J=" << J
		<< ", nbIters=" << nbIters << ", nbItersTh=" << nbItersTh
		<< ", nbSimuls=" << nbSimuls << ", skip=" << skip << "\n";
}

int Parameters::fromCommandLine(int argc, char **argv) {
	po::options_description opts("Options");
	opts.add_options()
		("nx", po::value<long>(&nx)->required(), "Number of sites in x direction")
		("ny", po::value<long>(&ny)->required(), "Number of sites in y direction")
		("simuls,s", po::value<long>(&nbSimuls)->default_value(1),
		 "Number simulations")
		("nbIters", po::value<long>(&nbIters)->required(),
		 "Number of iterations")
		("nbItersTh", po::value<long>(&nbItersTh)->default_value(0),
		 "Number of iterations of thermalization")
		("skip", po::value<long>(&skip)->default_value(-1),
		 "Number of iterations before computation (default: nx*ny)")
		("K", po::value<double>(&K)->required(), "Coupling constant for flow")
		("J", po::value<double>(&J)->default_value(0.),
		 "Coupling constant for orientation")
		("output,o", po::value<std::string>(&output)->default_value(
			DEFAULT_OUTPUT_FILE), "Output file")
        ("per,p", po::bool_switch(&periodic), "Use periodic boundary conditions")
        ("spins", po::bool_switch(&outputSpins), "Output the spins")
        ("obs", po::bool_switch(&outputObs), "Output the observables")
        ("correl,c", po::bool_switch(&computeCorrel),
		 "Compute correlations in loop")
		("dr", po::value<double>(&dr)->default_value(.25),
		 "Step for correlations")
        ("test", po::bool_switch(&test), "Test mode")
        ("lattice", po::bool_switch(&outputLattice), "Output lattice")
        ("verbose,v", po::bool_switch(&verbose), "Verbose mode")
        ("help,h", "Print help message and exit")
		;

	try {
		po::variables_map vars;
		po::store(po::parse_command_line(argc, argv, opts), vars);

        // Display help and exit
        if (vars.count("help")) {
			std::cout << "Usage: " << argv[0] << " options\n";
			std::cout << opts << std::endl;
            return 1;
        }

        po::notify(vars);

		//ntot = 2 * nx * ny; // Total number of sites
		//n_edges = 3 * nx * ny; // Total number of sites
		if (skip <= 0) 
			skip = nx * ny;

		if ((nx % 2 != 0)) {
			std::cerr
				<< "Warning: the linear numbers of sites in x should be even."
				<< std::endl;
			return 2;
		}

		withCouplings = (J != 0.);
	} catch (std::exception &e) {
		std::cerr << "Error: " << e.what() << std::endl;
		return 2;
	}

	return 0;
}
