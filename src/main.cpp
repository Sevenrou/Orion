/*
 * Orion Skycube Computing v1.0
 * Copyright (C) 2010 Thomas Kister
 *
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation; either version 2.1 of the License, or (at your
 * option) any later version.
 *
 * This library is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License
 * for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 */


#include <cstring>
#include <vector>

#include "readfile.h"
#include "arbrecube.h"
#include "getmeminfo.h"

#ifdef __APPLE__
	#include "os/mac_clock_gettime.h"
#endif // __APPLE__


#define SYNTAX	"UnifiedSkyCube [-s UpperOrLower] [-l] -a Algo Filename\n\n"\
				"  -s UPPERORLOWER\tString consisting of either u or l for each dimension\n"\
				"\t\t\tBy default it is l (lowest) for each dimension\n\n"\
				"  -label\t\tIndicates that the first column contains a label for\n"\
				"\t\t\teach element. These labels are used in the output to\n"\
				"\t\t\treplace the line number\n\n"\
				"  -a ALGO={depth,breadth,br_dom}\n"\
				"\t\t\tChooses the algorithm to use (by default: depth):\n"\
				"\t\t\t- depth : creates the tree recursively with closures\n"\
				"\t\t\t- breadth : creates the tree using a stack\n"\
				"\t\t\t- br_dom : breadth + domain optimizations\n\n"\
				"  -nolast\t\tDoes not compute the node of all dimensions first\n"\
				"\t\t\tApplies only to depth algorithm, no effect otherwise\n\n"\
				"  -h,--help\t\tDisplays this help message\n\n\n"\
				"Examples: UnifiedSkyCube -a depth -nolast -s uulluuluulluull dataset42.txt\n"\
				"          UnifiedSkyCube -a breadth dataset51.txt\n\n"


void OutputResult(	const ArbreCube& bouleau,
					const std::vector<bool>& FindLowest,
					ALGO SelectedAlgo,
					std::string& FichierResultat)
{
	for( std::vector<bool>::const_iterator iteUL = FindLowest.begin(); iteUL != FindLowest.end(); ++iteUL)
		FichierResultat += (*iteUL) ? "l" : "u";
	FichierResultat += ".res";

	std::ofstream OutputStream( FichierResultat.c_str());
	if( ! OutputStream.is_open())
	{
		std::cout << "Couldn't open file " << FichierResultat << " to store the result, sending to screen" << std::endl;

		if( SelectedAlgo == DEPTH)
			bouleau.AfficheClos( std::cout);
		else
			bouleau.AfficheResultat( std::cout);
	}
	else
	{
		std::cout << "Saving result into file " << FichierResultat << std::endl;

		if( SelectedAlgo == DEPTH)
			bouleau.AfficheClos( OutputStream);
		else
			bouleau.AfficheResultat( OutputStream);
	}
}


int main(int argc, char **argv)
{
	if( argc < 2)
	{
		std::cout << SYNTAX << std::endl;
		return 1;
	}

	std::vector<bool> FindLowest;
	ALGO SelectedAlgo = DEPTH;
	bool ComputeLastNode = true;
	bool HasLegend = false;

	long NumArg = 0;
	while( ++NumArg < argc - 1)
	{
		if( strcmp( argv[NumArg], "-s") == 0)
		{
			if( ++NumArg == argc - 1)
			{
				std::cout << SYNTAX << std::endl;
				return 1;
			}

			for( unsigned long j = 0; j < strlen( argv[NumArg]); j++)
				switch( argv[NumArg][j])
				{
					case 'l':	FindLowest.push_back(true);
								break;
					case 'u':	FindLowest.push_back(false);
								break;
					default :	std::cout << SYNTAX << std::endl;
								return 1;
				}
		}
		else if( strcmp( argv[NumArg], "-a") == 0)
		{
			if( ++NumArg == argc - 1)
			{
				std::cout << SYNTAX << std::endl;
				return 1;
			}

			if( strcmp( argv[NumArg], "depth") == 0)
				SelectedAlgo = DEPTH;
			else if( strcmp( argv[NumArg], "breadth") == 0)
				SelectedAlgo = BREADTH;
			else if( strcmp( argv[NumArg], "br_dom") == 0)
				SelectedAlgo = BR_DOM;
			else
			{
				std::cout << SYNTAX << std::endl;
				return 1;
			}
		}
		else if( strcmp( argv[NumArg], "-nolast") == 0)
		{
			ComputeLastNode = false;
		}
		else if( strcmp( argv[NumArg], "-label") == 0)
		{
			HasLegend = true;
		}
		else
		{
			std::cout << SYNTAX << std::endl;
			return 1;
		}
	}

	if( strcmp( argv[NumArg], "-h") == 0 || strcmp( argv[NumArg], "--help") == 0)
	{
		std::cout << SYNTAX << std::endl;
		return 1;
	}

	double* matrice = 0;
	long NombrePoints = 0;
	long NombreDimensions = 0;

	if( ! ParseFileIntoMatrix( argv[NumArg], FindLowest.size(), matrice, NombrePoints, NombreDimensions))
	{
		std::cout << "Error in data set parsing" << std::endl;
		return 1;
	}

	if( FindLowest.size() == 0)
		for( int j = 0; j < NombreDimensions; j++)
			FindLowest.push_back(true);

	ArbreCube bouleau( matrice, FindLowest, NombrePoints, NombreDimensions);

	// we start creating the result filename here to use the switch once
	std::string FichierResultat( argv[NumArg]);

	std::cout << "Selected algorithm ";
	switch( SelectedAlgo)
	{
		case DEPTH: 	if( ComputeLastNode)
						{
							std::cout << "Orion-Clos: Depth with last node first" << std::endl;
							FichierResultat += ".del.";
						}
						else
						{
							std::cout << "Depth" << std::endl;
							FichierResultat += ".dep.";
						}
						bouleau.DepthAlgo(ComputeLastNode);
						std::cout << "Processed nodes: " << bouleau.GetNbProcessedNodes() << std::endl;
						std::cout << "Closure nodes: " << bouleau.GetNbClos() << "/";
						std::cout <<  (1 << NombreDimensions) - 1 << std::endl;
						std::cout << "Skyline points found directly / by BNL: " << bouleau.GetNbSkylineFoundDirectly();
						std::cout << "/" << bouleau.GetNbSkylineFoundByBNL() << std::endl;
						break;
		case BREADTH:	std::cout << "Orion: Breadth" << std::endl;
						bouleau.BreadthAlgo(false);
						FichierResultat += ".bre.";
						break;
		case BR_DOM:	std::cout << "Orion-Tail: Breadth with domains" << std::endl;
						bouleau.BreadthAlgo(true);
						FichierResultat += ".brd.";
						break;
	}
	std::cout << "Type I nodes: " << bouleau.GetNbType1() << "/" <<  (1 << NombreDimensions) - 1 << std::endl;

	/////////////////////////////////////////////////////////////////////////////////////////
	// Stat of the day

	std::cout << "Peak memory usage: " << GetMemInfo() << std::endl;

#ifndef BENCH
	/////////////////////////////////////////////////////////////////////////////////////////
	// Storing result somewhere

	OutputResult( bouleau, FindLowest, SelectedAlgo, FichierResultat);

#endif // BENCH

	delete[] matrice;

	return 0;
}
