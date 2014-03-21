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


#include "readfile.h"


long StringTokenCounter( const std::string& Source, char delim)
{
	long NombreValeurs = 0;
	size_t pos1 = 0, pos2 = 0;

	while( ( pos2 = Source.find_first_of( delim, pos2 + 1)) != std::string::npos)
	{
		if( pos1 == pos2)
			return -pos1;
		pos1 = pos2 + 1;
		NombreValeurs++;
	}

	if( ! (pos1 == 0 && pos2 == std::string::npos))
		NombreValeurs++;

	return NombreValeurs;
}

void StringTokenizer(	std::istringstream& Source,
						double* matrice,
						const long NumeroPoint,
						const long NombrePoints,
						const long NombreDimensions,
						char delim,
						std::vector<std::string>* Labels)
{
	long i;

	if( Labels)
	{
		std::string label;
		size_t possep = Source.str().find_first_of( delim);
		Source.width(possep);
		Source >> label;
		Labels->push_back(label);
		Source.get();
	}

	for( i = 0; i < NombreDimensions; i++)
	{
		double pouet;
		std::string pouetTemp;
		Source >> pouetTemp;
		if(!pouetTemp.compare("inf")){
            		pouet = std::numeric_limits<double>::max();
		}else{
            		std::istringstream iss(pouetTemp);
            		iss >> pouet;
		}
		matrice[i*NombrePoints + NumeroPoint] = pouet;
		Source.get();
	}
}

bool ParseFileIntoMatrix(	const char* FileName,
							size_t NumDimExpected,
							double*& matrice,
							long& NombrePoints,
							long& NombreDimensions,
							std::vector<std::string>* Labels )
{
	std::ifstream FichierEntree( FileName);

	if( ! FichierEntree.is_open())
		return false;

	std::string ligne;
	getline( FichierEntree, ligne);

	int posdelim = ligne.find_first_of( " ,");
	if(posdelim == std::string::npos)
	{
		std::cout << "No delimiter (either space or comma) found in the first line" << std::endl;
		return false;
	}

	char delim = ligne[posdelim];
	std::cout << "Delimiter found: ";
	switch(delim)
	{
		case ' ':	std::cout << "space";
					break;
		case ',':	std::cout << "comma";
					break;
		default :	std::cout << "unknown!?";
	}
	std::cout << std::endl;

	NombreDimensions = StringTokenCounter( ligne, delim);

	// Is the first column for labels?
	if( Labels)
		NombreDimensions--;

	if( NombreDimensions == 0)
		return false;

	if( NumDimExpected > 0 && NumDimExpected != static_cast<size_t>(NombreDimensions))
	{
		std::cout	<< "Different number of dimensions between the -s argument (" << NumDimExpected
					<< ") and the dataset (" << NombreDimensions << ")" << std::endl;
		return false;
	}

	long NombreTemp;
	NombrePoints = 1;

	while( ! FichierEntree.eof())
	{
		getline( FichierEntree, ligne);
		NombreTemp = StringTokenCounter( ligne, delim);

		if( NombreTemp < 0)
		{
			std::cout << "Empty field line " << NombrePoints + 1 << " character " << -NombreTemp << std::endl;
			return false;
		}

		if( NombreTemp == 0)
		{
			std::cout << "Line " << NombrePoints + 1 << " empty" << std::endl;
			break;
		}

		if( Labels)
			NombreTemp--;

		if( NombreTemp != NombreDimensions)
		{
			std::cout << "Line " << NombrePoints + 1 << " : " << NombreTemp << " dimensions instead of " << NombreDimensions << std::endl;
			return false;
		}

		NombrePoints++;
	}

	std::cout << "Loading " << NombrePoints << " elements in " << NombreDimensions << " dimensions" << std::endl;

	matrice = new double[NombrePoints*NombreDimensions];

	if( Labels)
		Labels->reserve(NombrePoints);

	FichierEntree.close();
	FichierEntree.open( FileName);

	for( long i = 0; i < NombrePoints; i++)
	{
		getline( FichierEntree, ligne);
		std::istringstream iss( ligne);
		StringTokenizer( iss, matrice, i, NombrePoints, NombreDimensions, delim, Labels);
	}

	return true;
}
