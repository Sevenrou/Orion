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


long StringTokenCounter( const std::string& Source)
{
	long NombreValeurs = 0;
	size_t pos1 = 0, pos2 = 0;

	while( ( pos2 = Source.find_first_of( " ,", pos2 + 1)) != std::string::npos)
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
						const long NombreDimensions)
{
	long i;

	for( i = 0; i < NombreDimensions - 1; i++)
	{
		double pouet;
		Source >> pouet;
		matrice[i*NombrePoints + NumeroPoint] = pouet;
		Source.get();
	}

	Source >> matrice[i*NombrePoints + NumeroPoint];
}

bool ParseFileIntoMatrix( const char* FileName, size_t NumDimExpected, double*& matrice, long& NombrePoints, long& NombreDimensions)
{
	std::ifstream FichierEntree( FileName);

	if( ! FichierEntree.is_open())
		return false;

	std::string ligne;

	getline( FichierEntree, ligne);
	NombreDimensions = StringTokenCounter( ligne);

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
		NombreTemp = StringTokenCounter( ligne);

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

		if( NombreTemp != NombreDimensions)
		{
			std::cout << "Line " << NombrePoints + 1 << " : " << NombreTemp << " dimensions instead of " << NombreDimensions << std::endl;
			return false;
		}

		NombrePoints++;
	}

	std::cout << "Loading " << NombrePoints << " elements in " << NombreDimensions << " dimensions" << std::endl;

	matrice = new double[NombrePoints*NombreDimensions];

	FichierEntree.close();
	FichierEntree.open( FileName);

	for( long i = 0; i < NombrePoints; i++)
	{
		getline( FichierEntree, ligne);
		std::istringstream iss( ligne);
		StringTokenizer( iss, matrice, i, NombrePoints, NombreDimensions);
	}

	return true;
}