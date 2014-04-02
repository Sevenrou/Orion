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


#ifndef READFILE_H_
#define READFILE_H_


#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <limits>


bool ParseFileIntoMatrix(	const char* FileName, size_t NumDimExpected,
							double*& matrice,
							long& NombrePoints,
							long& NombreDimensions,
							std::vector<std::string>* Labels = 0 );


#endif /*READFILE_H_*/
