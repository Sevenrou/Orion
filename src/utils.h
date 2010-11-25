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


#ifndef UTILS_H_
#define UTILS_H_


#include <iostream>
#include <vector>
#include <ctime>

#ifdef __APPLE__
	#include "os/mac_clock_gettime.h"
#endif // __APPLE__


class Utils
{
public:
	static void GetTime( timespec& now);

	static timespec GetTime();

	static timespec GetDiffTime( const timespec& beg, const timespec& end);

	static void DisplayTime( const timespec& disp, std::ostream& Cout);
};


class Cnk
{
public:
	Cnk( unsigned long p_NbDim);

	unsigned long Value( unsigned long index);

private:
	const unsigned long NbDim;
	std::vector<unsigned long> VecCnk;
};


#endif // UTILS_H_
