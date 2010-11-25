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


#include "utils.h"


void Utils::GetTime( timespec& now)
{
	clock_gettime( CLOCK_MONOTONIC, &now);
}

timespec Utils::GetTime()
{
	timespec now;
	clock_gettime( CLOCK_MONOTONIC, &now);
	return now;
}

timespec Utils::GetDiffTime( const timespec& beg, const timespec& end)
{
	timespec retval;

	retval.tv_sec = end.tv_sec - beg.tv_sec;
	retval.tv_nsec = end.tv_nsec - beg.tv_nsec;
	if( retval.tv_nsec < 0)
	{
		retval.tv_nsec += 1000000000;
		retval.tv_sec--;
	}

	return retval;
}

void Utils::DisplayTime( const timespec& disp, std::ostream& Cout)
{
	Cout << disp.tv_sec << ".";
	Cout.width(9);
	Cout.fill('0');
	Cout << disp.tv_nsec << " second(s)";
}


Cnk::Cnk( unsigned long p_NbDim) :
	NbDim( p_NbDim),
	VecCnk( p_NbDim+1, 0)
{
	VecCnk[0] = 1;
	for( unsigned long i = 0; i < NbDim; i++)
		for( unsigned long j = NbDim; j > 0; j--)
			VecCnk[j] += VecCnk[j-1];
}

unsigned long Cnk::Value( unsigned long index)
{
	return VecCnk.at(index);
}
