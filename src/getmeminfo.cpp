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


#include <cstdio>
#include <sstream>

#include "getmeminfo.h"


unsigned long GetMemInfo()
{
	std::ostringstream oss2;
	oss2 << "ps h -p " << getpid() << " -o vsz";

	FILE* psret = popen( oss2.str().c_str(), "r");
	if( psret == NULL)
		return 0;

	unsigned long vmval;
	if( fscanf( psret, " %lu", &vmval) != 1)
		return 0;

	return vmval;
}
