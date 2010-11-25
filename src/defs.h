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


#ifndef DEFS_H_
#define DEFS_H_


#include <vector>
#include <list>
#include <set>


typedef std::set<long> DotSet;

struct ltdotset
{
	bool operator()(const DotSet& s1, const DotSet& s2) const
	{
		if( s1.size() < s2.size())
			return true;

		if( s1.size() > s2.size())
			return false;

		DotSet::const_iterator iteS1 = s1.begin();
		DotSet::const_iterator iteS2 = s2.begin();
		while( iteS1 != s1.end())
		{
			if( *iteS1 < *iteS2)
				return true;
			if( *iteS1 > *iteS2)
				return false;

			++iteS1;
			++iteS2;
		}
		return false;
	}
};

typedef std::set<DotSet, ltdotset> CombinedSkyline;

typedef std::list<DotSet> TempCombined;


// once the node is computed, store the result in more compact containers
typedef std::vector<long> CompactSet;
typedef std::vector<CompactSet> CompactComb;


// for closures (and not only?)
typedef std::vector<long> LatticePath;
typedef std::pair<DotSet,CombinedSkyline> HashKey;


// type of algorithm used to generate the tree
enum ALGO { DEPTH, BREADTH, BR_DOM };


struct DeleteObject
{
	template <typename T>
	void operator()(T *ptr) { delete ptr; }
};


#endif // DEFS_H_
