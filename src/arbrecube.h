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


#ifndef ARBRECUBE_H_
#define ARBRECUBE_H_


#include <map>
#include <boost/unordered_map.hpp>
#include <boost/pool/poolfwd.hpp>
#include <boost/pool/pool_alloc.hpp>

#include "defs.h"
#include "stx/btree_multimap.h"
#include "utils.h"


struct Noeud {
	Noeud( unsigned long NbEnfants, unsigned long PathSize) :
		EstType1( false),
		EstComplet( false),
		cs_D( 0),
		cs_I( 0)
	{
		enfants.reserve( NbEnfants);
		Chemin.reserve( PathSize);
	}

	~Noeud()
	{
		std::for_each( enfants.begin(), enfants.end(), DeleteObject());
	}

	void Fill_D_I( const DotSet& TempD, const CombinedSkyline& TempI);

	Noeud* parent;
	std::vector<Noeud*> enfants;
	long offset;
	bool EstType1;
	bool EstComplet;		// all points are in D

	std::vector<long> Chemin;

	CompactSet	cs_D;
	CompactComb	cs_I;

	DotSet depthD;
	CombinedSkyline depthI;
};

struct Closure {
	Closure(const LatticePath& path, long NombreDimensions);

	void AddElement(const LatticePath& path);

	std::vector<LatticePath> ClosedNodes;
	std::vector<LatticePath> Generators;

	std::vector<bool> PrunedDimensions;
};

struct ClosureHash : public std::unary_function<HashKey,size_t>
{
	std::size_t operator()(const HashKey& n) const
	{
		std::size_t seed = 0;
		boost::hash_combine(seed, boost::hash_range(n.first.begin(), n.first.end()));
		for(HashKey::second_type::const_iterator iteI = n.second.begin(); iteI != n.second.end(); ++iteI)
			boost::hash_combine(seed, boost::hash_range(iteI->begin(), iteI->end()));

		return seed;
	}
};

struct ClosureEq : public std::binary_function<HashKey,HashKey,bool>
{
	bool operator()(const HashKey& n1, const HashKey& n2) const
	{
		return n1.first == n2.first && n1.second == n2.second;
	}
};

typedef boost::unordered_map<HashKey,Closure*,ClosureHash,ClosureEq> HashClosure;


struct UnParent
{
	Noeud* Parent;
	long RemovedDim;
};

typedef std::vector<UnParent,boost::fast_pool_allocator<UnParent> > ParentsList;

struct UnDotSet
{
	const CompactSet* CombDotSet;
	long RemovedDim;
};

typedef std::vector<UnDotSet,boost::fast_pool_allocator<UnDotSet> > UnDotSetList;


enum PointOrderRelation { P1_DOM_P2, P2_DOM_P1, EQUIV, UNCOMP };


class ArbreCube
{
public:
	ArbreCube( const double* matrice_p, const std::vector<bool>& FindLowest, long NombrePoints_p, long NombreDimensions_p);
	~ArbreCube();

	bool LoadClosures( std::istream& Cin);

	void DepthAlgo( bool ComputeLast);
	void BreadthAlgo( bool UseClosure);

	void AfficheResultat( std::ostream& Cout_Resultat, std::vector<std::string>* Labels) const;
	void AfficheClos( std::ostream& Cout, std::vector<std::string>* Labels) const;

	// makes sense only for depth algo
	unsigned long GetNbProcessedNodes()	const		{ return Compteur; }
	size_t GetNbClos() const;

	unsigned long GetNbType1() const				{ return Type1Count; }
	uint64_t GetNbSkylineFoundDirectly() const		{ return NbSkylineFoundDirectly; }
	uint64_t GetNbSkylineFoundByBNL() const			{ return NbSkylineFoundTotal; }

	static std::ostream& PrintPath( const std::vector<long>& Chemin, std::ostream& Cout);

	template<typename T>
	static void PrintD(const T& TempD, std::ostream& Cout);
	template<typename TI, typename TD>
	static void PrintI(const TI& TempI, std::ostream& Cout);

	void HashStat( std::ostream& Cout) const;

private:
	void GenereBTrees();
	void GenereDimension1(bool depth);
	void GenereDimensionProf( const Noeud& ParentNoeud, const DotSet& pcsD, const CombinedSkyline& pcsI, Closure* ParentGroup = 0);
	void GenereDimensionN( unsigned long DimNumber);

	/*
	 * Determines whether TempNoeud is of type 1 and fills it accordingly
	 */
	template<class _D, class _I, class D2, class I2>
	void ManageType1( Noeud* TempNoeud,
								DotSet& TempD,
								CombinedSkyline& TempI,
								const _D& pcsD,
								const _I& pcsI,
								const D2& newDimD,
								const I2& newDimI);

	// depth-related methods
		/*
		 * This step processes the parent node's sets of combined points to remove those
		 * - already in the distinct list of the current node
		 * - that can be distinguished by the newly added dimension
		 */
		void Depth_Step_2_1( DotSet& TempD,
							CombinedSkyline& TempI,
							const Noeud* ParentNoeud,
							const Noeud* TempNewDim);

		/*
		 * This step processes the new dimension's node's sets of combined points to remove those
		 * - already in the distinct list of the current node
		 * - that can be distinguished by the newly added dimension
		 */
		void Depth_Step_2_2( DotSet& TempD,
							CombinedSkyline& TempI,
							const Noeud* TempNewDim,
							const std::vector<long>& Path);

		/*
		 * This method is used in Step 3.1 to unite elements which belong
		 * to the computed intervals
		 */
		void Union( DotSet& sk,
					stx::btree_multimap<double,long>::const_iterator& iteLowerBound,
					stx::btree_multimap<double,long>::const_iterator& iteUpperBound);

		/*
		 * Compute the node of all dimensions
		 * Useful for closures since this node will by definition alway be a closure
		 * Can accelerate a bit the creation of the left-most branch of the tree
		 */
		const Noeud* ComputeLastNode();

	// breadth-related methods
		/*
		 * This method creates a list of all the direct parents of the current node
		 * Each element of a list comprises the parent node and the missing dimension compared to the actual node
		 */
		void CreateParentsList( Noeud* TempNoeud, ParentsList& ListeComposantes);

		/*
		 * This step simplifies the list of combined sets by removing those
		 * having an element already in the current node's distinct list
		 */
		void Breadth_Step_2_1(	UnDotSetList& VecUDS,
								const ParentsList& ListeComposantes,
								const Noeud* TempNoeud,
								const DotSet& TempD);

		/*
		 * This step checks for sets of combined points to see whether the newly
		 * added dimension can help distinguish them. This means that we're trying
		 * to remove some points that become dominant with a new dimension added
		 */
		void Breadth_Step_2_2(	const UnDotSetList& VecUDS,
								DotSet& TempD,
								CombinedSkyline& TempI);

		/*
		 * This step takes each element of every dimension within the range defined by the known skyline points
		 * Every dimension set is intersected with each other
		 * Remaining elements are added in node's skyline
		 */
		void RangeBNL( const Noeud* TempNoeud, DotSet& TempD, CombinedSkyline& TempI);

		/*
		 * Visit every dimension and for each of them define an interval from which to pick up elements
		 * If we are using a breadth-first algorithm, we do intersections
		 * If we are using a depth-first algorithm, we do unions
		 */
		void Step_3_1(	DotSet& TempDotSet,
						const DotSet& TempD,
						const CombinedSkyline& TempI,
						const std::vector<long>& Path);

		/*
		 * Copies the contents of Stockage into D and I of TempNoeud
		 */
		void Step_3_3( const std::list<DotSet>& Stockage, DotSet& TempD, CombinedSkyline& TempI);

		/*
		 * This method is used in Step 3.1 to compute the intersection
		 * of elements which belong to the computed intervals
		 */
		void Intersection(	DotSet& sk,
							stx::btree_multimap<double,long>::const_iterator& iteLowerBound,
							stx::btree_multimap<double,long>::const_iterator& iteUpperBound);


	/*
	 * Takes the list of points retrieved from step 3.1 and remove those already existing in the
	 * distinct and combined set lists of elements of the current node
	 */
	void RemoveIncluded( DotSet& InterDotSet,
						const DotSet& TempD,
						const CombinedSkyline& TempI);

	/*
	 * Returns true if the intersection is empty, false otherwise
	 */
	template<class In, class In2>
	bool AreDisjoint( In first, In last, In2 first2, In2 last2);

	/*
	 * Compares P1 and P2 on the dimensions provided by Chemin. Returns:
	 * - P1_DOM_P2 if P1 dominates P2
	 * - P2_DOM_P1 if P2 dominates P1
	 * - EQUIV if they are combined
	 * - UNCOMP if they aren't comparable
	 */
	PointOrderRelation ComparePoints( long P1, long P2, const std::vector<long>& Chemin);

	/*
	 * Compares all elements from TempDotSet and ResultDotSet between them
	 * Chemin provides the list of dimensions on which to compare the elements
	 * Stores the result in ResultDotSet
	 */
	template<class In>
	void BNL( In first, In last, const std::vector<long>& Chemin, std::list<DotSet>& ResultDotSet);

	// domain-related methods
		/*
		 * Main method that loads SP and calls Evaluate for every element of a selected B+-tree
		 */
		void TakeTheBus( DotSet& TempD, CombinedSkyline& TempI, const Noeud* TempNoeud);

		/*
		 * Evaluate compares the element to those of SP by the way of the sum its values on the selected dimensions
		 */
		long Evaluate(	long NumPoint,
						std::multimap<double,long> &SP,
						const std::vector<long>& Chemin,
						DotSet& TempD,
						CombinedSkyline& TempI);

		/*
		 * Computes the sum of an element's values on the selected dimensions
		 */
		double Sum( long NumPoint, const std::vector<long>& Chemin);

		/*
		 * Returns the domain size of a dimension (max - min)
		 */
		double getDomaineSize(int dimension);

	// puts all elements from D and I into Omega
	template<class In, class In2>
	void Flatten( In first, In last, In2 first2, In2 last2, DotSet& Result);

	// closure-related methods
		/*
		 * Returns true if Skyline(TempNoeud) == Skyline(N), for any N element of VecNoeuds
		 * Returns true if all elements of VecNoeuds are null
		 * Returns false otherwise
		 *
		 * Distinct and combined properties matter here
		 */
		template<class T>
		bool IncludedIn( const Noeud* TempNoeud, const T& VecNoeuds);

		/*
		 *
		 */
		template<class InputIterator>
		bool PathIncludedIn(InputIterator begin, InputIterator end, const std::vector<long>& myPath);

		/*
		 * Returns the node for which Skyline(TempNoeud) == Skyline(N), for any N element of VecNoeuds
		 * Returns null otherwise
		 *
		 * Distinct and combined properties matter here
		 */
		Closure* FindClosure(const HashKey& keyToFind, const HashClosure& VecNoeuds);

	// display methods
	void AfficheLargeur(std::vector<Noeud*>& Pile,
						std::ostream& Cout,
						std::vector<std::string>* Labels) const;
	void AfficheNoeud(	const Noeud *TempNoeud,
						std::vector<Noeud*>& Pile,
						std::ostream& Cout,
						std::vector<std::string>* Labels) const;

	template<class _D, class _I>
	void AfficheSkyline(const _D& cs_D,
						const _I& cs_I,
						bool isComplete,
						std::ostream& Cout,
						std::vector<std::string>* Labels) const;

	ALGO CurrentAlgo;
	const double* const matrice;
	const long NombrePoints;
	const long NombreDimensions;
	const std::vector<bool>& FindLowest;

	std::vector<int> FindLowest1;
	std::vector<int> FindLowest2;

	Cnk MyCnk;
	std::vector<double> domainSize;

	std::vector<stx::btree_multimap<double,long> > VecBtree;

	HashClosure MesNoeudClos;

	Noeud* racine;
	const Noeud* lastNode;
	std::vector<Noeud*> PileNm1;
	unsigned long Compteur;					// processed nodes
	unsigned long Type1Count;
	uint64_t NbSkylineFoundDirectly;
	uint64_t NbSkylineFoundTotal;
	bool isLastNodeFirst;	 				// whether to compute the last node first (only relevant for depth)
};


template<class In, class In2>
bool ArbreCube::AreDisjoint( In first, In last, In2 first2, In2 last2) {
	if( first == last || first2 == last2)
		return true;

	while( first != last && first2 != last2)
		if( *first < *first2)
			++first;
		else
		{
			if( *first > *first2)
				++first2;
			else
				return false;
		}

	return true;
}


template<class In>
void ArbreCube::BNL( In first, In last, const std::vector<long>& Chemin, std::list<DotSet>& ResultDotSet)
{
	std::vector<long>::const_iterator iteChemin;
	std::list<DotSet>::iterator iteResult;

	DotSet EmptyDotSet;
	bool MustInsert;

	for( ; first != last; ++first)
	{
		MustInsert = true;

		for( iteResult = ResultDotSet.begin(); iteResult != ResultDotSet.end(); )
		{
			PointOrderRelation ResultComp = ComparePoints( *first, *((*iteResult).begin()), Chemin);

			if( ResultComp == P2_DOM_P1)
			{
				MustInsert = false;
				break;
			}
			if( ResultComp == EQUIV)
			{
				(*iteResult).insert( *first);
				MustInsert = false;
				break;
			}
			if( ResultComp == P1_DOM_P2)
				ResultDotSet.erase( iteResult++);
			else
				++iteResult;
		}

		if( MustInsert)
			ResultDotSet.insert( ResultDotSet.begin(), EmptyDotSet)->insert( *first);
	}
}


template<class In, class In2>
void ArbreCube::Flatten( In first, In last, In2 first2, In2 last2, DotSet& Result)
{
	for( ; first != last; ++first)
		Result.insert( *first);

	for( ; first2 != last2; ++first2)
		Result.insert( first2->begin(), first2->end());
}

template<class InputIterator>
bool ArbreCube::PathIncludedIn(InputIterator begin, InputIterator end, const std::vector<long>& myPath) {
	for( ; begin != end; ++begin) {
		if(std::includes((*begin)->Chemin.begin(), (*begin)->Chemin.end(), myPath.begin(), myPath.end()))
			return true;
	}
	return false;
}

template<typename T>
void ArbreCube::PrintD(const T& TempD, std::ostream& Cout) {
	if( TempD.empty())
		return;

	typename T::const_iterator iteDs = TempD.begin();
	for(;;)
	{
		Cout << 'e' << (*iteDs) + 1;
		if( ++iteDs != TempD.end())
			Cout << ',';
		else
			return;
	}
}


template<typename TI,typename TD>
void ArbreCube::PrintI(const TI& TempI, std::ostream& Cout) {
	if( TempI.empty())
		return;

	typename TI::const_iterator iteSk = TempI.begin();
	for(;;)
	{
		typename TD::const_iterator iteDs = (*iteSk).begin();
		for(;;)
		{
			Cout << 'e' << *iteDs + 1;
			if( ++iteDs != (*iteSk).end())
				Cout << "-";
			else
				break;
		}

		if( ++iteSk != TempI.end())
			Cout << ',';
		else
			break;
	}
}

template<class _D, class _I, class D2, class I2>
void ArbreCube::ManageType1( Noeud* TempNoeud,
							DotSet& TempD,
							CombinedSkyline& TempI,
							const _D& pcsD,
							const _I& pcsI,
							const D2& newDimD,
							const I2& newDimI)
{
	if( ! pcsD.empty())
	{
		long D_Value = *(pcsD.begin());

		if( ! newDimD.empty())
		{
			if( D_Value == *(newDimD.begin()))
			{
				TempD.insert( D_Value);
				TempNoeud->EstType1 = true;
			}
		}
		else
		{
			if( std::binary_search( (*(newDimI.begin())).begin(),
									(*(newDimI.begin())).end(),
									D_Value) )
			{
				TempD.insert( D_Value);
				TempNoeud->EstType1 = true;
			}
		}
	}
	else
	{
		if( ! newDimD.empty())
		{
			long D_Value = *(newDimD.begin());

			if( std::binary_search( (*(pcsI.begin())).begin(),
									(*(pcsI.begin())).end(),
									D_Value) )
			{
				TempD.insert( D_Value);
				TempNoeud->EstType1 = true;
			}
		}
		else
		{
			DotSet TempDotSet;
			set_intersection(	(pcsI.begin())->begin(),
					(pcsI.begin())->end(),
					(newDimI.begin())->begin(),
					(newDimI.begin())->end(),
					inserter( TempDotSet, TempDotSet.begin()) );

			// type 1
			if( ! TempDotSet.empty())
			{
				TempNoeud->EstType1 = true;
				if( TempDotSet.size() > 1)
					TempI.insert( TempDotSet);
				else
					TempD = TempDotSet;
			}
		}
	}

	if( TempNoeud->EstType1)
		Type1Count++;
}


template<class _D, class _I>
void ArbreCube::AfficheSkyline( const _D& cs_D,
								const _I& cs_I,
								bool isComplete,
								std::ostream& Cout,
								std::vector<std::string>* Labels ) const
{
	if(isComplete)
	{
		if( Labels)
		{
			Cout << (*Labels)[0];
			for( long i = 1; i < NombrePoints; i++)
				Cout << ',' << (*Labels)[i];
		}
		else
		{
			Cout << "e0";
			for( long i = 1; i < NombrePoints; i++)
				Cout << ',' << 'e' << i;
		}
		return;
	}

	typename _D::const_iterator iteDs;
	if( ! cs_D.empty())
	{
		for( iteDs = cs_D.begin();;)
		{
			if( Labels)
				Cout << (*Labels)[*iteDs];
			else
				Cout << 'e' << (*iteDs);
			if( ++iteDs != cs_D.end())
				Cout << ',';
			else
				break;
		}
	}

	typename _I::const_iterator iteSk;
	if( ! cs_I.empty())
	{
		if( ! cs_D.empty())
			Cout << ',';

		for( iteSk = cs_I.begin();;)
		{
			for( iteDs = (*iteSk).begin();;)
			{
				if( Labels)
					Cout << (*Labels)[*iteDs];
				else
					Cout << 'e' << (*iteDs);
				if( ++iteDs != (*iteSk).end())
					Cout << "-";
				else
					break;
			}

			if( ++iteSk != cs_I.end())
				Cout << ',';
			else
				break;
		}
	}
}


#endif /*ARBRECUBE_H_*/
