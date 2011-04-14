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


#include "arbrecube.h"


std::ostream& operator<<(std::ostream &Cout, const std::vector<long>& Path) {
	for( std::vector<long>::const_iterator reve = Path.begin(); reve != Path.end(); ++reve)
		Cout << static_cast<unsigned char>(*reve + 65);
	return Cout;
}


Closure::Closure(const LatticePath& path, long NombreDimensions) :
	PrunedDimensions( NombreDimensions, false)
{
	ClosedNodes.push_back(path);
	Generators.push_back(path);

	// Also update the dimensions to prune
	for(std::vector<long>::const_iterator it = path.begin(); it != path.end(); ++it)
		PrunedDimensions[*it] = true;
#ifdef DEBUG_CLOS
	std::cout << "New closure with " << path << std::endl;
#endif
}

void Closure::AddElement(const LatticePath& path)
{
#ifdef DEBUG_CLOS
	std::cout << "\t\tAddElement: " << path << std::endl;
#endif

	std::vector<LatticePath>::const_iterator iteClos = ClosedNodes.begin();
	for( ; iteClos != ClosedNodes.end(); ++iteClos)
	{
#ifdef DEBUG_CLOS
		std::cout << "\t\t(*iteClos): " << *iteClos << std::endl;
#endif
		if( std::includes(	iteClos->begin(), iteClos->end(),
							path.begin(), path.end()))
			break;
	}

	// The node belongs to the closure but isn't "closed" so it's a new closed node
	if( iteClos == ClosedNodes.end())
	{
		ClosedNodes.push_back(path);

		// Also update the dimensions to prune
		for(std::vector<long>::const_iterator it = path.begin(); it != path.end(); ++it)
			PrunedDimensions[*it] = true;
#ifdef DEBUG_CLOS
		std::cout << "Added as close" << std::endl;
#endif
		return;
	}

	// The node is "closed" so it could be a generator
	std::vector<LatticePath>::iterator iteNoeuds = Generators.begin();
	while( iteNoeuds != Generators.end())
	{
#ifdef DEBUG_CLOS
		std::cout << "\t\tGen: " << *iteNoeuds << std::endl;
#endif
		// ok it isn't a generator, we just ignore it
		if( iteNoeuds->size() < path.size() &&
			std::includes( path.begin(), path.end(),
							iteNoeuds->begin(), iteNoeuds->end() ) )
		{
#ifdef DEBUG_CLOS
			std::cout << "\t\tIsn't gen" << std::endl;
#endif
			return;
		}

		if( iteNoeuds->size() > path.size() &&
			std::includes( iteNoeuds->begin(), iteNoeuds->end(),
							path.begin(), path.end() ) )
		{
#ifdef DEBUG_CLOS
			std::cout << "\t\tErase " << *iteNoeuds << std::endl;
#endif
			iteNoeuds = Generators.erase(iteNoeuds);
		}
		else
			++iteNoeuds;
	}
#ifdef DEBUG_CLOS
	std::cout << "\t\tAdded as generator" << std::endl;
#endif
	Generators.push_back( path);
}


void Noeud::Fill_D_I( const DotSet& TempD, const CombinedSkyline& TempI)
{
	cs_D.reserve( TempD.size());
	cs_D.assign( TempD.begin(), TempD.end());

	cs_I.reserve( TempI.size());
	CombinedSkyline::const_iterator iteSk;
	for( iteSk = TempI.begin(); iteSk != TempI.end(); ++iteSk)
	{
		cs_I.push_back( CompactSet(0));
		cs_I.back().reserve( (*iteSk).size());
		cs_I.back().assign( (*iteSk).begin(), (*iteSk).end());
	}
}


ArbreCube::ArbreCube( const double* matrice_p, const std::vector<bool>& FindLowest_p, long NombrePoints_p, long NombreDimensions_p) :
	matrice( matrice_p),
	NombrePoints( NombrePoints_p),
	NombreDimensions( NombreDimensions_p),
	FindLowest( FindLowest_p),
	MyCnk( NombreDimensions_p),
	racine( 0),
	lastNode( 0),
	Compteur( 0),
	Type1Count( 0),
	NbSkylineFoundDirectly( 0),
	NbSkylineFoundTotal( 0)
{
	MesNoeudClos.max_load_factor(0.75);

	for(std::vector<bool>::const_iterator it = FindLowest.begin(); it != FindLowest.end(); ++it) {
		FindLowest1.push_back(*it ? 1 : 2);
		FindLowest2.push_back(*it ? 2 : 1);
	}
}

ArbreCube::~ArbreCube()
{
	delete racine;
	for( HashClosure::const_iterator iteClos = MesNoeudClos.begin(); iteClos != MesNoeudClos.end(); ++iteClos)
		delete iteClos->second;
	//std::for_each( MesNoeudClos.begin(), MesNoeudClos.end(), DeleteObject());
	if( lastNode)
		delete lastNode;
}

void ArbreCube::DepthAlgo( bool ComputeLast)
{
	if( racine != 0)			// the tree has already been generated
		throw std::exception();

	CurrentAlgo = DEPTH;
	isLastNodeFirst = ComputeLast;

	std::cout << "Generating B+-Trees... ";
#ifndef DEBUG
	timespec TempusFugit = Utils::GetTime();
#endif // !DEBUG

	GenereBTrees();

#ifndef DEBUG
	timespec TempusFugit2 = Utils::GetTime();
	TempusFugit2 = Utils::GetDiffTime( TempusFugit, TempusFugit2);
	Utils::DisplayTime( TempusFugit2, std::cout);
#endif // !DEBUG
	std::cout << std::endl;

#ifndef DEBUG
	std::cout << "Generating N-Dimension spaces..." << std::endl;
#endif // !DEBUG

	GenereDimension1(true);
	Compteur = NombreDimensions;

	// Since the node of all dimensions is always a closure we compute and add it now
	if( isLastNodeFirst)
	{
		lastNode = ComputeLastNode();
		Compteur++;
	}

	std::vector<Noeud*>::const_iterator iteChildren;
	for( iteChildren = racine->enfants.begin(); iteChildren != racine->enfants.end(); ++iteChildren)
		GenereDimensionProf( **iteChildren, (*iteChildren)->depthD, (*iteChildren)->depthI);

#ifndef DEBUG
	TempusFugit2 = Utils::GetTime();
	TempusFugit2 = Utils::GetDiffTime( TempusFugit, TempusFugit2);
	std::cout << "Computation done in ";
	Utils::DisplayTime( TempusFugit2, std::cout);
	std::cout << " seconds" << std::endl;
#endif // !DEBUG
}

void ArbreCube::BreadthAlgo( bool UseClosure)
{
	if( racine != 0)			// the tree has already been generated
		throw std::exception();

	CurrentAlgo = UseClosure ? BR_DOM : BREADTH;

	std::cout << "Generating B+-Trees... ";
#ifndef DEBUG
	timespec TempusFugit = Utils::GetTime();
#endif // !DEBUG

	GenereBTrees();

#ifndef DEBUG
	timespec TempusFugit2 = Utils::GetTime();
	TempusFugit2 = Utils::GetDiffTime( TempusFugit, TempusFugit2);
	Utils::DisplayTime( TempusFugit2, std::cout);
#endif // !DEBUG
	std::cout << std::endl;

	PileNm1.clear();
	PileNm1.reserve( MyCnk.Value(1));

#ifndef DEBUG
	std::cout << "Generating N-Dimension spaces with N =";
#endif // !DEBUG

	GenereDimension1(false);
	unsigned long NumDim = 2;
	while( ! PileNm1.empty())
		GenereDimensionN( NumDim++);

#ifndef DEBUG
	TempusFugit2 = Utils::GetTime();
	TempusFugit2 = Utils::GetDiffTime( TempusFugit, TempusFugit2);
	std::cout << " in ";
	Utils::DisplayTime( TempusFugit2, std::cout);
	std::cout << " second" << std::endl;
#endif // !DEBUG
}

void ArbreCube::GenereBTrees()
{
	stx::btree_multimap<double,long> TempBTree;

	for( long i = 0; i < NombreDimensions; i++)
	{
		VecBtree.push_back( TempBTree);				// we copy an empty B+-Tree in the vector

		for( long j = 0; j < NombrePoints; j++)
			VecBtree[i].insert( matrice[ i * NombrePoints + j ], j);

		double domValue = getDomaineSize(i);
		domainSize.push_back(domValue);
	}
}

void ArbreCube::GenereDimension1(bool depth)
{
	racine = new Noeud(NombreDimensions,0);
	racine->parent = 0;
	racine->offset = -1;

	Noeud* TempNoeud;

	if( CurrentAlgo == BREADTH || CurrentAlgo == BR_DOM)
#ifndef DEBUG
		std::cout << " 1" << std::flush;
#else
		std::cout << "Generating N-Dimension spaces with N = 1" << std::endl;
#endif // !DEBUG

	for( long i = 0; i < NombreDimensions; i++)
	{
		DotSet TempDotSet;
		CombinedSkyline TempComb;

		TempNoeud = new Noeud(NombreDimensions-i-1,1);
		TempNoeud->parent = racine;
		TempNoeud->offset = i;
		TempNoeud->Chemin.push_back( i);
		TempNoeud->EstType1 = true;
		racine->enfants.push_back( TempNoeud);

		std::pair<stx::btree_multimap<double,long>::const_iterator,stx::btree_multimap<double,long>::const_iterator> Res;

		if( FindLowest[i])
			Res = VecBtree[i].equal_range( VecBtree[i].begin().key());					// the first element holds the smallest value
		else
			Res = VecBtree[i].equal_range( VecBtree[i].rbegin().key());				// the last element holds the biggest value
		TempDotSet.clear();
		while( Res.first != Res.second)
		{
			TempDotSet.insert( Res.first.data());
			++Res.first;
		}

		if( TempDotSet.size() > 1)
		{
			TempComb.insert(TempDotSet);
			TempDotSet.clear();
		}

		if(depth) {
			TempNoeud->depthD = TempDotSet;
			TempNoeud->depthI = TempComb;
		}
		else
			TempNoeud->Fill_D_I( TempDotSet, TempComb);
		PileNm1.push_back( TempNoeud);
	}
	Compteur = NombreDimensions;
	Type1Count = NombreDimensions;
}

void ArbreCube::GenereDimensionProf(const Noeud& ParentNoeud, const DotSet& pcsD, const CombinedSkyline& pcsI, Closure* ParentGroup)
{
#ifdef DEBUG
		std::cout << "Noeud actuel : " << ParentNoeud.Chemin << std::endl;
#endif

	Closure* GroupFound = ParentGroup;
	HashKey ParentKey(pcsD,pcsI);
	if(GroupFound == 0)
		GroupFound = FindClosure(ParentKey, MesNoeudClos);

	// if we computed the last node first, we don't do anything if we reach it in the descent
	if( ParentNoeud.Chemin.size() != static_cast<size_t>(NombreDimensions - 1) || !isLastNodeFirst)
	{
		for( long i = ParentNoeud.offset + 1; i < NombreDimensions; i++)
		{
			Noeud TempNoeud(NombreDimensions-i-1,ParentNoeud.Chemin.size()+1);
			TempNoeud.offset = i;
			TempNoeud.Chemin = ParentNoeud.Chemin;
			TempNoeud.Chemin.push_back(i);
#ifdef DEBUG
				std::cout << "Je suis: " << TempNoeud.Chemin << std::endl;
#endif
#ifdef DEBUG_CLOS
			std::cout << ParentGroup << std::endl;
#endif
			if( GroupFound != 0 && GroupFound->PrunedDimensions[i]) {
				TempNoeud.EstType1 = ParentNoeud.EstType1;
				GenereDimensionProf(TempNoeud, pcsD, pcsI, GroupFound);
			}
			else {
				Noeud* TempNewDim;
				DotSet TempD;
				CombinedSkyline TempI;

				TempNewDim = racine->enfants[i];

				// type 1 ?
				if( ParentNoeud.EstType1)
					ManageType1( &TempNoeud, TempD, TempI, pcsD, pcsI, TempNewDim->depthD, TempNewDim->depthI);

				// type 2
				if( TempNoeud.EstType1 == true)
				{
					NbSkylineFoundDirectly += TempD.size();
					if( TempI.size() == 1)
						NbSkylineFoundDirectly += TempI.begin()->size();
				}
				else
				{
					// step 1
					TempD.insert( pcsD.begin(), pcsD.end());
					TempD.insert( TempNewDim->depthD.begin(), TempNewDim->depthD.end());
					NbSkylineFoundDirectly += TempD.size();

					DotSet bigDotSet;
					const DotSet emptyDotSet;
					std::list<DotSet> Result;
					DotSet::const_iterator iteP = pcsD.begin();
					DotSet::const_iterator iteN = TempNewDim->depthD.begin();
					for( long i = 0; i < NombrePoints; i++)
					{
						if( iteP != pcsD.end() && i == *iteP)
						{
							Result.push_back(emptyDotSet);
							Result.back().insert(i);
							++iteP;
						}
						else
						{
							if( iteN != TempNewDim->depthD.end() && i == *iteN)
							{
								Result.push_back(emptyDotSet);
								Result.back().insert(i);
								++iteN;
							}
							else
								bigDotSet.insert(i);
						}
					}

					BNL( bigDotSet.begin(), bigDotSet.end(), TempNoeud.Chemin, Result);
					for(std::list<DotSet>::const_iterator iteR = Result.begin(); iteR != Result.end(); ++iteR)
					{
						if( (*iteR).size() == 1)
							TempD.insert(*(*iteR).begin());
						else
							TempI.insert(*iteR);
					}
				}

				//TempNoeud->Fill_D_I( TempD, TempI);

				NbSkylineFoundTotal += TempD.size();
				for( CombinedSkyline::const_iterator itI = TempI.begin(); itI != TempI.end(); ++itI)
					NbSkylineFoundTotal += (*itI).size();

				Compteur++;
				if( (Compteur & 0xFFFFFE00) == Compteur)
				{
					std::cout << "\rGenerated " << Compteur << " nodes" << std::endl;
					std::cout.flush();
				}

				GenereDimensionProf(TempNoeud, TempD, TempI, 0);
			}
		}
	}
#ifdef DEBUG
	std::cout << "Je suis: " << ParentNoeud.Chemin << std::endl;
#endif

	// we got a parent group meaning this node is in the same group than its parent
	// so we don't bother trying to add it to a group, it's neither a closure nor a generator
	if(ParentGroup != 0)
		return;

	// maybe this node belongs to a group from one of its legacy
	if(!GroupFound)
		GroupFound = FindClosure(ParentKey, MesNoeudClos);

	if(GroupFound)
		GroupFound->AddElement(ParentNoeud.Chemin);
	else
		MesNoeudClos[ParentKey] = new Closure( ParentNoeud.Chemin, NombreDimensions);
}

void ArbreCube::GenereDimensionN( unsigned long DimNumber)
{
	std::vector<Noeud*>::iterator iteParents;
	std::vector<Noeud*> TempPile;
	TempPile.reserve( MyCnk.Value(DimNumber));

	Noeud* TempNoeud;
	Noeud* TempParent;
	Noeud* TempNewDim;

	ParentsList ListeComposantes;
	ParentsList::const_iterator iteComp;

	UnDotSetList VecUDS;

#ifndef DEBUG
	std::cout << ' ' << (*(PileNm1.begin()))->offset + 2 << std::flush;
#else
	std::cout << "Generating N-Dimension spaces with N = " << (*(PileNm1.begin()))->offset + 2 << std::endl;
#endif // DEBUG

	for( iteParents = PileNm1.begin(); iteParents != PileNm1.end(); ++iteParents)
	{
		for( long i = (*iteParents)->offset + 1; i < NombreDimensions; i++)
		{
			TempNoeud = new Noeud(NombreDimensions-i-1,(*iteParents)->Chemin.size()+1);
			TempNoeud->offset = i;
			TempNoeud->parent = *iteParents;
			(*iteParents)->enfants.push_back( TempNoeud);
			TempNoeud->Chemin = (*iteParents)->Chemin;
			TempNoeud->Chemin.push_back( i);

			TempParent = *iteParents;
			TempNewDim = racine->enfants[i];

			DotSet TempD;
			CombinedSkyline TempI;

			if( TempNoeud->parent->EstComplet)
				TempNoeud->EstComplet = true;
			else
			{
				// type 1 ?
				if( TempNoeud->parent->EstType1)
					ManageType1( TempNoeud, TempD, TempI, TempParent->cs_D, TempParent->cs_I, TempNewDim->cs_D, TempNewDim->cs_I);

				// type 2
				if( TempNoeud->EstType1 == false)
				{
					ListeComposantes.clear();
					CreateParentsList( TempNoeud, ListeComposantes);

					// step 1
					for( iteComp = ListeComposantes.begin(); iteComp != ListeComposantes.end(); ++iteComp)
						TempD.insert( (*iteComp).Parent->cs_D.begin(), (*iteComp).Parent->cs_D.end());

					if( TempD.size() != static_cast<size_t>(NombrePoints))
					{
						// step 2
						Breadth_Step_2_1( VecUDS, ListeComposantes, TempNoeud, TempD);
						Breadth_Step_2_2( VecUDS, TempD, TempI);

						if( TempD.size() != static_cast<size_t>(NombrePoints))
						{
							// step 3
							if( CurrentAlgo == BREADTH)
								RangeBNL( TempNoeud, TempD, TempI);
							else
							{
								assert( CurrentAlgo == BR_DOM);
								TakeTheBus( TempD, TempI, TempNoeud);
							}
						}
					}
				}
			}

			TempNoeud->Fill_D_I( TempD, TempI);

			TempPile.push_back( TempNoeud);
		}
	}

	TempPile.pop_back();	// remove the last element since it's always a leaf
	PileNm1.swap( TempPile);
}

const Noeud* ArbreCube::ComputeLastNode()
{
	Noeud* TempNoeud = new Noeud(0,NombreDimensions);
	TempNoeud->offset = NombreDimensions - 1;
	TempNoeud->parent = 0;
	for( long i = 0; i < NombreDimensions; i++)
		TempNoeud->Chemin.push_back(i);

	DotSet TempDotSet;
	TempDotSet.insert(0);

	TempCombined TempSkyline;
	TempSkyline.push_back( TempDotSet);

	for( long j = 1; j < NombrePoints; j++)
	{
		bool KeepIt = true;

		for( TempCombined::iterator iteTempSk = TempSkyline.begin(); iteTempSk != TempSkyline.end(); )
		{
			PointOrderRelation ResultComp = ComparePoints( j, *((*iteTempSk).begin()), TempNoeud->Chemin);

			if( ResultComp == P2_DOM_P1)
			{
				KeepIt = false;
				break;
			}
			if( ResultComp == EQUIV)
			{
				KeepIt = false;
				(*iteTempSk).insert( j);
				break;
			}
			if( ResultComp == P1_DOM_P2)
				TempSkyline.erase( iteTempSk++);
			else
				++iteTempSk;
		}

		if( KeepIt)
		{
			TempDotSet.clear();
			TempDotSet.insert( j);
			TempSkyline.push_back( TempDotSet);
		}
	}

	TempDotSet.clear();
	CombinedSkyline TempI;

	TempCombined::iterator iteSk = TempSkyline.begin();
	while( iteSk != TempSkyline.end())
	{
		if( (*iteSk).size() == 1)
			TempDotSet.insert( *((*iteSk).begin()));
		else
			TempI.insert( *iteSk);
		++iteSk;
	}

	//TempNoeud->Fill_D_I( TempDotSet, TempI);

	MesNoeudClos[HashKey(TempDotSet, TempI)] = new Closure(TempNoeud->Chemin, NombreDimensions);

	return TempNoeud;
}

void ArbreCube::CreateParentsList( Noeud* TempNoeud, ParentsList& ListeComposantes)
{
	const std::vector<long>& Chemin = TempNoeud->Chemin;

#ifdef DEBUG
	std::vector<long>::const_iterator itc;
	std::cout << "Chemin : ";
	for( itc = Chemin.begin(); itc != Chemin.end(); ++itc)
		std::cout << *itc << ",";
	std::cout << std::endl;
#endif // DEBUG

	UnParent TempUnParent;

	// TODO change naive to smart processing (must be faster recursively)
	for( size_t i = 0; i < Chemin.size(); i++)
	{
		TempUnParent.Parent = racine;
		for( size_t j = 0; j < Chemin.size(); j++)
			if( i != j)
				TempUnParent.Parent = TempUnParent.Parent->enfants[ Chemin[j] - TempUnParent.Parent->offset - 1];
			else
				TempUnParent.RemovedDim = Chemin[j];

		ListeComposantes.push_back( TempUnParent);
	}
}

void ArbreCube::Breadth_Step_2_1(	UnDotSetList& VecUDS,
									const ParentsList& ListeComposantes,
									const Noeud* TempNoeud,
									const DotSet& TempD)
{
	UnDotSet TempUDS;

	VecUDS.clear();
	for( ParentsList::const_iterator iteComp = ListeComposantes.begin(); iteComp != ListeComposantes.end(); ++iteComp)
	{
		for( CompactComb::const_iterator iteSkyline = (*iteComp).Parent->cs_I.begin(); iteSkyline != (*iteComp).Parent->cs_I.end(); ++iteSkyline)
		{
			if( AreDisjoint( TempD.begin(), TempD.end(), iteSkyline->begin(), iteSkyline->end()))
			{
				TempUDS.CombDotSet = &(*iteSkyline);
				TempUDS.RemovedDim = (*iteComp).RemovedDim;
				VecUDS.push_back( TempUDS);
			}
		}
	}
}

void ArbreCube::Breadth_Step_2_2(	const UnDotSetList& VecUDS,
									DotSet& TempD,
									CombinedSkyline& TempI)
{
	double TempMin;
	double TempVal;
	DotSet TempDotSet;
	CompactSet::const_iterator iteDotSet;
	CompactSet ListeMin;				// not a set, just using fast_lloc
	bool MustFindLowest;

	for( UnDotSetList::const_iterator iteVecUDS = VecUDS.begin(); iteVecUDS != VecUDS.end(); ++iteVecUDS)
	{
		iteDotSet = (*iteVecUDS).CombDotSet->begin();
		TempMin = matrice[ *iteDotSet + (*iteVecUDS).RemovedDim * NombrePoints];
		ListeMin.clear();
		ListeMin.push_back( *iteDotSet);
		++iteDotSet;
		MustFindLowest = FindLowest[(*iteVecUDS).RemovedDim];

		// we look for the minimum value on the new dimension
		while( iteDotSet != (*iteVecUDS).CombDotSet->end())
		{
			TempVal = matrice[ *iteDotSet + (*iteVecUDS).RemovedDim * NombrePoints];
			if( (MustFindLowest && TempVal < TempMin) || (!MustFindLowest && TempVal > TempMin))
			{
				TempMin = TempVal;
				ListeMin.clear();
				ListeMin.push_back( *iteDotSet);
			}
			else
			{
				if( TempVal == TempMin)
					ListeMin.push_back( *iteDotSet);
			}
			++iteDotSet;
		}

		if( ListeMin.size() == 1)
			TempD.insert( ListeMin[0]);
		else
		{
			TempDotSet.clear();
			TempDotSet.insert( ListeMin.begin(), ListeMin.end());
			TempI.insert( TempDotSet);
		}
	}
}

void ArbreCube::RangeBNL( const Noeud* TempNoeud, DotSet& TempD, CombinedSkyline& TempI)
{
	DotSet TempDotSet;

	Step_3_1( TempDotSet, TempD, TempI, TempNoeud->Chemin);

#ifdef DEBUG
	std::cout << "\t\t\tEtape 3 : Points finaux: ";
	for( DotSet::const_iterator iteDD = TempDotSet.begin(); iteDD != TempDotSet.end(); ++iteDD)
		std::cout << 'P' << *iteDD + 1 << ',';
	std::cout << std::endl;
#endif // DEBUG

	if( TempDotSet.empty())
		return;

	std::list<DotSet> Stockage;
	DotSet EmptyDotSet;

	for( DotSet::const_iterator iteDs = TempD.begin(); iteDs != TempD.end(); ++iteDs)
		Stockage.insert(Stockage.begin(), EmptyDotSet)->insert( *iteDs);
	for( CombinedSkyline::const_iterator iteDs = TempI.begin(); iteDs != TempI.end(); ++iteDs)
		Stockage.push_back( *iteDs);
	BNL( TempDotSet.begin(), TempDotSet.end(), TempNoeud->Chemin, Stockage);

	Step_3_3( Stockage, TempD, TempI);
}

// TODO This method should never be used again with Depth algorithm, stand-by for conditionals cleanup
void ArbreCube::Step_3_1(	DotSet& TempDotSet,
							const DotSet& TempD,
							const CombinedSkyline& TempI,
							const std::vector<long>& Path)
{
	DotSet::const_iterator iteDotSet;
	CombinedSkyline::const_iterator iteSkyline;
	double TempMin;
	double TempMax;
	double TempVal;
	stx::btree_multimap<double,long>::const_iterator iteLowerBound;
	stx::btree_multimap<double,long>::const_iterator iteUpperBound;
	DotSet InterDotSet;

#ifdef DEBUG
	stx::btree_multimap<double,long>::const_iterator iteLowerBoundDebug;
#endif // DEBUG

	for( std::vector<long>::const_iterator iteChemin = Path.begin(); iteChemin != Path.end(); ++iteChemin)
	{
		long Dimension_Number = *iteChemin;
		iteDotSet = TempD.begin();
		iteSkyline = TempI.begin();

		if( ! TempD.empty())
		{
			TempMin = TempMax = matrice[ Dimension_Number * NombrePoints + *iteDotSet ];
			++iteDotSet;
		}
		else
		{
			TempMin = TempMax = matrice[ Dimension_Number * NombrePoints + *((*iteSkyline).begin()) ];
			++iteSkyline;
		}

		for( ; iteDotSet != TempD.end(); ++iteDotSet)
		{
			TempVal = matrice[ Dimension_Number * NombrePoints + *iteDotSet ];
			if( TempVal < TempMin)
				TempMin = TempVal;
			else
			{
				if( TempVal > TempMax)
					TempMax = TempVal;
			}
		}

		for( ; iteSkyline != TempI.end(); ++iteSkyline)
		{
			TempVal = matrice[ Dimension_Number * NombrePoints + *((*iteSkyline).begin()) ];
			if( TempVal < TempMin)
				TempMin = TempVal;
			else
			{
				if( TempVal > TempMax)
					TempMax = TempVal;
			}
		}

#ifdef DEBUG
		std::cout << "\t\t\tEtape 3 : Dimension " << static_cast<unsigned char>(*iteChemin + 65) << ", min: " << TempMin << ", max: " << TempMax << std::endl;
#endif // DEBUG

		iteLowerBound = VecBtree[Dimension_Number].lower_bound( TempMin);
		iteUpperBound = VecBtree[Dimension_Number].upper_bound( TempMax);

		if( iteChemin == Path.begin())
		{
			TempDotSet.clear();
			for( ; iteLowerBound != iteUpperBound; ++iteLowerBound)
				TempDotSet.insert( iteLowerBound.data());
#ifdef DEBUG
			std::cout << "\t\t\tEtape 3 : Intersection : ";
			for( DotSet::const_iterator iteDD = TempDotSet.begin(); iteDD != TempDotSet.end(); ++iteDD)
				std::cout << 'P' << *iteDD + 1 << ',';
			std::cout << std::endl;
#endif // DEBUG
		}
		else
		{
#ifdef DEBUG
			InterDotSet.clear();
			iteLowerBoundDebug = iteLowerBound;
			for( ; iteLowerBoundDebug != iteUpperBound; ++iteLowerBoundDebug)
				InterDotSet.insert( iteLowerBound.data());

			std::cout << "\t\t\tEtape 3 : Intersection : ";
			for( DotSet::const_iterator iteDD = InterDotSet.begin(); iteDD != InterDotSet.end(); ++iteDD)
				std::cout << 'P' << *iteDD + 1 << ',';
			std::cout << std::endl;
#endif // DEBUG

			if( CurrentAlgo == BREADTH)
				Intersection( TempDotSet, iteLowerBound, iteUpperBound);
			else
			{
				assert( CurrentAlgo == DEPTH);
				Union( TempDotSet, iteLowerBound, iteUpperBound);
			}
		}

		if( CurrentAlgo == BREADTH && TempDotSet.empty())
			break;
	}

	RemoveIncluded( TempDotSet, TempD, TempI);
}


void ArbreCube::Step_3_3( const std::list<DotSet>& Stockage, DotSet& TempD, CombinedSkyline& TempI)
{
	assert( ! Stockage.empty());

	TempD.clear();
	TempI.clear();

	for( std::list<DotSet>::const_iterator iteCSL = Stockage.begin(); iteCSL != Stockage.end(); ++iteCSL)
	{
		if( (*iteCSL).size() == 1)
			TempD.insert( *((*iteCSL).begin()));
		else
			TempI.insert( *iteCSL);
	}
}

PointOrderRelation CPtab[4] = { EQUIV , P1_DOM_P2 , P2_DOM_P1 , UNCOMP };

PointOrderRelation ArbreCube::ComparePoints( long P1, long P2, const std::vector<long>& Chemin)
{
	double v1;
	double v2;
	uint8_t Status = 0;
	std::vector<long>::const_iterator iteChemin;

	for( iteChemin = Chemin.begin(); Status != 3 && iteChemin != Chemin.end(); ++iteChemin)
	{
		v1 = matrice[ *iteChemin * NombrePoints + P1 ];
		v2 = matrice[ *iteChemin * NombrePoints + P2 ];

		if( v1 < v2)
			Status |= FindLowest1[*iteChemin];
		else
			if( v1 > v2)
				Status |= FindLowest2[*iteChemin];
	}

	return CPtab[Status];
}

void ArbreCube::Intersection(	DotSet& sk,
								stx::btree_multimap<double,long>::const_iterator& iteLowerBound,
								stx::btree_multimap<double,long>::const_iterator& iteUpperBound)
{

#ifdef DEBUG
	if( sk.empty())
		throw std::exception();
#endif // DEBUG

	if( iteLowerBound == iteUpperBound)
	{
		sk.clear();
		return;
	}

	std::vector<long> VecPoints;

	while( iteLowerBound != iteUpperBound)
	{
		VecPoints.push_back( iteLowerBound->second);
		++iteLowerBound;
	}
	std::sort( VecPoints.begin(), VecPoints.end());

	DotSet::const_iterator iteSK = sk.begin();
	std::vector<long>::const_iterator iteVP = VecPoints.begin();

	if( *iteSK < *iteVP)
		iteSK = sk.lower_bound( *iteVP);
	sk.erase( sk.begin(), iteSK);

	while( iteSK != sk.end() && iteVP != VecPoints.end())
	{
		if( *iteSK < *iteVP)
			sk.erase(iteSK++);
		else
		{
			if( *iteSK == *iteVP)
				++iteSK;
			++iteVP;;
		}
	}

	if( iteSK != sk.end())
		sk.erase( iteSK, sk.end());
}


void ArbreCube::Union(	DotSet& sk,
						stx::btree_multimap<double,long>::const_iterator& iteLowerBound,
						stx::btree_multimap<double,long>::const_iterator& iteUpperBound)
{
	if( iteLowerBound == iteUpperBound)
		return;

	for( ; iteLowerBound != iteUpperBound; ++iteLowerBound)
		sk.insert( iteLowerBound->second);
}


void ArbreCube::RemoveIncluded( DotSet& InterDotSet,
								const DotSet& TempD,
								const CombinedSkyline& TempI)
{
	DotSet::iterator iteIDS = InterDotSet.begin();
	DotSet::const_iterator iteTN = TempD.begin();

	while( iteIDS != InterDotSet.end() && iteTN != TempD.end())
	{
		if( *iteIDS < *iteTN)
			++iteIDS;
		else
		{
			if( *iteIDS == *iteTN)
				InterDotSet.erase( iteIDS++);
			++iteTN;
		}
	}

	for( CombinedSkyline::const_iterator iteSk = TempI.begin(); iteSk != TempI.end(); ++iteSk)
	{
		for( iteIDS = InterDotSet.begin(); iteIDS != InterDotSet.end();)
		{
			if( (*iteSk).count( *iteIDS) == 1)
				InterDotSet.erase( iteIDS++);
			else
				++iteIDS;
		}
	}
}


void ArbreCube::TakeTheBus( DotSet& TempD, CombinedSkyline& TempI, const Noeud* TempNoeud)
{
	// Get the maximal domain value
	long dim = TempNoeud->Chemin[0];
	double value = domainSize[dim];
	std::vector<long>::const_iterator iteChemin;
	for( iteChemin = TempNoeud->Chemin.begin()+1; iteChemin != TempNoeud->Chemin.end(); ++iteChemin)
	{
		double vtmp = domainSize[*iteChemin];
		if( vtmp >= value)
		{
			value = vtmp;
			dim = *iteChemin;
		}
	}

	// Prepare our list of Skylines already found in the previous steps
	std::multimap<double,long> SP;
	DotSet Omega;
	Flatten( TempD.begin(), TempD.end(), TempI.begin(), TempI.end(), Omega);

	for( DotSet::const_iterator iteDs = Omega.begin(); iteDs != Omega.end(); ++iteDs)
		SP.insert( std::pair<double,long>(Sum(*iteDs,TempNoeud->Chemin), *iteDs));

	const stx::btree_multimap<double,long>& TempBTree = VecBtree[dim];
	long count = 0;

	if( FindLowest[dim])
	{
		for( stx::btree_multimap<double,long>::const_iterator itPoints = TempBTree.begin(); itPoints != TempBTree.end(); ++itPoints)
		{
			// If the element isn't skyline, evaluate it
			if( Omega.find(itPoints.data()) == Omega.end())
				count += Evaluate( itPoints.data(), SP, TempNoeud->Chemin, TempD, TempI);
		}
	}
	else
	{
		for( stx::btree_multimap<double,long>::const_reverse_iterator itPoints = TempBTree.rbegin(); itPoints != TempBTree.rend(); ++itPoints)
		{
			if( Omega.find(itPoints.data()) == Omega.end())
				count += Evaluate( itPoints.data(), SP, TempNoeud->Chemin, TempD, TempI);
		}
	}

#ifdef DEBUG
	/*std::cout << "Final Elements in SP: " << std::endl;
	for (std::multimap<double, long>::iterator it = SP.begin();it != SP.end();++it)
		std::cout << "  [" << (*it).first << ", P" << (*it).second+1 << "]" << std::endl;*/
	std::cout << "Number of comparisons: " << count << std::endl;
	//std::cout << "Final Result for the current node: " << std::endl;
	//AfficheSkyline(TempNoeud,std::cout);
	//std::cout << std::endl;
#endif //DEBUG
}


long ArbreCube::Evaluate(	long NumPoint,
							std::multimap<double,long> &SP,
							const std::vector<long>& Chemin,
							DotSet& TempD,
							CombinedSkyline& TempI )
{
	double Fq = Sum(NumPoint,Chemin);

//#ifdef DEBUG
//	std::cout << "Evaluating point P" << NumPoint+1 << std::endl;
//#endif //DEBUG

	long count = 0;
	for (std::multimap<double, long>::const_iterator it = SP.begin();it != SP.end();++it)
	{
		if(Fq < (*it).first) // Direct insertion as a Skyline
		{
			SP.insert(std::pair<double, long>(Fq,NumPoint));
			// Insert directly in the distinct skyline points
			TempD.insert(NumPoint);
			return 0;
		}
		else
		{
			count++;
			PointOrderRelation ResultComp = ComparePoints((*it).second,NumPoint, Chemin);
			// Domination test p < q or p EQUIV q
			if( ResultComp == P1_DOM_P2)
			{
				//#ifdef DEBUG
				//std::cout << "P" << (*it).second+1 << "<P" << NumPoint+1 << std::endl;
				//#endif //DEBUG
				return count;
			}
			else if (ResultComp == EQUIV)
			{
				//#ifdef DEBUG
				//std::cout << "P" << (*it).second+1 << " EQUIV P" << NumPoint+1 << std::endl;
				//#endif //DEBUG
				for( CombinedSkyline::iterator iteDs = TempI.begin(); iteDs != TempI.end(); ++iteDs)
				{
					DotSet ds = *iteDs;
					if(ds.find((*it).second) != ds.end())
					{
						// FOUND & INSERT
						ds.insert(NumPoint);
						TempI.erase(iteDs);
						TempI.insert(ds);
						SP.insert(std::pair<double, long>(Fq,NumPoint));
						return count;
					}
				}
				//#ifdef DEBUG
				//std::cout << "Creating new non-distinct set for P" << NumPoint+1 << " and P"<< (*it).second+1 << std::endl;
				//#endif //DEBUG

				DotSet ds1;
				ds1.insert((*it).second);
				ds1.insert(NumPoint);
				TempI.insert(ds1);
				// Enlever P du indistinct set
				TempD.erase((*it).second);
				SP.insert(std::pair<double, long>(Fq,NumPoint));
				return count;
			}

		}

	}

//#ifdef DEBUG
//	std::cout << "Inserting to SP at the end" << std::endl;
//	std::cout << "Number of comparison tests " << count << std::endl;
//#endif //DEBUG

	SP.insert(std::pair<double, long>(Fq,NumPoint));
	TempD.insert(NumPoint);
	return count;
}


double ArbreCube::Sum( long NumPoint, const std::vector<long>& Chemin)
{
	double Total = 0;

	for( size_t i = 0; i < Chemin.size(); i++)
	{
		if( FindLowest[Chemin[i]])
			Total += matrice[ Chemin[i] * NombrePoints + NumPoint];
		else
			Total -= matrice[ Chemin[i] * NombrePoints + NumPoint];
	}

	return Total;
}


double ArbreCube::getDomaineSize(int dimension)
{
	double minValue = VecBtree[dimension].begin().key();
	double maxValue = VecBtree[dimension].rbegin().key();

	return maxValue - minValue;
}


Closure* ArbreCube::FindClosure( const HashKey& keyToFind, const HashClosure& VecNoeuds)
{
	if( VecNoeuds.empty())
		return 0;

	HashClosure::const_iterator iteClos = VecNoeuds.find(keyToFind);
	if( iteClos != VecNoeuds.end())
		return iteClos->second;

	return 0;
}


void ArbreCube::AfficheResultat( std::ostream& Cout, std::vector<std::string>* Labels) const
{
	std::vector<Noeud*> Pile;
	for( size_t i = 0; i < racine->enfants.size(); i++)
		AfficheNoeud( racine->enfants[i], Pile, Cout, Labels);

	while( ! Pile.empty())
		AfficheLargeur( Pile, Cout, Labels);
}


void ArbreCube::AfficheLargeur( std::vector<Noeud*>& Pile,
								std::ostream& Cout,
								std::vector<std::string>* Labels) const
{
	std::vector<Noeud*> TempPile;

	for( size_t i = 0; i < Pile.size(); i++)
		AfficheNoeud( Pile[i], TempPile, Cout, Labels);

	Pile.swap(TempPile);
}


void ArbreCube::AfficheNoeud(	const Noeud *TempNoeud,
								std::vector<Noeud*>& Pile,
								std::ostream& Cout,
								std::vector<std::string>* Labels) const
{
	if( TempNoeud == 0)
		return;

	PrintPath( TempNoeud->Chemin, Cout);

	Cout << " : ";
	AfficheSkyline( TempNoeud->cs_D, TempNoeud->cs_I, TempNoeud->EstComplet, Cout, Labels);
	Cout << std::endl;

	for( size_t i = 0; i < TempNoeud->enfants.size(); i++)
		Pile.push_back( TempNoeud->enfants[i]);
}


void ArbreCube::AfficheClos( std::ostream& Cout, std::vector<std::string>* Labels) const
{
	for( HashClosure::const_iterator iteClos = MesNoeudClos.begin(); iteClos != MesNoeudClos.end(); ++iteClos)
	{
		// closures
		std::vector<LatticePath>::const_iterator itClosedNodes = iteClos->second->ClosedNodes.begin();
		PrintPath( *itClosedNodes, Cout);
		for(++itClosedNodes ; itClosedNodes != iteClos->second->ClosedNodes.end(); ++itClosedNodes)
		{
			Cout << ',';
			PrintPath( *itClosedNodes, Cout);
		}

		// skyline
		Cout << " : ";
		AfficheSkyline( iteClos->first.first, iteClos->first.second, false, Cout, Labels);

		// generators
		Cout << " : ";
		std::vector<LatticePath>::const_iterator itgen = iteClos->second->Generators.begin();
		PrintPath( *itgen, Cout);
		for( ++itgen; itgen != iteClos->second->Generators.end(); ++itgen)
		{
			Cout << ',';
			PrintPath( *itgen, Cout);
		}

		Cout << std::endl;
	}
}


std::ostream& ArbreCube::PrintPath( const std::vector<long>& Chemin, std::ostream& Cout)
{
	for( std::vector<long>::const_iterator reve = Chemin.begin(); reve != Chemin.end(); ++reve)
		Cout << 'd' << (*reve);
	return Cout;
}


void ArbreCube::HashStat(std::ostream& Cout) const
{
	Cout << MesNoeudClos.bucket_count() << " buckets :";
	for( size_t nbBuckets = 0; nbBuckets < MesNoeudClos.bucket_count(); ++nbBuckets)
		Cout << " " << MesNoeudClos.bucket_size(nbBuckets);
	Cout << std::endl;
}


size_t ArbreCube::GetNbClos() const {
	size_t total = 0;
	for( HashClosure::const_iterator iteL = MesNoeudClos.begin(); iteL != MesNoeudClos.end(); ++iteL) {
		total += iteL->second->ClosedNodes.size();
	}
	return total;
}


