#ifndef ALIGNEDGENEEDGE_H
#define ALIGNEDGENEEDGE_H

#include "base_header.h"
#include "gene.h"
#include "petal.h"
#include <map>
#include "time.h"
class AlignedGeneEdge
{
public:

	static map<  pair < pair < string, string>, pair <string, string> >,float > WeightCache;
	static map <  pair<string,string>,float> GOWeightCache;

	static vector<AlignedGeneEdge*> AllEdges;

	//A means petalA.  A1 is interactor 1 in petal A.  
	Gene * A1;
	Gene * A2;
	Gene * B1;
	Gene * B2;
	int Weight;
	AlignedGeneEdge(Gene* petalAgene1, Gene* petalAgene2, Gene* petalBgene1, Gene* petalBgeneB2);
	
	
	float SetWeight();
	int adjIndex;
	static double InformationContentOfGenes( Gene*, Gene*);
	static void LoadGOWeightCache();
	static void PrintTime();

    static void SetGOAlgorithm(const string& algorithm);

    static double GetCachedGOScore(Gene* geneA, Gene* geneB);
	static double LoadGONormalization();


    static string GOAlgorithm;
    static double GONormalization;
	static double GONormalize;
};
#endif
