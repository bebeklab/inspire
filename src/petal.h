#ifndef PETAL_H
#define PETAL_H
#include "base_header.h"
#include "gene.h"
//#include <boost/format.hpp>
//#include <boost/lexical_cast.hpp>

#include <sqlite3.h> 

class Petal
{
public:
	bool GenerateGO;
	static sqlite3* con;
		
	std::vector<Gene*> GenesList;
	//GeneMap gives us the index of the gene in the GenesList and in the Matrix. Perfect way to look up a gene in cosntant time
	std::map<std::string, int> GenesMap;
	string Name;
	int size;
	int interactionSize;
	//Adjacency matrix.  0 is no interaction, 1 is predicted, 2 is confirmed
	int ** matrix;
	Petal(std::string Name, int);
	~Petal(); // destructor function
	int AddGene(string Name);
	int getGeneIndex(string Name);
	void AddInteraction(std::string Name1, std::string Name2, string InteractionType);
	void PrintPetalGenes();
	void PrintPetalInteractions();
	int getSize(){return size;};
	static int LoadPetals(std::string PetalAName, std::string PetalBName ,  vector<Petal*> *Petals ) ;
	int checkInteract(std::string Name1, std::string Name2);

	

	static int OpenDatabase(
	    const std::string& databaseFile
	);

	static void CloseDatabase();
	
	//static vector<string> trimGo( vector<string>);
	static vector<string> removeDuplicatesFromStringVector( vector<string> vec);
	static void printStringVector( vector<string> vec );
	
	
	
	
};
#endif // PETAL_H
