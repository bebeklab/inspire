#ifndef GENE_H
#define GENE_H

#include "base_header.h"



class Petal;
class Gene
{
public:
    std::string Name;
	//Petal it belongs to
	Petal * ParentPetal;
	std::vector<string> GO;
	std::vector<string> Pfam;
	std::vector<string> Uniprot; 
    Gene(string Name,Petal * const petal);
	void AddGO(string);
	void AddPfam(string);
	void AddUniprot(string);
	std::string  getName(){return this->Name;}; 
	int GenesListIndex;

	// Operator Overloading 
    // Comparison Operators == and !=
	bool operator==( Gene& other)  ;
	bool operator!=( Gene& other)  ;
};

#endif // GENE_H
