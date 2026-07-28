#include "gene.h"

Gene::Gene(string Name,Petal * const petal)
{
	//This is for the sync'd depth first traversal
        this->Name = Name;
        this->ParentPetal = petal;
        this->GO =  vector<string>();
        this->Pfam = vector<string>();
		this->Uniprot = vector<string>();
		
		this->GenesListIndex = -1;
				
}


void Gene::AddGO( string term) {
	//cout << "Adding GO term " << term << " to " << this->Name << endl;
	this->GO.push_back(term);
}
void Gene::AddPfam( string term) {
	//cout << "Adding Pfam term" << this->Name << endl;
	this->Pfam.push_back(term);
}


void Gene::AddUniprot(string preferredAccession){
	this->Uniprot.push_back(preferredAccession);
}

// define overloaded  operator are defined  
bool Gene::operator==( Gene& other) {
// with string
//	return ((this->Name).compare(other.Name));
// with index
	return (this->GenesListIndex==other.GenesListIndex);
	
}
bool Gene::operator!=( Gene& other) {
	return !( (*this) == (other));
}	
