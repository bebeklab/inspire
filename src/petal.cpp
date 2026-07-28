#include "petal.h"

#include <sqlite3.h>

#include <cstdlib>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

using namespace std;

sqlite3* Petal::con = NULL;


namespace
{

	void CheckSQLiteResult(
	    int result,
	    sqlite3* database,
	    const string& operation)
	{
	    if (result != SQLITE_OK)
	    {
	        string message = operation;
	
	        if (database != NULL)
	        {
	            message += ": ";
	            message += sqlite3_errmsg(database);
	        }
	
	        throw runtime_error(message);
	    }
	}
	
	void BindText(
	    sqlite3_stmt* statement,
	    int parameter,
	    const string& value)
	{
	    int result = sqlite3_bind_text(
	        statement,
	        parameter,
	        value.c_str(),
	        -1,
	        SQLITE_TRANSIENT
	    );
	
	    if (result != SQLITE_OK)
	    {
	        throw runtime_error(
	            "Could not bind text to SQLite statement"
	        );
	    }
	}
	
	string ColumnText(
	    sqlite3_stmt* statement,
	    int column)
	{
	    const unsigned char* value =
	        sqlite3_column_text(statement, column);
	
	    if (value == NULL)
	    {
	        return "";
	    }
	
	    return reinterpret_cast<const char*>(value);
	}
	
	sqlite3_stmt* PrepareStatement(
	    sqlite3* database,
	    const string& sql)
	{
	    sqlite3_stmt* statement = NULL;
	
	    int result = sqlite3_prepare_v2(
	        database,
	        sql.c_str(),
	        -1,
	        &statement,
	        NULL
	    );
	
	    CheckSQLiteResult(
	        result,
	        database,
	        "Could not prepare SQLite statement"
	    );
	
	    return statement;
	}

}



Petal::Petal(std::string Name, int SizeOfPetal){
	
	printf("new petal %s (%d)\n", Name.c_str(),SizeOfPetal);
	this->GenerateGO = false;
	this->size = SizeOfPetal;
	this->interactionSize=0;
	this->GenesList = vector<Gene *>();
//	this->GenesMap;
	this->Name = Name;
	
	if (SizeOfPetal <= 0){ 
		this->matrix = NULL;
		return;
	}

	
	//Allocate matrix
	//this->matrix = (int**) malloc( sizeof(int*) * SizeOfPetal);
	// Allocate and initialize the interaction matrix.
	this->matrix = static_cast<int**>(malloc(sizeof(int*) * SizeOfPetal));
	if (this->matrix == NULL){
		throw bad_alloc();
	}

	for (int i = 0; i < SizeOfPetal; i++)
	{
	    this->matrix[i] =
	        static_cast<int*>(
	            calloc(
	                SizeOfPetal,
	                sizeof(int)
	            )
	        );

	    if (this->matrix[i] == NULL)
	    {
	        for (int j = 0; j < i; j++)
	        {
	            free(this->matrix[j]);
	        }

	        free(this->matrix);
	        this->matrix = NULL;

	        throw bad_alloc();
	    }
	}
 
}

Petal::~Petal(){
	for (int i =0; i < this->size; i++ ){
		free(this->matrix[i]);
	}
	free(this->matrix );
};

int Petal::getGeneIndex(std::string Name){
	
	return this->GenesMap[Name];
	
}


int Petal::checkInteract(std::string Name1, std::string Name2){
	
	return this->matrix[ this->GenesMap[Name1] ][this->GenesMap[Name2] ];
	
}


//Return index
int Petal::AddGene(std::string Name)
{
	if ( this->GenesMap.count(Name) == 0 ) {
		Gene * gene = new Gene( Name, this);
		this->GenesList.push_back( gene );
		int index = this->GenesList.size() - 1;
		gene->GenesListIndex=index;
		this->GenesMap[Name] = index;
	}
	return this->GenesMap[Name];
}

void Petal::AddInteraction(std::string Name1, std::string Name2, string InteractionType)
{
	int i = this->AddGene( Name1);
	int j = this->AddGene( Name2);
	
	//Add to the matrix
	if ( InteractionType == "pp"){
		this->matrix[i][j] = 2;	
		this->matrix[j][i] = 2;	
	}else if (  InteractionType == "pd" ) {
		this->matrix[i][j] = 1;
		this->matrix[j][i] = 1;
	}
	//cout << "Adding Interaction " << i<< " " << j << " "<< " of type" << this->matrix[i][j]<<endl;
	
}
void Petal::PrintPetalInteractions(){
	int counter=0;
	cout<<CYAN<<"Petal Interactions"<<green<<endl;
	map<std::string,int>::iterator it;
	map<std::string,int>::iterator it2;
	for ( it=(this->GenesMap).begin() ; it != (this->GenesMap).end(); it++ ){
		for ( it2=(this->GenesMap).begin() ; it2 != (this->GenesMap).end(); it2++ ){
			if((*it).second>=(*it2).second){
				if (this->matrix[(*it).second][(*it2).second]==1){
					cout<< ++counter <<" " << (*it).first << "\tpd\t"<<  (*it2).first <<endl;
				}
				else if(this->matrix[(*it).second][(*it2).second]==2){
					cout<< ++counter <<" " << (*it).first << "\tpp\t"<<  (*it2).first <<endl;
				}
			}
		}
	}
	cout<< NC<< endl;
}



void Petal::PrintPetalGenes (){
	cout<<GREEN<<"Petal Name: "<<YELLOW<<this->Name<<endl<<GREEN<<"Petal Size: "<<YELLOW<<(this->size)<<red<<endl;

	for (vector<Gene *>::iterator it = this->GenesList.begin(); it!=this->GenesList.end(); ++it) {
		cout << (*it)->getName() << endl;
	}
	cout<<NC<<"---- "<<endl;

}


int Petal::OpenDatabase(
    const std::string& databaseFile)
{
    if (Petal::con != NULL)
    {
        sqlite3_close(Petal::con);
        Petal::con = NULL;
    }

    int result =
        sqlite3_open_v2(
            databaseFile.c_str(),
            &Petal::con,
            SQLITE_OPEN_READONLY,
            NULL
        );

    if (result != SQLITE_OK)
    {
        cerr
            << "Could not open SQLite database: "
            << databaseFile
            << endl;

        if (Petal::con != NULL)
        {
            cerr
                << sqlite3_errmsg(Petal::con)
                << endl;

            sqlite3_close(Petal::con);
            Petal::con = NULL;
        }

        return -1;
    }

    sqlite3_busy_timeout(
        Petal::con,
        30000
    );

    return 0;
}

void Petal::CloseDatabase()
{
    if (Petal::con != NULL)
    {
        sqlite3_close(Petal::con);
        Petal::con = NULL;
    }
}

int Petal::LoadPetals(
    std::string PetalAName,
    std::string PetalBName,
    vector<Petal*>* Petals)
{
    if (Petal::con == NULL)
    {
        cerr << "SQLite database is not open." << endl;
        return -1;
    }

    if (Petals == NULL || Petals->size() < 2)
    {
        cerr << "Petals vector must have two positions." << endl;
        return -1;
    }

    string petalNames[2];

    petalNames[0] = PetalAName;
    petalNames[1] = PetalBName;

    try
    {
        cout << "Begin PETAL setup" << endl;

        /*
         * Process the two requested petals.
         */
        for (int petalNumber = 0;
             petalNumber < 2;
             petalNumber++)
        {
            const string& petalName =
                petalNames[petalNumber];

            /*
             * ---------------------------------------------------------
             * 1. Get the number of genes in this petal.
             * ---------------------------------------------------------
             *
             * The Python database builder stored this number in
             * petal.gene_count.
             */
            sqlite3_stmt* countStatement =
                PrepareStatement(
                    Petal::con,
                    "SELECT gene_count "
                    "FROM petal "
                    "WHERE petal = ?1"
                );

            BindText(
                countStatement,
                1,
                petalName
            );

            int result =
                sqlite3_step(countStatement);

            if (result != SQLITE_ROW)
            {
                sqlite3_finalize(countStatement);

                cerr
                    << "Petal was not found in SQLite: "
                    << petalName
                    << endl;

                return -1;
            }

            int geneCount =
                sqlite3_column_int(
                    countStatement,
                    0
                );

            sqlite3_finalize(countStatement);

            if (geneCount <= 0)
            {
                cerr
                    << "Petal has no genes: "
                    << petalName
                    << endl;

                return -1;
            }

            /*
             * Allocate the Petal and its interaction matrix.
             */
            (*Petals)[petalNumber] =
                new Petal(
                    petalName,
                    geneCount
                );

            Petal* ptrPetal =
                (*Petals)[petalNumber];

            /*
             * ---------------------------------------------------------
             * 2. Add every gene belonging to the petal.
             * ---------------------------------------------------------
             *
             * We add genes before adding interactions. This ensures that
             * GenesList, GenesMap, and matrix indexes are stable.
             */
            sqlite3_stmt* geneStatement =
                PrepareStatement(
                    Petal::con,
                    "SELECT gene_symbol "
                    "FROM petal_gene "
                    "WHERE petal = ?1 "
                    "ORDER BY gene_symbol"
                );

            BindText(
                geneStatement,
                1,
                petalName
            );

            while (
                (result = sqlite3_step(geneStatement))
                == SQLITE_ROW)
            {
                string geneName =
                    ColumnText(
                        geneStatement,
                        0
                    );

                ptrPetal->AddGene(
                    geneName
                );
            }

            if (result != SQLITE_DONE)
            {
                string errorMessage =
                    sqlite3_errmsg(Petal::con);

                sqlite3_finalize(geneStatement);

                throw runtime_error(
                    "Error loading genes for petal " +
                    petalName +
                    ": " +
                    errorMessage
                );
            }

            sqlite3_finalize(geneStatement);

            /*
             * The number read from petal_gene must match petal.gene_count.
             */
            if (
                static_cast<int>(
                    ptrPetal->GenesList.size()
                )
                != geneCount)
            {
                throw runtime_error(
                    "Gene count does not match for petal " +
                    petalName
                );
            }

            /*
             * ---------------------------------------------------------
             * 3. Load the unique petal edges.
             * ---------------------------------------------------------
             *
             * The new database stores all networks in one edge table.
             * The petal column identifies the network.
             */
            sqlite3_stmt* edgeStatement =
                PrepareStatement(
                    Petal::con,
                    "SELECT "
                    "    gene1, "
                    "    gene2, "
                    "    interaction_type "
                    "FROM edge "
                    "WHERE petal = ?1 "
                    "ORDER BY gene1, gene2"
                );

            BindText(
                edgeStatement,
                1,
                petalName
            );

            int interactionCounter = 0;

            while (
                (result = sqlite3_step(edgeStatement))
                == SQLITE_ROW)
            {
                string gene1 =
                    ColumnText(
                        edgeStatement,
                        0
                    );

                string gene2 =
                    ColumnText(
                        edgeStatement,
                        1
                    );

                string interactionType =
                    ColumnText(
                        edgeStatement,
                        2
                    );

                ptrPetal->AddInteraction(
                    gene1,
                    gene2,
                    interactionType
                );

                interactionCounter++;
            }

            if (result != SQLITE_DONE)
            {
                string errorMessage =
                    sqlite3_errmsg(Petal::con);

                sqlite3_finalize(edgeStatement);

                throw runtime_error(
                    "Error loading edges for petal " +
                    petalName +
                    ": " +
                    errorMessage
                );
            }

            sqlite3_finalize(edgeStatement);

            ptrPetal->interactionSize =
                interactionCounter;

            cout
                << "# "
                << interactionCounter
                << " interactions in "
                << petalName
                << endl;

            /*
             * ---------------------------------------------------------
             * 4. Load the preferred UniProt accession for each gene.
             * ---------------------------------------------------------
             *
             * Your Gene class stores UniProt accessions in:
             *
             *     vector<string> Uniprot
             *
             * The new database normally has exactly one preferred
             * accession for each gene.
             */
            for (
                vector<Gene*>::iterator geneIterator =
                    ptrPetal->GenesList.begin();
                geneIterator !=
                    ptrPetal->GenesList.end();
                ++geneIterator)
            {
                Gene* gene =
                    *geneIterator;

                sqlite3_stmt* uniprotStatement =
                    PrepareStatement(
                        Petal::con,
                        "SELECT uniprot_accession "
                        "FROM gene_uniprot "
                        "WHERE gene_symbol = ?1 "
                        "  AND is_preferred = 1"
                    );

                BindText(
                    uniprotStatement,
                    1,
                    gene->Name
                );

                result =
                    sqlite3_step(
                        uniprotStatement
                    );

                if (result != SQLITE_ROW)
                {
                    sqlite3_finalize(
                        uniprotStatement
                    );

                    throw runtime_error(
                        "No preferred UniProt accession for gene " +
                        gene->Name
                    );
                }

                string uniprotAccession =
                    ColumnText(
                        uniprotStatement,
                        0
                    );

                gene->Uniprot.push_back(
                    uniprotAccession
                );

                /*
                 * Verify that the gene has only one preferred accession.
                 */
                result =
                    sqlite3_step(
                        uniprotStatement
                    );

                if (result == SQLITE_ROW)
                {
                    sqlite3_finalize(
                        uniprotStatement
                    );

                    throw runtime_error(
                        "Multiple preferred UniProt accessions for gene " +
                        gene->Name
                    );
                }

                if (result != SQLITE_DONE)
                {
                    string errorMessage =
                        sqlite3_errmsg(Petal::con);

                    sqlite3_finalize(
                        uniprotStatement
                    );

                    throw runtime_error(
                        "Error loading UniProt accession for gene " +
                        gene->Name +
                        ": " +
                        errorMessage
                    );
                }

                sqlite3_finalize(
                    uniprotStatement
                );

                /*
                 * -----------------------------------------------------
                 * 5. Load Pfam domains using the preferred accession.
                 * -----------------------------------------------------
                 */
                sqlite3_stmt* pfamStatement =
                    PrepareStatement(
                        Petal::con,
                        "SELECT pfam_accession "
                        "FROM protein_pfam "
                        "WHERE uniprot_accession = ?1 "
                        "ORDER BY pfam_accession"
                    );

                BindText(
                    pfamStatement,
                    1,
                    uniprotAccession
                );

                while (
                    (result =
                        sqlite3_step(
                            pfamStatement
                        ))
                    == SQLITE_ROW)
                {
                    string pfamAccession =
                        ColumnText(
                            pfamStatement,
                            0
                        );

                    gene->AddPfam(
                        pfamAccession
                    );
                }

                if (result != SQLITE_DONE)
                {
                    string errorMessage =
                        sqlite3_errmsg(Petal::con);

                    sqlite3_finalize(
                        pfamStatement
                    );

                    throw runtime_error(
                        "Error loading Pfam domains for gene " +
                        gene->Name +
                        ": " +
                        errorMessage
                    );
                }

                sqlite3_finalize(
                    pfamStatement
                );
            }

            /*
             * GO terms are no longer loaded into Gene::GO here.
             *
             * The Python program already calculated the pairwise GO
             * scores and stored them in go_pair_score.
             */
            ptrPetal->GenerateGO = false;

            cout
                << "Done adding UniProt and Pfam annotations for "
                << petalName
                << "."
                << endl;
        }

        cout
            << "Petal Sizes: "
            << (*Petals)[0]->size
            << ", "
            << (*Petals)[1]->size
            << endl;

        /*
         * Keep your existing diagnostic output.
         */
        (*Petals)[0]->PrintPetalGenes();
        (*Petals)[0]->PrintPetalInteractions();

        (*Petals)[1]->PrintPetalGenes();
        (*Petals)[1]->PrintPetalInteractions();

        cout
            << "Done creating Petals."
            << endl;

        return
            (*Petals)[0]->size
            *
            (*Petals)[1]->size;
    }
    catch (const exception& error)
    {
        cerr
            << "Could not load petals: "
            << error.what()
            << endl;

        /*
         * Clean up any Petal objects that were created before
         * the error occurred.
         */
        if ((*Petals)[0] != NULL)
        {
            delete (*Petals)[0];
            (*Petals)[0] = NULL;
        }

        if ((*Petals)[1] != NULL)
        {
            delete (*Petals)[1];
            (*Petals)[1] = NULL;
        }

        return -1;
    }
}


vector<string> Petal::removeDuplicatesFromStringVector( vector<string> vec) {

	std::map<std::string, int> stringCount;
	
	vector<string> toReturn;
	for (vector<string>::iterator it = vec.begin(); it<vec.end(); ++it) {
		if ( stringCount.count( (*it) ) == 0) {
			toReturn.push_back( (*it));
			stringCount[ (*it) ] = 1;	
		}
	}
	return toReturn;

}

void Petal::printStringVector( vector<string> vec ) {


	cout <<RED<<"Begin printing vector"<<NC <<endl;
	for (vector<string>::iterator it = vec.begin(); it<vec.end(); ++it) {
		cout << (*it)<<endl;
	}
	cout <<RED<<"End printing vector"<<NC <<endl;

}


