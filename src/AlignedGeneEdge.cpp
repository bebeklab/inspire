#include "AlignedGeneEdge.h"
#include "base_header.h"
#include <map>
#include <iostream> 
#include <math.h>
using namespace std;

map<  pair < pair < string, string>, pair <string, string> >, float > AlignedGeneEdge::WeightCache;
// map<  string,int > AlignedGeneEdge::GenesAssociated;
map <  pair<string,string>,float> AlignedGeneEdge::GOWeightCache;
vector<AlignedGeneEdge*> AlignedGeneEdge::AllEdges;

string AlignedGeneEdge::GOAlgorithm = "inspire_legacy_v1";
double AlignedGeneEdge::GONormalize = 1.0;
//double AlignedGeneEdge::GONormalize = 23.7367;
extern time_t start_time;
extern int VERBOSE;
 

void AlignedGeneEdge::SetGOAlgorithm(
    const string& algorithm)
{
    if (
        algorithm != "inspire_legacy_v1" &&
        algorithm != "inspire_legacy_exact_v1" &&
        algorithm != "resnik_bp_v1")
    {
        throw runtime_error(
            "Unknown GO scoring algorithm: " + algorithm
        );
    }

    GOAlgorithm = algorithm;
}


AlignedGeneEdge::AlignedGeneEdge(Gene* petalAgene1, Gene* petalAgene2, Gene* petalBgene1, Gene* petalBgene2){
	this->A1 = petalAgene1;
	this->A2 = petalAgene2;
	this->B1 = petalBgene1;
	this->B2 = petalBgene2;
	this->Weight = 0;
	this->adjIndex=-1;
}

double AlignedGeneEdge::GetCachedGOScore(
    Gene* geneA,
    Gene* geneB)
{
    if (Petal::con == NULL)
    {
        throw runtime_error(
            "SQLite database is not open"
        );
    }

    if (geneA == NULL || geneB == NULL)
    {
        throw runtime_error(
            "Null Gene passed to GetCachedGOScore"
        );
    }

    if (geneA->Uniprot.empty())
    {
        throw runtime_error(
            "No UniProt accession loaded for gene " +
            geneA->Name
        );
    }

    if (geneB->Uniprot.empty())
    {
        throw runtime_error(
            "No UniProt accession loaded for gene " +
            geneB->Name
        );
    }

    string protein1 =
        geneA->Uniprot[0];

    string protein2 =
        geneB->Uniprot[0];

    if (protein2 < protein1)
    {
        string temporary =
            protein1;

        protein1 =
            protein2;

        protein2 =
            temporary;
    }

    sqlite3_stmt* statement =
        NULL;

    const char* sql =
        "SELECT score "
        "FROM go_pair_score "
        "WHERE protein1 = ?1 "
        "  AND protein2 = ?2 "
        "  AND algorithm = ?3 "
        "LIMIT 1";

    int result =
        sqlite3_prepare_v2(
            Petal::con,
            sql,
            -1,
            &statement,
            NULL
        );

    if (result != SQLITE_OK)
    {
        throw runtime_error(
            string(
                "Could not prepare GO score query: "
            ) +
            sqlite3_errmsg(Petal::con)
        );
    }

    result =
        sqlite3_bind_text(
            statement,
            1,
            protein1.c_str(),
            -1,
            SQLITE_TRANSIENT
        );

    if (result != SQLITE_OK)
    {
        sqlite3_finalize(statement);

        throw runtime_error(
            "Could not bind first UniProt accession"
        );
    }

    result =
        sqlite3_bind_text(
            statement,
            2,
            protein2.c_str(),
            -1,
            SQLITE_TRANSIENT
        );

    if (result != SQLITE_OK)
    {
        sqlite3_finalize(statement);

        throw runtime_error(
            "Could not bind second UniProt accession"
        );
    }

    result =
        sqlite3_bind_text(
            statement,
            3,
            GOAlgorithm.c_str(),
            -1,
            SQLITE_TRANSIENT
        );

    if (result != SQLITE_OK)
    {
        sqlite3_finalize(statement);

        throw runtime_error(
            "Could not bind GO algorithm"
        );
    }

    result =
        sqlite3_step(statement);

    if (result != SQLITE_ROW)
    {
        string errorMessage =
            "Missing cached GO score for " +
            geneA->Name +
            " (" +
            protein1 +
            ") and " +
            geneB->Name +
            " (" +
            protein2 +
            "), algorithm " +
            GOAlgorithm;

        sqlite3_finalize(statement);

        throw runtime_error(
            errorMessage
        );
    }

    double score =
        sqlite3_column_double(
            statement,
            0
        );

    sqlite3_finalize(statement);

    return score;
}

// returns true if left<right
bool pfamStringCompare( const string &left, const string &right ){
	for( string::const_iterator lit = left.begin(), rit = right.begin(); lit != left.end() && rit != right.end(); ++lit, ++rit )
        if( tolower( *lit ) < tolower( *rit ) )
           return true;
        else if( tolower( *lit ) > tolower( *rit ) )
           return false;
     if( left.size() < right.size() )
        return true;
     return false;
}
 
//void AlignedGeneEdge::printAGE(){
//	cout<< "{" << this->A1  << "  " << this->B1 << "}"<< endl 
//		<< "{" << this->A2  << "  " << this->B2 << "}"<< endl ;
//}
//
//void printAGE(AlignedGeneEdge* A){
//	cout<< "{" << A->A1  << "  " << A->B1 << "}"<< endl 
//		<< "{" << A->A2  << "  " << A->B2 << "}"<< endl ;
//}
//
//
// void printAGE(AlignedGeneEdge* A,AlignedGeneEdge* B){
//	cout<< endl;
//	//"{" << A->A1->Name  << "  " << A->B1->Name << "}"<< " --- "
//	//	<< "{" << B->A1->Name  << "  " << B->B1->Name << "}"<< endl 
//	//	<< "{" << A->A2->Name  << "  " << A->B2->Name << "}"<< " --- "
//	//	<< "{" << B->A2->Name  << "  " << B->B2->Name << "}"<< endl ;
// }


float AlignedGeneEdge::SetWeight( ) {
	
	if (AlignedGeneEdge::WeightCache.count( std::make_pair( std::make_pair(this->A1->Name, this->A2->Name), std::make_pair(this->B1->Name, this->B2->Name))) == 0 ) {
		//cout  << "Weighting:" << this->A1->Name << "-"<< this->B1->Name << "---" << this->Weight << "---" << this->A2->Name << "--" << this->B2->Name <<  endl;
		//calculate the weight for the two gene alignement
 
 		float go_weight = 0;
		float A1B1_score = 0;
		float A2B2_score = 0;
		float pfam_score = 0;
		float unioncount = 0;
		float intercount = 0;
		int totalCount=(this->A1->Pfam).size()+(this->B1->Pfam).size();
		vector<string>::iterator it;
		vector<string>::iterator it2;
				
		if ( (this->A1->Pfam).size() !=0    ||  (this->B1->Pfam).size() !=0){
			// a1b1 score is the intesection of pfams over union of pfam terms.
			//A1B1_score     = len(( A1_pfam &  B1_pfam)) / float(len ( A1_pfam | B1_pfam ))
			//cout << "Begin sorting Pfam Vectors" << endl;
			sort( (this->A1->Pfam).begin(), (this->A1->Pfam).end(), pfamStringCompare );
			sort( (this->B1->Pfam).begin(), (this->B1->Pfam).end(), pfamStringCompare );
			//cout << "End sorting Pfam Vectors" << endl;
			
			std::vector<string> mergedPfam(totalCount);
			unioncount=0;
			intercount=0;
			merge( (this->A1->Pfam).begin(), (this->A1->Pfam).end(),
				   (this->B1->Pfam).begin(), (this->B1->Pfam).end(), mergedPfam.begin());
			//cout << "End merging Pfam Vectors" << endl;
			// this cleans up the vector
			it2 = unique (mergedPfam.begin(), mergedPfam.end()); // 
		  	mergedPfam.resize( it2 - mergedPfam.begin() );       // 
			
			unioncount=mergedPfam.size();
			intercount=totalCount-unioncount;
			A1B1_score=float(intercount/unioncount);
			//cout << "A1B1_Score: " << A1B1_score << endl;
		}
		
		totalCount=(this->A2->Pfam).size()+(this->B2->Pfam).size();
		//cout << "Begin Repeat for Child node" << endl;
		if ( (this->A2->Pfam).size() !=0    ||  (this->B2->Pfam).size() !=0){
			// a1b1 score is the intesection of pfams over union of pfam terms.
			//A1B1_score     = len(( A1_pfam &  B1_pfam)) / float(len ( A1_pfam | B1_pfam ))
			//cout << "Begin sorting Pfam Vectors" << endl;
			sort( (this->A2->Pfam).begin(), (this->A2->Pfam).end(), pfamStringCompare );
			sort( (this->B2->Pfam).begin(), (this->B2->Pfam).end(), pfamStringCompare );
			//cout << "End sorting Pfam Vectors" << endl;
			
			std::vector<string> mergedPfam(totalCount);
			unioncount=0;
			intercount=0;
			merge( (this->A2->Pfam).begin(), (this->A2->Pfam).end(),
				   (this->B2->Pfam).begin(), (this->B2->Pfam).end(), mergedPfam.begin());
			//cout << "End merging Pfam Vectors" << endl;
  			// using predicate comparison:
			// this cleans up the vector
			it2 = unique (mergedPfam.begin(), mergedPfam.end()); // 
		  	mergedPfam.resize( it2 - mergedPfam.begin() );       // 
			
			unioncount=mergedPfam.size();
			intercount=totalCount-unioncount;
			A2B2_score=float(intercount/unioncount);
			//cout << "A2B2_Score: " << A2B2_score << endl;
			
		}
		pfam_score=A1B1_score*A2B2_score;
		//cout << "Done repeat for Child node: Total pfam score: " << pfam_score << endl;
		//This is the new GO way of doing it

		double InformationContentOf1 =
		    AlignedGeneEdge::InformationContentOfGenes(
		        this->A1,
		        this->B1
		    )
		    /
		    AlignedGeneEdge::GONormalize;

		double InformationContentOf2 =
		    AlignedGeneEdge::InformationContentOfGenes(
		        this->A2,
		        this->B2
		    )
		    /
		    AlignedGeneEdge::GONormalize;


		//cout<< "InformationContentOf1: "<<InformationContentOf1<<endl;
		//cout<< "InformationContentOf2: "<<InformationContentOf2<<endl;
		//go_weight = (InformationContentOf1 + InformationContentOf2)/2;
		go_weight = (InformationContentOf1*InformationContentOf2);

		
		
		//cout << this->A1->Name << "-" << this->A2->Name <<endl;
		//cout << this->B1->Name << "-" << this->B2->Name <<endl;
		//cout << "GO:" << go_weight << endl;
		//cout << "pfam:" << pfam_score << endl;
		this->Weight=(int)(100*(go_weight+pfam_score));		
		//this->WeightCalculated=true;		

// comment next
	//	cout<< this->A1->Name << "~"<< this->B1->Name << "---" << this->Weight << "---" << this->A2->Name << "~" << this->B2->Name <<  endl;


		AlignedGeneEdge::WeightCache[ std::make_pair( std::make_pair(this->A1->Name, this->A2->Name), std::make_pair(this->B1->Name, this->B2->Name))] = this->Weight;
		
		
		// TODO ... add sql for fig 6...
		
		//this->Weight = 2.0;
		//this->WeightCalculated = true;
		//cout << "WE dib;t have it!" <<endl;
	}else {
		//cout << "WE have it!" << this->Weight <<endl;
		this->Weight = AlignedGeneEdge::WeightCache[ std::make_pair( std::make_pair(this->A1->Name, this->A2->Name), std::make_pair(this->B1->Name, this->B2->Name))];
		
	}
	// this->IndividualWeight = AlignedGeneEdge::WeightCache[ std::make_pair( std::make_pair(this->A1->Name, this->A2->Name), std::make_pair(this->B1->Name, this->B2->Name))];
	return this->Weight;

}

 
void AlignedGeneEdge::LoadGOWeightCache()
{
    if (Petal::con == NULL)
    {
        throw runtime_error(
            "SQLite database is not open"
        );
    }

    /*
     * Clear scores left from a previous algorithm or run.
     */
    AlignedGeneEdge::GOWeightCache.clear();

    sqlite3_stmt* statement = NULL;

    const char* sql =
        "SELECT "
        "    g1.gene_symbol, "
        "    g2.gene_symbol, "
        "    s.score "
        "FROM go_pair_score AS s "
        "INNER JOIN gene_uniprot AS g1 "
        "    ON g1.uniprot_accession = s.protein1 "
        "   AND g1.is_preferred = 1 "
        "INNER JOIN gene_uniprot AS g2 "
        "    ON g2.uniprot_accession = s.protein2 "
        "   AND g2.is_preferred = 1 "
        "WHERE s.algorithm = ?1";

    int result =
        sqlite3_prepare_v2(
            Petal::con,
            sql,
            -1,
            &statement,
            NULL
        );

    if (result != SQLITE_OK)
    {
        throw runtime_error(
            string(
                "Could not prepare GO cache query: "
            ) +
            sqlite3_errmsg(Petal::con)
        );
    }

    result =
        sqlite3_bind_text(
            statement,
            1,
            AlignedGeneEdge::GOAlgorithm.c_str(),
            -1,
            SQLITE_TRANSIENT
        );

    if (result != SQLITE_OK)
    {
        sqlite3_finalize(statement);

        throw runtime_error(
            "Could not bind GO algorithm"
        );
    }

    int loadedCount = 0;

    while (
        (result = sqlite3_step(statement))
        == SQLITE_ROW)
    {
        const unsigned char* geneAText =
            sqlite3_column_text(
                statement,
                0
            );

        const unsigned char* geneBText =
            sqlite3_column_text(
                statement,
                1
            );

        if (
            geneAText == NULL
            ||
            geneBText == NULL)
        {
            sqlite3_finalize(statement);

            throw runtime_error(
                "GO cache contains a null gene symbol"
            );
        }

        string geneA =
            reinterpret_cast<const char*>(
                geneAText
            );

        string geneB =
            reinterpret_cast<const char*>(
                geneBText
            );

        float score =
            static_cast<float>(
                sqlite3_column_double(
                    statement,
                    2
                )
            );

        /*
         * Preserve the behavior of your old cache loader:
         * store the score in both directions.
         */
        AlignedGeneEdge::GOWeightCache[
            make_pair(geneA, geneB)
        ] = score;

        AlignedGeneEdge::GOWeightCache[
            make_pair(geneB, geneA)
        ] = score;

        loadedCount++;
    }

    if (result != SQLITE_DONE)
    {
        string errorMessage =
            sqlite3_errmsg(Petal::con);

        sqlite3_finalize(statement);

        throw runtime_error(
            "Error while loading GO weight cache: " +
            errorMessage
        );
    }

    sqlite3_finalize(statement);

    if (loadedCount == 0)
    {
        throw runtime_error(
            "No GO scores were found for algorithm " +
            AlignedGeneEdge::GOAlgorithm
        );
    }

    cout
        << "Loaded "
        << loadedCount
        << " GO pair scores for "
        << AlignedGeneEdge::GOAlgorithm
        << endl;
}


double AlignedGeneEdge::LoadGONormalization()
{
    if (Petal::con == NULL)
    {
        throw runtime_error(
            "SQLite database is not open"
        );
    }

    sqlite3_stmt* statement = NULL;

    const char* sql =
        "SELECT MAX(score) "
        "FROM go_pair_score "
        "WHERE algorithm = ?1";

    int result =
        sqlite3_prepare_v2(
            Petal::con,
            sql,
            -1,
            &statement,
            NULL
        );

    if (result != SQLITE_OK)
    {
        throw runtime_error(
            string(
                "Could not prepare GO normalization query: "
            ) +
            sqlite3_errmsg(Petal::con)
        );
    }

    result =
        sqlite3_bind_text(
            statement,
            1,
            AlignedGeneEdge::GOAlgorithm.c_str(),
            -1,
            SQLITE_TRANSIENT
        );

    if (result != SQLITE_OK)
    {
        sqlite3_finalize(statement);

        throw runtime_error(
            "Could not bind GO algorithm for normalization"
        );
    }

    result =
        sqlite3_step(statement);

    if (
        result != SQLITE_ROW
        ||
        sqlite3_column_type(
            statement,
            0
        ) == SQLITE_NULL)
    {
        sqlite3_finalize(statement);

        throw runtime_error(
            "Could not find a normalization score for " +
            AlignedGeneEdge::GOAlgorithm
        );
    }

    double maximumScore =
        sqlite3_column_double(
            statement,
            0
        );

    sqlite3_finalize(statement);

    if (maximumScore <= 0.0)
    {
        throw runtime_error(
            "GO normalization score must be greater than zero"
        );
    }

    AlignedGeneEdge::GONormalize =
        maximumScore;

    cout
        << "GO normalization for "
        << AlignedGeneEdge::GOAlgorithm
        << ": "
        << AlignedGeneEdge::GONormalize
        << endl;

    return AlignedGeneEdge::GONormalize;
}





void AlignedGeneEdge::PrintTime(){
	
//	end = clock();
	long int a=(time (NULL)) - start_time;
	std::cout <<yellow<< "["<< a/60 <<":"<<a%60<<"] "<<NC;	
	
}

double AlignedGeneEdge::InformationContentOfGenes( Gene* A, Gene* B){
    if (A == NULL || B == NULL){
        throw runtime_error(
            "Null Gene passed to InformationContentOfGenes"
        );
    }

    map< pair<string, string>, float >::const_iterator found =
		 AlignedGeneEdge::GOWeightCache.find( make_pair( A->Name, B->Name ) );
    if (found == AlignedGeneEdge::GOWeightCache.end()){
        throw runtime_error(
            "Missing cached GO score for " + A->Name + " and " +
            B->Name + " using algorithm " + AlignedGeneEdge::GOAlgorithm
        );
    }

    return found->second;
}


