// GURKAN BEBEK @ Copyright
// July 2026
// No use of this code is allowed 

/* on mac os x

SQLITE_PREFIX="$(brew --prefix sqlite)"


Debug build:
SQLITE_PREFIX="$(brew --prefix sqlite)" && clang++ -std=c++11 -O0 -g -Wall -Wextra -Wpedantic -fsanitize=address,undefined -fno-omit-frame-pointer -I"$SQLITE_PREFIX/include" -L"$SQLITE_PREFIX/lib" gene.cpp petal.cpp AlignedGeneEdge.cpp main.cpp -lsqlite3 -o INSPIRE_debug


production build:
SQLITE_PREFIX="$(brew --prefix sqlite)" && clang++ -std=c++11 -O3 -DNDEBUG -Wall -Wextra -Wpedantic -I"$SQLITE_PREFIX/include" -L"$SQLITE_PREFIX/lib" gene.cpp petal.cpp AlignedGeneEdge.cpp main.cpp -lsqlite3 -o INSPIRE

	
*/


#include "petal.h"
#include <limits>
#include "base_header.h"
#include "AlignedGeneEdge.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <time.h>
#include <math.h>
//using boost::format;
using namespace std;

 
//#define antiCONSTANT_k 255
//#define CONSTANT_k 0.006375

#define antiCONSTANT_k 200
#define CONSTANT_k 0.005	 // = 1/200
#define MISSMATCH 0.5

#define INF 999999
#define NONZERO 1212123
#define NOT_ADDED 563463
#define ADDED 11241
#define MAX 250000  /// this is pretty much |V_1|x|V_2| of AGE 
#define ZERO -1
			//cout << (float)(tree_size * 100 * tree_size)/((max_tree-1)* treecost) <<endl;
//			cout << normalizeScor << maximumtreecostetreecost(tree_size, matreecosttreecostx_tree, treecost) <<endl;

 
time_t start_time;
int VERBOSE;
float normalizeScore( int tree_size, float max_tree, float treecost) {
	//return (float)(tree_size * 100 * tree_size)/((max_tree-1)* treecost);
	return (float)(tree_size * 100 * tree_size)/((max_tree -1)* treecost);
}

void printAGE(AlignedGeneEdge* A,AlignedGeneEdge* B){
	cout << cyan <<"\t\t{A1:" << (*A).A1->Name  << "::B1:" << (*A).B1->Name << "}"<< " --- "
	<< "{C1:" << (*B).A1->Name  << "\tD1:" << (*B).B1->Name << "}"<< endl 
	<< "\t\t{A2:" << (*A).A2->Name  << "\tB2:" << (*A).B2->Name << "}"<< " --- "
	<< "{C2:" << (*B).A2->Name  << "\tD2:" << (*B).B2->Name << "}"<< NC<< endl ;
}

void printAGE(AlignedGeneEdge* A ){
	cout << cyan <<"\t\t[A1. " << (*A).A1->Name  << "::A2. " << (*A).A2->Name << "]"<< " ---  "
		<< "[B1. " << (*A).B1->Name  << "::B2. " << (*A).B2->Name << "}"<< NC<< endl;
 }



int FindAGE(Petal* Petal1, Petal* Petal2, vector<AlignedGeneEdge*> * AGE){
	
	// for each Petal1 gene:
	int counterPetal1=0;
	int counterPetal2=0;
	vector<AlignedGeneEdge *>::iterator itAGE;
	cout<<GREEN<<"AGE for "<<cyan<<Petal1->Name<<" - "<< Petal2->Name <<NC<< endl;
	cout<<purple<<"initial size for AGE: "<<(*AGE).size()<<NC<<endl;
	
	// get edge from Petal 1.
	for (int i=0;i < Petal1->size; i++){
		for (int j=0;j < Petal1->size; j++){
			//should be j=i if edge direction doesnt matter 
			//should be j=0 if edge direction matters 
			if( Petal1->matrix[i][j]!=0){
				counterPetal1++;
				// get edge from Petal 2.
				for (int i2=0;i2 < Petal2->size; i2++){
					for (int j2=0;j2 < Petal2->size; j2++){
						if( Petal2->matrix[i2][j2]!=0){
							
							AlignedGeneEdge* pAGE = new AlignedGeneEdge(Petal1->GenesList[i],Petal1->GenesList[j],
																		Petal2->GenesList[i2],Petal2->GenesList[j2]);
							//cout <<"test" << Petal1->GenesList[i]->Name << "," << Petal2->GenesList[i2]->Name << endl;
							//cout << "test" << Petal1->GenesList[j]->Name << "," << Petal2->GenesList[j2]->Name << endl;
																		//Gene * A1, Gene * A2, Gene * B1, Gene * B2);
							++counterPetal2;
							pAGE->SetWeight();

//yaz						cout<< RED<<counterPetal2<<" "<<green<<pAGE.Weight<<purple<<std::endl;
							
							(*AGE).push_back( pAGE );
						}
					}
				}
			}					
		}
	}
	cout<< "counterPetal1 "<<counterPetal1<<endl;
	cout<< "counterPetal2 "<<counterPetal2/counterPetal1<<endl;
	cout<< "counterPetal1*counterPetal2 "<<counterPetal2<<endl;
	
	cout<<"(*AGE).size()="<<(*AGE).size()<<endl;
	
//	ptrPetal->GenesList[index]
	
	
	
	cout<<NC<<"---- "<<endl;
	
	return 0;
}
void writeAGEgraph(char const  *s,unsigned char*** A,int size){
	int j=0;
	int k=0;
	int count=0;
	ofstream myfile;
	myfile.open( s);
	myfile<<"# " <<size<<endl;
	for( j=0;j<size;j++){
		for(k=0;k<size;k++){
			if(( (*A)[j][k]) >0){
				count++;
				if(j>k)
					myfile <<j<<"\t"<<k<<"\t" << ((*A)[j][k]) <<endl;
				else
					myfile <<j<<"\t"<<k<<"\t" << ((*A)[k][j]) <<endl;
			}
		}
	}
	cout<< "# non zero entries:" << count<<endl;
	myfile.close();
}

void writeAGEgraph3(char const  *s,unsigned short*** A,int size){
	int j=0;
	int k=0;
	int count=0;
	ofstream myfile;
	myfile.open( s);
	myfile<<"# " <<size<<endl;
	for( j=0;j<size;j++){
		for(k=0;k<size;k++){
			if(( (*A)[j][k]) >0){
				count++;
				myfile <<j<<"\t"<<k<<"\t" //<< scientific
				 << ((*A)[j][k]) 
				<<endl;
			}
		}
	}
	cout<< "# non zero entries:" << count<<endl;
	myfile.close();
}

void writeAGEgraph2(char const  *s,float*** A,int size){
	int j=0;
	int k=0;
	int count=0;
	ofstream myfile;
	myfile.open( s);
	myfile<<"# " <<size<<endl;
	for( j=0;j<size;j++){
		for(k=0;k<size;k++){
		//	if(( (*A)[j][k]) <0){
				count++;
				myfile <<j<<"\t"<<k<<"\t" //<< scientific
				 << ((*A)[j][k]) 
				<<endl;
		//	}
		}
	}
	cout<< "# non zero entries:" << count<<endl;
	myfile.close();
}



 


//n is size of prim's tree
//src_node is where we start from
int prims_spanningtree(unsigned int src_node, unsigned char*** AGEgraph, int*** tree,unsigned int n, vector<AlignedGeneEdge*>*  AGE)
{

	float* d;
	int*   visited;
	int*   parent;

	int k, stcost;
	int minumum;
	unsigned int u=0;
	unsigned int p=0;
	unsigned int q=0;
	unsigned int v=0;
//	float** cost;

//*	cout<< "Allocs..."<<endl;
	d = (float*) malloc(MAX * sizeof(float)) ; 
	if (d == NULL) { cout<<"COULDNT ALLOC MEM for d"<<endl;exit( -1);}
	visited = (int*) malloc(MAX * sizeof(int));
	if (visited == NULL) { cout<<"COULDNT ALLOC MEM for visited"<<endl;exit( -1);}
	parent  = (int*) malloc(MAX * sizeof(int));
	if (parent == NULL) { cout<<"COULDNT ALLOC MEM for parent"<<endl;exit( -1);}

//	cost = (float**) malloc( n * sizeof(float *));
//	if (cost == NULL) { cout<<"COULDNT ALLOC MEM for cost"<<endl;exit( -1);}

	int a;
	//duplicate cost array
	/* for each pointer, allocate storage for an array of ints */
//	for (p = 0; p < n; p++){
//		cost[p] = (float*) malloc(n * sizeof(float));
//		if ((cost[p]) == NULL) { cout<<"COULDNT ALLOC MEM for cost[p]"<<endl;exit( -1);}
		
		//		for (q = 0; q < n; q++){
//			if(p>=q)
//				a=((*AGEgraph)[p][q]);	
//			else
//				a=((*AGEgraph)[q][p]);	

		//	cout<< YELLOW<< p << " "<< q <<" " <<a<<endl;
//			cost[p][q]=0;
//			if (a!=0)
//				cost[p][q]=(float)25500/a;
//			else
//				cost[p][q]=INF;
//			cout<< GREEN << p << " "<< q <<" " << cost[p][q] <<endl;
//		}
//	}
//	writeAGEgraph2("newAdj.txt",&cost,n);

//*	cout<< green << "Start prim...."<<NC<<endl;

//*	cout<< "source node: ";
//*	printAGE( ((*AGE).at(src_node) ));
	
	for (p = 0; p < n; p++)
	{
		parent[p] = src_node;

//		a=cost[src_node][p];
		// this line is replaced with
		if(src_node>=p)
			a=((*AGEgraph)[src_node][p]);	
		else
			a=((*AGEgraph)[p][src_node]);	
		if(a!=0)
			a=antiCONSTANT_k*100/a;
		else
			a=INF;

 		d[p] = a;
		if ( a == 0 ) {
		cout<<((*AGE).at(p))->A1->Name<<endl;
		cout<<((*AGE).at(p))->A2->Name<<endl;
		cout<<((*AGE).at(p))->B1->Name<<endl;
		cout<<((*AGE).at(p))->B2->Name<<endl;
			exit(1);

		}
		visited[p] = 0;
//*		cout<< purple<<  p  <<" "<< d[p] << " " <<endl;
//*		if(d[p] != INF){
//*			printAGE( ((*AGE).at(p) ));
//*		}
		
	}
	u = src_node;
//*	int count=0;
	
	for (v = 0; v < n; v++){
		if ( (((*AGE).at(v))->A1 == ((*AGE).at(u))->A1 &&
			 ((*AGE).at(v))->A2 == ((*AGE).at(u))->A2 )||
			 (((*AGE).at(v))->B1 == ((*AGE).at(u))->B1 &&
			 ((*AGE).at(v))->B2 == ((*AGE).at(u))->B2 ) ||

			 (((*AGE).at(v))->A1 == ((*AGE).at(u))->A2 &&
				 ((*AGE).at(v))->A2 == ((*AGE).at(u))->A1 )||
				 (((*AGE).at(v))->B1 == ((*AGE).at(u))->B2 &&
				 ((*AGE).at(v))->B2 == ((*AGE).at(u))->B1 )
			){
//*				cout<<count<<RED  <<"v "<<v<< " ++++"; printAGE( ((*AGE).at(v)) );
//				cout<<count++<<RED<<"u "<<u<< " ++++"; printAGE( ((*AGE).at(u)) );
  			visited[v]=1;


		}
	}

//I think this is the actual traversal
	visited[src_node] = 1;
	stcost = 0;
	k = 0;
	for (p = 0; p < (n-1); p++){

//*		cout<< yellow<<"u:"<<u <<"\t parent[p]:"<< parent[p] << ""<<endl;
		
		minumum = INF; // large #
	
		for (q = 0; q < n; q++){
			if (!visited[q] &&   (d[q] < minumum)){
				minumum = d[q];
				u = q;
			}
//*			else
//				cout <<cyan<< q << " "<< visited[q] <<" "<< d[q] << " "<< endl;
		}
		if(minumum==INF){
//*			cout<< " min=INF "<<endl;
			break;
		}
//*		cout<< red <<"u:"<<u <<"\t parent[p]:"<< parent[p] << ""<<endl;
		
		visited[u] = 1;
		
		// also visit the whole column and row 
		// GCL commenting tihs out
	  	for (unsigned int v = 0; v < n; v++){
	  		if (

/* This is now redundant, because of the following for loop
 * (((*AGE).at(v))->A1 == ((*AGE).at(u))->A1 &&
				 ((*AGE).at(v))->A2 == ((*AGE).at(u))->A2 )||
				 (((*AGE).at(v))->B1 == ((*AGE).at(u))->B1 &&
				 ((*AGE).at(v))->B2 == ((*AGE).at(u))->B2 ) ||
*/

/*
				 (((*AGE).at(v))->A1 == ((*AGE).at(u))->A2 &&
					 ((*AGE).at(v))->A2 == ((*AGE).at(u))->A1 )||
					 (((*AGE).at(v))->B1 == ((*AGE).at(u))->B2 &&
					 ((*AGE).at(v))->B2 == ((*AGE).at(u))->B1 )
				){
*/
//This deletes the other direction
				 (((*AGE).at(v))->A1 == ((*AGE).at(u))->A2 &&
					 ((*AGE).at(v))->A2 == ((*AGE).at(u))->A1 )&&
					 (((*AGE).at(v))->B1 == ((*AGE).at(u))->B2 &&
					 ((*AGE).at(v))->B2 == ((*AGE).at(u))->B1 )
				){
//* 					cout<<count<<RED<<"v "<<v<< " ++++"; printAGE( ((*AGE).at(v)) );
//*					cout<<count++<<RED<<"u "<<u<< " ++++"; printAGE( ((*AGE).at(u)) );
					visited[v]=1;
	  		}

		}		
		//"visit" any nodes that do not agree with our most recent alignment
	  	for (unsigned int v = 0; v < n; v++){
			//So, if A1 =A1, then B1 HAS TO BE equal to B1, otherwise we have a conflicting alignment
	  		if ( 
				(((*AGE).at(v))->A1 == ((*AGE).at(u))->A1 &&
				(((*AGE).at(v))->B1 != ((*AGE).at(u))->B1)) ||
				//If A1=A2, then  B1=B2 must be true
	  		      (((*AGE).at(v))->A1 == ((*AGE).at(u))->A2 &&
				(((*AGE).at(v))->B1 != ((*AGE).at(u))->B2)) ||
				
				//If B1=B1 then A1=A1
	  		        (((*AGE).at(v))->B1 == ((*AGE).at(u))->B1 &&
				(((*AGE).at(v))->A1 != ((*AGE).at(u))->A1)) ||

				//If B1=B2 then A1=A2
				 (((*AGE).at(v))->B1 == ((*AGE).at(u))->B2 &&
				 (((*AGE).at(v))->A1 != ((*AGE).at(u))->A2  )) ||


				//And now for the other side
				(((*AGE).at(v))->A2 == ((*AGE).at(u))->A2 &&
				(((*AGE).at(v))->B2 != ((*AGE).at(u))->B2)) ||
				//If A2=A1, then  B2=B1 must be true
	  		      (((*AGE).at(v))->A2 == ((*AGE).at(u))->A1 &&
				(((*AGE).at(v))->B2 != ((*AGE).at(u))->B1)) ||
				
				//If B2=B2 then A2=A2
	  		        (((*AGE).at(v))->B2 == ((*AGE).at(u))->B2 &&
				(((*AGE).at(v))->A2 != ((*AGE).at(u))->A2)) ||

				//If B2=B1 then A2=A1
				 (((*AGE).at(v))->B2 == ((*AGE).at(u))->B1 &&
				 (((*AGE).at(v))->A2 != ((*AGE).at(u))->A1  ))
){
// 					cout<<count<<RED<<"v "<<v<< " ++++"; printAGE( ((*AGE).at(v)) );
//					cout<<count++<<RED<<"u "<<u<< " ++++"; printAGE( ((*AGE).at(u)) );
					visited[v]=1;
	  		}

		}		

		
		stcost = stcost + d[u];
		(*tree)[k][1] = parent[u];
		(*tree)[k][2] = u;
		(*tree)[k][3] = d[u];
		(*tree)[k++][4] = ADDED;
		
//*		cout<<NC<< k  << ". added " << parent[u]<< " to " << u <<" with cost " << d[u] <<endl;
		for (unsigned int v = 0; v < n; v++)
		{
			int cost_u_v;
			if(u>=v)
				cost_u_v=((*AGEgraph)[u][v]);	
			else
				cost_u_v=((*AGEgraph)[v][u]);	
			if(cost_u_v!=0)
				cost_u_v=antiCONSTANT_k*100/cost_u_v;
			else
				cost_u_v=INF;

//			if (!visited[v] && d[v]>0 && (cost[u][v] < d[v]) ) // ignore 0s in the matrix
			if (!visited[v] && d[v]>0 && (cost_u_v < d[v]) ) // ignore 0s in the matrix
			{
				d[v] = cost_u_v;
				parent[v] = u;
			}
		}	
    }

	free(d);
	free(visited);
	free(parent);

//	for (p = 0; p < n; p++)
//		free (cost[p]);
//	free(cost);	
	if (k==0)
		return numeric_limits<int>::max();
	return (stcost);
}


int prims_display(float stcost,int*** tree,int n,vector<AlignedGeneEdge*>* AGE, bool print)
{
	unsigned int i,j=0;
	
	for (i = 0; i <(unsigned int)(n-1); i++){
		if(( (*tree)[i][4]) == ADDED){
			if ( print) {
				cout << (*tree)[i][1] << " " << (*tree)[i][2] <<" cost: " << (*tree)[i][3] << endl;
				printAGE( (*AGE).at((*tree)[i][1]));
				printAGE( (*AGE).at((*tree)[i][2]));
			}
			j++;
 		}
    }
	if ( print)
		cout
		    << "\nDisplayed tree score: "
		    << stcost
		    << endl
		    << endl;
	return j;
}


int main(int argc, char *argv[])
{
	
	if (argc < 4)
	{
	    cout
	        << "Usage: "
	        << argv[0]
	        << " DATABASE PETAL_A PETAL_B [GO_ALGORITHM]"
	        << endl;

	    return -1;
	}

	if (
	    Petal::OpenDatabase(
	        argv[1]
	    )
	    != 0)
	{
	    return -1;
	}

	string goAlgorithm =
	    "inspire_legacy_v1";

	if (argc >= 5)
	{
	    goAlgorithm =
	        argv[4];
	}

	AlignedGeneEdge::SetGOAlgorithm(goAlgorithm);

	AlignedGeneEdge::LoadGONormalization();
	//AlignedGeneEdge::GONormalize = 23.7367; 
	cout
	    << "GO normalization in use: "
	    << AlignedGeneEdge::GONormalize
	    << endl;
	
	AlignedGeneEdge::LoadGOWeightCache();
	
	
	
  	start_time = time (NULL);
	time_t end_time;
	VERBOSE=0;//0: no verbose - >0: yes verbose

 	
	int ape_count;
	int treecost=0;
	//int edge_type_count[3];
 
	cout<< "Create two petals"<<endl;
	// Create two petals
	vector<Petal*> Petals  = vector<Petal*>(2);
	// Load petals 
	// TODO this should be one at a time with limited parameters. (hence, 'ape_count' should be the cartesian product.)
	ape_count = Petal::LoadPetals(argv[2],argv[3],&Petals);
	cout<< "petals are loaded"<<endl;
 	
	// empty network...
	if(ape_count==-1){
		cout<< "+++ One of the networks has no interaction +++ NO CIGAR +++"<<endl;
		cout << "$RESULT$\t" << (Petals[0])->Name << "\t"<< (Petals[1])->Name <<  "\tnorm_score:\t" << "-0\t\tNO_RESULT"<<endl;
		return 0;
		
	}
	
	
	
	// Create APE or Aligned Gene Edges
	// AGE is a vector of size 0 initially
	vector<AlignedGeneEdge*>  AGE(0);


	
	// FindAGE Fills out AGE with Petals 0 and 1.
	FindAGE((Petals[0]), (Petals[1]), &AGE);
	
	int age_size = AGE.size();
	//allocate space for modified prim start postion.
	int* startNode=(int*) malloc( age_size*sizeof(int));

	if (startNode == NULL) { cout<<"COULDNT ALLOC MEM for startNode"<<endl;exit( -1);}

	for (int i=0;i<age_size;i++)
		startNode[i]=0;
		
	// construct AGE graph  
	
	// create array
	cout<<YELLOW<< "AGE nodes: "<<age_size<<NC<<endl;

	/*  allocate storage for an array of pointers */
	cout<< "Allocate AGEgraph :  "<< age_size<<" x "<<age_size<<"\t"<< age_size*sizeof(unsigned char)*age_size<<" bytes"<<endl;

	unsigned char** AGEgraph = (unsigned char**) malloc( (age_size) * sizeof(unsigned char *));
	if (AGEgraph == NULL) { cout<<"COULDNT ALLOC MEM for AGEgraph"<<endl;exit( -1);}

	/* for each pointer, allocate storage for an array of ints */
	for (  int i = 0; i < age_size; i++) {
		AGEgraph[i] = (unsigned char*) malloc(i+1 * sizeof(unsigned char));
		if ((AGEgraph[i]) == NULL) { cout<<"COULDNT ALLOC MEM for AGEgraph[i]"<<endl;exit( -1);}
		for(  int j=0; j<i;j++){
			AGEgraph[i][j]=0;	
		}
	}
	//spanning tree allocation
	int** tree = (int**) malloc( (age_size) * sizeof(int *));
	if (tree == NULL) { cout<<"COULDNT ALLOC MEM for tree"<<endl;exit( -1);}

	/* for each pointer, allocate storage for an array of ints */
	for (  int i = 0; i < age_size; i++) {
		tree[i] = (int*) malloc(5 * sizeof(int));
		if ((tree[i]) == NULL) { cout<<"COULDNT ALLOC MEM for tree[i]"<<endl;exit( -1);}
		
		for(unsigned int j=0; j<5;j++){
			tree[i][j]=0;	
			
		}
		tree[i][4]=NOT_ADDED;
	}

	// add links
	// # check <A> edge alignemnt / <B> Node duplication

	vector<AlignedGeneEdge*>::iterator it;		
	vector<AlignedGeneEdge*>::iterator it2;		
	int i=0;
	for(it=AGE.begin(); it!=AGE.end();++it){
		if((*it)->adjIndex==-1){
			(*it)->adjIndex=i;
		}
		else{
			cout<<"SCREW"<<endl;
			return 0;
		}
		i++;
	}

	end_time  = time (NULL);
 	printf("\n[%ld min %ld sec] Nodes are done in \n",(end_time-start_time)/60,(end_time-start_time)%60 );
	
	i=0;
	int j=0;
	// for each APE
	for(it=AGE.begin(); it!=AGE.end();++it){

		// ---------------------- 		# <A>
		// CASE A: edge extension 		# A1-[A2/a1]-a2 (d1-d2-d3)
		// ---------------------- 		# B1-[B2/b1]-b2 (d1-d2-d3)
		// for each A1~B1 and A2~B2 
		j=0;
		// find A3~B3 such that 

		for(it2=AGE.begin(); it2!=AGE.end();it2++){
			// cout<<i<<"\t"<<*it<<"\t"<<j<<"\t"<<*it2<<"\t"<<AGEgraph[i][j]<<endl;
			//first: A3~B3 is same as A2~B2
			//second: A3 != A1 and B3 != B1

			if(i!=j){	
				if (  (((*it)->A2) == ((*it2)->A1) 
				    && ((*it)->B2) == ((*it2)->B1) 
				    && ((*it)->A1) != ((*it2)->A2) 
				    && ((*it)->B1) != ((*it2)->B2) 
					&& (*it)!=(*it2)) 

/*					||   (((*it)->A1) == ((*it2)->A2) 
					    && ((*it)->B1) == ((*it2)->B2) 
					    && ((*it)->A2) != ((*it2)->A1) 
					    && ((*it)->B2) != ((*it2)->B1) 
						&& (*it)!=(*it2))
*/								
					){
        	
					if( ((*it)->Weight)==0 && ((*it2)->Weight)==0 ){
						//AGEgraph[i][j]=0;
						//AGEgraph[j][i]=0;
					}
					else{
						int ceiled=ceil((((*it)->Weight)*((*it2)->Weight))*CONSTANT_k);
						
						#ifdef DEBUG
						if(ceiled>=255){
							cout<<"ERROR 4 - unsigned char maxed out "<<endl;return 4;
						}
						#endif /* DEBUG */
						
						
						//if (AGEgraph[i][j] >0){
						//	cout<<RED<<"ERROR 5 - 1 "<<endl;
						//}	
 				//	//	if(ceiled > AGEgraph[i][j])
						if(i>j){
							AGEgraph[i][j]=(unsigned char) ceiled;
						//	 cout<< red  <<i<<" "<<j<<" "<<yellow<< "AGEgraph["<<i<<"]["<<j<<"]: " <<( int) AGEgraph[i][j]<<"\t"<< "((*it)->Weight): "<<((*it)->Weight) <<"\t"<<"((*it2)->Weight): "<< ((*it2)->Weight) << "\tceiled: "<<  ceiled << endl;	
	
						}
						else{
							AGEgraph[j][i]=(unsigned char) ceiled;
						//	 cout<< red  <<j<<" "<<i<<" "<<yellow<< "AGEgraph["<<j<<"]["<<i<<"]: " <<( int) AGEgraph[j][i]<<"\t"<< "((*it)->Weight): "<<((*it)->Weight) <<"\t"<<"((*it2)->Weight): "<< ((*it2)->Weight) << "\tceiled: "<<  ceiled << endl;	
	
	
						}				
					
							
					 
						startNode[i]=1;
						startNode[j]=1;
						
					//	cout<< green<<i<<" "<<j<<" "<<yellow<<  scientific<< AGEgraph[i][j]  << endl;
					//	cout<< green<<i<<" "<<j<<" "<<yellow<<  scientific<< AGEgraph[i][j]  << endl;
						#ifdef DEBUG
						if (i!= (*it)->adjIndex ) {cout<<"ERROR 3"<<endl;return 3;}
						if (j!= (*it2)->adjIndex) {cout<<"ERROR 3"<<endl;return 3;}
						#endif /* DEBUG */
//						printAGE( (AGE.at(i)) ); // cout<< ((*it)->A1)->Name << ((*it)->A2)->Name;
//						printAGE( (AGE.at(j)) ); // cout<< ((*it2)->A1)->Name << ((*it2)->A2)->Name;
						
					//	cout<< "writing AGE Graph to ["<<i<<"]["<<j<<"] "<< AGEgraph[i][j] << "\t AGE.at(i,j)->adjIndex shoud be ="<< AGE.at(i)->adjIndex <<" "<<AGE.at(j)->adjIndex <<endl;
        	
					}
				}
//			}//endif(i!=j)
			// cout<<i<<"\t"<<*it<<"\t"<<j<<"\t"<<*it2<<"\t"<<AGEgraph[i][j]<<endl;
			
			// cout<<i<<"\t"<<j<<"\t"<<AGEgraph[i][j]<<endl;
			
			
//			if(i>j){	
				else if ( (  ((*it)->A2) == ((*it2)->A1) 
				    && ((*it)->A1) != ((*it2)->A2) 
					&& (*it)!=(*it2) 
					&& ((*it)->B1) != ((*it2)->B1) 
					&& ((*it)->B1) != ((*it2)->B2) 
					&& ((*it)->B1) != ((*it2)->B2)	
					&& ((*it)->B2) != ((*it2)->B2) 
					&& (Petals[1])->checkInteract( ((*it)->B2)->Name, ((*it2)->B1)->Name) 
					) || 
					(  ((*it)->B2) == ((*it2)->B1) 
					    && ((*it)->B1) != ((*it2)->B2) 
						&& (*it)!=(*it2) 
						&& ((*it)->A1) != ((*it2)->A1) 
						&& ((*it)->A1) != ((*it2)->A2) 
						&& ((*it)->A1) != ((*it2)->A2)	
						&& ((*it)->A2) != ((*it2)->A2) 
						&& (Petals[0])->checkInteract( ((*it)->A2)->Name, ((*it2)->A1)->Name) 
						)
						
					) {
				
					if( ((*it)->Weight)==0 && ((*it2)->Weight)==0 ){
						//AGEgraph[i][j]=0;
						//AGEgraph[j][i]=0;
					//	cout<< red  <<i<<" "<<j<<" "<<yellow<<  scientific<< AGEgraph[i][j]<<endl;					

					}
					else{
						
						int ceiled=ceil((((*it)->Weight)*((*it2)->Weight))*MISSMATCH*CONSTANT_k);
						#ifdef DEBUG
						if( ceiled >= 255){
							cout<<"ERROR 4 - unsigned char maxed out "<<endl;return 4;
						}	
						#endif /* DEBUG */
						
						
						//if (AGEgraph[i][j] >0){
						//	cout<<RED<<"ERROR 5 - 2 "<<endl;
						//}
						
				//		if (ceiled>AGEgraph[i][j])
						
						if(i>j){
							AGEgraph[i][j]=(unsigned char) ceiled;
//							 cout<< red  <<i<<" "<<j<<" "<<yellow<< "AGEgraph["<<i<<"]["<<j<<"]: " <<( int) AGEgraph[i][j]<<"\t"<< "((*it)->Weight): "<<((*it)->Weight) <<"\t"<<"((*it2)->Weight): "<< ((*it2)->Weight) << "\tceiled: "<<  ceiled << endl;	
	
						}
						else{
							AGEgraph[j][i]=(unsigned char) ceiled;
//							 cout<< red  <<j<<" "<<i<<" "<<yellow<< "AGEgraph["<<j<<"]["<<i<<"]: " <<( int) AGEgraph[j][i]<<"\t"<< "((*it)->Weight): "<<((*it)->Weight) <<"\t"<<"((*it2)->Weight): "<< ((*it2)->Weight) << "\tceiled: "<<  ceiled << endl;	
						}
						

						#ifdef DEBUG
						if (i!= (*it)->adjIndex ) {cout<<"ERROR 3"<<endl;return 3;}
						if (j!= (*it2)->adjIndex) {cout<<"ERROR 3"<<endl;return 3;}
						#endif /* DEBUG */
					}
				}// else if
   		//	}//endif(i!=j)
//
		//	if(i>j){	
				else if ( ( ((*it)->A1) != ((*it2)->A1) 
					&& ((*it)->A1) != ((*it2)->A2) 
					&& ((*it)->A1) != ((*it2)->A2)	
					&& ((*it)->A2) != ((*it2)->A2) 
					&& (Petals[0])->checkInteract( ((*it)->A2)->Name, ((*it2)->A1)->Name) 
					&& (*it)!=(*it2) 
					&& ((*it)->B1) != ((*it2)->B1) 
					&& ((*it)->B1) != ((*it2)->B2) 
					&& ((*it)->B1) != ((*it2)->B2)	
					&& ((*it)->B2) != ((*it2)->B2) 
					&& (Petals[1])->checkInteract( ((*it)->B2)->Name, ((*it2)->B1)->Name) )	) {
					if( ((*it)->Weight)==0 && ((*it2)->Weight)==0 ){
						//AGEgraph[i][j]=0;
						//AGEgraph[j][i]=0;
// 						cout<< red  <<i<<" "<<j<<" "<<yellow<<  scientific<< AGEgraph[i][j]<<endl;					
					}
					else{
						
						int ceiled=ceil((((*it)->Weight)*((*it2)->Weight))*MISSMATCH*CONSTANT_k);
						#ifdef DEBUG
						if( ceiled >= 255){
							cout<<"ERROR 4 - unsigned char maxed out "<<endl;return 4;
						}	
						#endif /* DEBUG */
					//	if (AGEgraph[i][j]>0){
					//		cout<<RED<<"ERROR 5 - 3"<<endl;
					//	}
					   if(i>j){
					   	AGEgraph[i][j]=(unsigned char) ceiled;
					   //	 cout<< red  <<i<<" "<<j<<" "<<yellow<< "AGEgraph["<<i<<"]["<<j<<"]: " <<( int) AGEgraph[i][j]<<"\t"<< "((*it)->Weight): "<<((*it)->Weight) <<"\t"<<"((*it2)->Weight): "<< ((*it2)->Weight) << "\tceiled: "<<  ceiled << endl;	
                       
					   }
					   else{
					   	AGEgraph[j][i]=(unsigned char) ceiled;
					   	// cout<< red  <<j<<" "<<i<<" "<<yellow<< "AGEgraph["<<j<<"]["<<i<<"]: " <<( int) AGEgraph[j][i]<<"\t"<< "((*it)->Weight): "<<((*it)->Weight) <<"\t"<<"((*it2)->Weight): "<< ((*it2)->Weight) << "\tceiled: "<<  ceiled << endl;	
					   }
					
						#ifdef DEBUG
						if (i!= (*it)->adjIndex ) {cout<<"ERROR 3"<<endl;return 3;}
						if (j!= (*it2)->adjIndex) {cout<<"ERROR 3"<<endl;return 3;}
						#endif /* DEBUG */

					}
				}
	//			end_time  = time (NULL);
	//			printf("\n[%ld min %ld sec] j=%d \n",(end_time-start_time)/60,(end_time-start_time)%60 ,j);
			

			}//endif(i!=j)


	
			j++;
	
		}
//		cout<< "writing AGE Graph to [0][0] "<< AGEgraph[0][0] << endl;
//		cout<< "writing AGE Graph to [" << age_size-1<< "][" <<age_size-1<< "]"<< AGEgraph[age_size-1][age_size-1] <<endl;
			
//	/		end_time  = time (NULL);
//		 	printf("\n[%ld min %ld sec]  i=%d \n",(end_time-start_time)/60,(end_time-start_time)%60 ,i);
		
		i++;
	}
	end_time  = time (NULL);
 	printf("\n[%ld min %ld sec] Edges are done in \n",(end_time-start_time)/60,(end_time-start_time)%60 );

//	char const *s1 = "AGEadj2.txt";

//	 " <<s1<< AGEgraph[0][0]<< "\t AGE.at(0)->adjIndex"<< AGE.at(0)->adjIndex ;
//	cout<< "writing AGE Graph to " << age_size-1<< " " 
//								   << AGEgraph[age_size-1][age_size-1]
//								   << "\t AGE.at(age_size-1)->adjIndex"<< AGE.at(age_size-1)->adjIndex ;
//	printAGE( (AGE.at(0)) );

/*
	// USE THIS TO GET ADJ MATRIX OUTPUT
	cout<< "writing AGE Graph to" <<s1<<endl;
	writeAGEgraph(s1,&AGEgraph,age_size);
*/

// to plot use 
// on the laptop
/// splot2.pl  -style points  -plotcmds /Users/gurkan/Research/blast/blast/drosophila_data/src/heatmap.gnuplot  -using "1:2:3"  -showgnuplot -outfile AGEadj2.png AGEadj2.txt

	  
	// ClusterSearchPhase($cluster_out,$s1_ppi, $s2_ppi,$APE_out, $APEmap_out,$overlap_rate, $k,$node_score, $score_limit);
	
	// allocate space for spanning tree
	
	// search using   prims_spanningtree(int src_node, int*** cost, int*** tree,int n)
		
	/// MULTIPLE VERSION
	//int	minumum=numeric_limits<int>::max();
	float	maximum=0;
	int start_node=-1;
	j=0;
	for (i =0;i<(age_size);i++){

		//// addition to here 
		//if not exact node match count is 0 do the following -
		///end if
		///else just use the exact matches to seed faster 
		if(startNode[i]==1){

// 			cout<< j << "th RUN"<<endl;
			int max_tree=min((Petals[0])->interactionSize, (Petals[1])->interactionSize);
			treecost=prims_spanningtree(i, &AGEgraph, &tree, (unsigned int)(age_size), &AGE );
			if ( treecost == 0 ) {
				exit(1);
			}
			int tree_size=prims_display(treecost, &tree, age_size, &AGE, false );
			float normalizedTreeCost =(float) normalizeScore(tree_size, max_tree, treecost);
			// 1. Initialize d[s] with 0, P[s] with 0, and
			//if (treecost< minumum){
			if (normalizedTreeCost > maximum){
				//minumum=normalizedTreeCost;
				maximum=normalizedTreeCost;
				start_node=i;
				
 				cout<<i<< "th run results ("<< j++ <<") [new minumum]:" << maximum<<endl;
				prims_display(normalizedTreeCost, &tree, age_size, &AGE, true );	
				
			}	
			// only display if it is worth something
//			if(treecost < (numeric_limits<int>::max())){

			// THIS is even less visualization:


			//reset tree tracking.
			for ( int k = 0; k < age_size; k++) {
			 
				tree[k][4]=NOT_ADDED;
			}
//			cout<<"."<<endl;
		}
		///end else
	}
//	cout<< "minumum: "<< treecost<< " "<< start_node<<endl;
	if(start_node==-1){
		cout<< "+++ None of the starting nodes returned significant results +++ NO CIGAR +++"<<endl;
		cout << "$RESULT$\t" << (Petals[0])->Name << "\t"<< (Petals[1])->Name <<  "\tnorm_score:\t" << "-0\t\tNO_RESULT"<<endl;
		//return 0;
	}
	else{
		treecost=prims_spanningtree(start_node, &AGEgraph, &tree, (unsigned int)(age_size), &AGE );
	
		//SINGLE VERSION	treecost=prims_spanningtree(0, &AGEgraph, &tree, (unsigned int)(age_size), &AGE );
	
		int tree_size=prims_display(treecost, &tree, age_size, &AGE, false );
		
		
		//	normalized score
		cout <<"Petal1:"<<(Petals[0])->interactionSize<<"\tPetal2:"<<(Petals[1])->interactionSize<<"\ttreecost: "<<treecost<<"\t tree size:"<<tree_size <<endl;
		int max_tree;
		if(treecost !=0){
			// TAKE THE UNION instead of min.
			max_tree=min((Petals[0])->interactionSize, (Petals[1])->interactionSize);
			// OLD
			// cout << "$RESULT$\t" << (Petals[0])->Name << "\t"<< (Petals[1])->Name <<  "\tnorm_score:\t" << (float)(( max_tree - 1)*100 )*((float)tree_size/max_tree)/treecost <<endl;

			// NEW  (size/max)(100*size/score)
			cout << "$RESULT$\t" << (Petals[0])->Name << "\t"<< (Petals[1])->Name <<  "\tnorm_score:\t";
			//cout << (float)(tree_size * 100 * tree_size)/((max_tree-1)* treecost) <<endl;
			cout << normalizeScore(tree_size, max_tree, treecost) <<endl;
		}
		else{
			cout << "$RESULT$\t" << (Petals[0])->Name << "\t"<< (Petals[1])->Name <<  "\tnorm_score:\t" << "-0\t\tNO_RESULT"<<endl;
		}
		cout<<NC;
	}

	/* now for each pointer, free its array of ints */
	for ( int i = 0; i < age_size ; i++) {
		free(AGEgraph[i]);
		free(tree[i]);
	}
	free(tree);
	free(AGEgraph);
	free(startNode);

	end_time  = time (NULL);
 	printf("\n[%ld min %ld sec] All done in \n",(end_time-start_time)/60,(end_time-start_time)%60 );
	return 0;
}



