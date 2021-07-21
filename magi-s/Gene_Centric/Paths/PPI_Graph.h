#ifndef _PPI_GRAPH
#define _PPI_GRAPH

#include <stdio.h>
#include <stdlib.h>

const int geneNameLen=50;
const int maxNumNode=50000;
extern int numNodes;
extern int coExpresGeneNum;
extern int totalSevereMutInCases;
extern int totalMissenseMutInCases;
extern int totalSevereMutInControl;
extern int totalLengthGenes;
//extern float log_sumNumNodes;
extern float coexpressionPValue[102];// [0, 0.01, 0.02, ..., 0.99, 1]: the co-expression p-value for every edge (for very coExpresionValue > x 
extern float meanCoExpression;
extern float varianceCoExpression;
typedef struct PPI_Node{
	int nodeId;
	char nodeName[geneNameLen]; // protein/gene name

	int numSevereMutInCases;
	int numMissenseMutInCases;
	int numSevereMutInControl;
	
	double prob;

	double weightCases;// The weight assigned for the cases I assume is double 
	double potentialWeightCases; // TODO if the user chooses max for multiple seeds, this is where to store the weight for a gene when you're checking if it's the maximum score or not 
	int weightControl;//The weight assigned to control I assume is integer
	int degree; //number of connections in the PPI network
	double log_length;
	int length;
	int *neighbours; // id of list of nodes each node is connected to
}PPI_Node;

typedef struct coExpresionGeneHash{
	int nodeId;// The gene id in the PPI network table
	char geneName[geneNameLen];
	int hashId; // The coexpression table Id
}coExpresionGeneHash;

extern coExpresionGeneHash coExpresionGeneHashTable[maxNumNode];

extern PPI_Node listNodes[maxNumNode];
extern float coExpresionMatrix[maxNumNode][maxNumNode];

int createPPI_Graph(FILE*, FILE *); 
//int assignScoreCases(FILE *);
//int assignScoreControls(FILE *);
int createCoExpresionMatix(FILE *);
int assignScoreToBothControlandCases(FILE *, FILE *, FILE *, FILE *);
int assignScorePrecalculated(FILE *, FILE *);
int createCoExpresionGeneHash(FILE *);
int isConnectedPPI(int, int);
double log_N_Choose_M(int , int);
bool isSubGraphConnectedComponent(int *, int); // a list of node Ids and the number of noces
int calculateCoExpressionPValue();
float tScorePValue(float *, int);
#endif
