#ifndef CLUSTERIN_PHASE
#define CLUSTERIN_PHASE
#include "PPI_Graph.h"
//0.467581 0.16842
float minAvgExpr=0.415;
float minAvgExpr2=0.01;
//float minAvgEdgeDensity=0.116734;
float minAvgEdgeDensity=0.085;
int upperLimitSize; // we don't want clusters bigger than this size (it is an upper limit)
int lowerLimitSize;
float minAvgExprPerSize[100];
float minAvgEdgeDensityPerSize[100];

const int maxClusterSize=50000;
int minClusterWeight;
typedef struct clustersSelected{
	int nodeId[200];
	//int nodeIdCount[200];
	int sizeCluster;
	float totalScore;//-log(p) : for each set of nodes you calculate it 
	float seedScore; // thr -log(P) score form seeds (summation of -log(p))
	float averageCoexpresion;//average pairwise coexpresion
	float averageEdgeDensity;
	
	int numSevereMutInCases;
	int numMissenseMutInCases;
	int numServeMutInControl;
	bool used;
}clustersSelected;

typedef struct adjClusterList{
	int *adjClusterIds;
	int numberAdjClusters;
}adjClusterList;

adjClusterList *clusterGraphAdjRep;

int **clusterGraphMatrixRep;

clustersSelected listClusterSelected[maxClusterSize];
int totalClusters;
char potNewPath[100][geneNameLen];

#endif
