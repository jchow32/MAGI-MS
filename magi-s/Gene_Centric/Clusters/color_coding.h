#ifndef _color_coding
#define _color_coding
#include "PPI_Graph.h"

int numColor; // the is the length of the path we are looking to find(colors start from 1)
const int maxNumColor=20;
const int weightContorl=100; // maximum weight of control allowd (maximum size of knapsack
//const float coExpressHardCut=0.37;

double ****DPMatrix;
int colorAssigned[maxNumNode];

int colorBinaryRep[maxNumColor+1];//color 0 is not there. Starts by color1

int initialize();
bool hasColor(int, int); // given the color list A (an integer which repersents the color set). check if a partocular color B is there or not.
int addColor(int, int);// given the color list A, add color B. return the new list
int removeColor(int, int);// given the color list A, remove color B. return the new list
int randomColor();
#endif
