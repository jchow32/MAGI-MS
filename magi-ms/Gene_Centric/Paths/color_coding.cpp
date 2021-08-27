#include <stdio.h>
#include "PPI_Graph.h"
#include "color_coding.h"
#include <math.h>
#include <time.h>
#include <string.h>
//int maxControlWeight;
int listOfGenesPicked[100];
int countListOfGenesPicked;
const int maxNumServeMutInControl=10;
int controlServMutAllowed=0;
bool hasColor(int A, int B)
{

	if ((A & colorBinaryRep[B])==colorBinaryRep[B])
	return true;
	else return false;
}

int addColor(int A, int B)
{
	return (A | colorBinaryRep[B]);
}

int removeColor(int A, int B)
{
	return (A & (~colorBinaryRep[B]));
}


int initialize()
{
	for (int count=1; count<numColor+1; count++)
	{
		colorBinaryRep[count]=(int) pow(2,count-1);
	}

	DPMatrix = (double ****) malloc(numColor*sizeof(int ***));


	for (int coExpresionLength=0; coExpresionLength<numColor; coExpresionLength++)
	{
		DPMatrix[coExpresionLength] = (double ***) malloc((pow(2, numColor)+1)*sizeof(double **));
	}
	
	for (int  coExpresionLength=0; coExpresionLength<numColor; coExpresionLength++)
	{
		for (int colorList=1; colorList<pow(2, numColor)+1; colorList++)
		{
			DPMatrix[coExpresionLength][colorList] = (double **) malloc(20000*sizeof(double*));
		}

	}


	
	for (int coExpresionLength=0; coExpresionLength<numColor; coExpresionLength++)
	{
		for (int colorList=1; colorList<pow(2,numColor)+1; colorList++)
		{
			for (int nodeId=0; nodeId<20000; nodeId++)
			{
				DPMatrix[coExpresionLength][colorList][nodeId] = (double *) malloc(maxNumServeMutInControl*sizeof(double));
			}
		}
	}


	for (int coExpresionLength=0; coExpresionLength<numColor; coExpresionLength++)
	{
		for (int colorList=1; colorList<pow(2,numColor)+1; colorList++)
		{
			for (int nodeId=0; nodeId<20000; nodeId++)
			{
				for (int controlSever=0; controlSever<maxNumServeMutInControl; controlSever++)
				{
					DPMatrix[coExpresionLength][colorList][nodeId][controlSever]=-1;
				}
			}
		}
	}

}


int freeMatrixes()
{
	
		
	for (int coExpresionLength=0; coExpresionLength<numColor; coExpresionLength++)
	{
		for (int colorList=1; colorList<pow(2,numColor)+1; colorList++)
		{
			for (int nodeId=0; nodeId<20000; nodeId++)
			{
				free(DPMatrix[coExpresionLength][colorList][nodeId]);
			}
		}
	}


	for (int coExpresionLength=0; coExpresionLength<numColor; coExpresionLength++)
	{
		for (int colorList=1; colorList<pow(2,numColor)+1; colorList++)
		{
			free(DPMatrix[coExpresionLength][colorList]);
		}
	}

	
	for (int coExpresionLength=0; coExpresionLength<numColor; coExpresionLength++)
	{
			free(DPMatrix[coExpresionLength]);
	}

	free(DPMatrix);


}



int randomColor()
{
	srand(time(NULL));
	for (int count=0; count<numNodes; count++)
	{
		colorAssigned[count]=((random()%(numColor))+1);
		//printf("%i\n", colorAssigned[count]);
	} 
}

bool oneColorsOrMore(int list)//returns true if the list is only consist of one color
{
	for (int count=1; count<numColor+1; count++)
	{
		if (list==colorBinaryRep[colorBinaryRep[count]])
			return true;
	}
return false;
}

double dynamicProgrammingFillMatrix(int coExpresionLength, int colorList, int nodeId, int controlServeMut)
{
double maxWeight=-2000;
double tempWeight=-2000;
int edgeCoExpresionLength;
	if (DPMatrix[coExpresionLength][colorList][nodeId][controlServeMut]!=-1)
	{
		//printf("L111\n");
		return DPMatrix[coExpresionLength][colorList][nodeId][controlServeMut];
	}

	if (colorList==colorBinaryRep[colorAssigned[nodeId]] && coExpresionLength==numColor-1 && controlServeMut==listNodes[nodeId].weightControl)
	{
		//printf("L116\n");
		DPMatrix[coExpresionLength][colorList][nodeId][controlServeMut]=listNodes[nodeId].weightCases;
		return listNodes[nodeId].weightCases;
	}		


	if (colorList==colorBinaryRep[colorAssigned[nodeId]] && coExpresionLength!=numColor-1 && controlServeMut!=listNodes[nodeId].weightControl)
	{
		//printf("L125\n");
		DPMatrix[coExpresionLength][colorList][nodeId][controlServeMut]=-2000;
		return DPMatrix[coExpresionLength][colorList][nodeId][controlServeMut];
	}
	if (oneColorsOrMore(colorList) && colorList!=colorBinaryRep[colorAssigned[nodeId]])
	{
		//printf("L131\n");
		DPMatrix[coExpresionLength][colorList][nodeId][controlServeMut]=-2000;
		return DPMatrix[coExpresionLength][colorList][nodeId][controlServeMut];
	}

	if (hasColor(colorList, colorAssigned[nodeId])==false || controlServeMut<listNodes[nodeId].weightControl)///THIS MUST BE FIXED
	//if (hasColor(colorList, colorAssigned[nodeId])==false)
	{	
		//printf("L139\n");
		DPMatrix[coExpresionLength][colorList][nodeId][controlServeMut]=-2000;
		return DPMatrix[coExpresionLength][colorList][nodeId][controlServeMut];
	}
	for (int countDegree=0; countDegree<listNodes[nodeId].degree; countDegree++)
	{
		

		//HERE CALCULATE THE edgeCoExpresionLength
		if (coExpresionMatrix[nodeId][listNodes[nodeId].neighbours[countDegree]]>coExpressHardCut)
			edgeCoExpresionLength=0;
		else edgeCoExpresionLength=1;

		if (hasColor(colorList, colorAssigned[listNodes[nodeId].neighbours[countDegree]]) && colorAssigned[listNodes[nodeId].neighbours[countDegree]]!=colorAssigned[nodeId] && coExpresionLength-edgeCoExpresionLength> (numColor-2) && controlServeMut-listNodes[nodeId].weightControl>=0)
		//if (hasColor(colorList, colorAssigned[listNodes[nodeId].neighbours[countDegree]]) && colorAssigned[listNodes[nodeId].neighbours[countDegree]]!=colorAssigned[nodeId] && controlW-listNodes[nodeId].weightControl>=0)
		{
			//printf("L155\n");
			tempWeight = dynamicProgrammingFillMatrix(coExpresionLength-edgeCoExpresionLength , removeColor(colorList, colorAssigned[nodeId]), listNodes[nodeId].neighbours[countDegree], controlServeMut-listNodes[nodeId].weightControl)+listNodes[nodeId].weightCases;
			
			if (tempWeight>maxWeight)
			{
				maxWeight=tempWeight;		
			}
		}
	}	
	

DPMatrix[coExpresionLength][colorList][nodeId][controlServeMut]=maxWeight;
return maxWeight;

}


float calculateAvgEdgeDensity()
{
float sumTotal=0;
	for (int count=0; count<countListOfGenesPicked; count++)
	{
		for (int count2=0; count2<countListOfGenesPicked; count2++)
		{
			if (count!=count2)
			{
				sumTotal=sumTotal+isConnectedPPI(listOfGenesPicked[count], listOfGenesPicked[count2]);
			}
		}
	}
//printf("%f %f\n", sumTotal, (float)sumTotal/(float)(countListOfGenesPicked*(countListOfGenesPicked-1) ));
return (float)sumTotal/(float)(countListOfGenesPicked*(countListOfGenesPicked-1));
}







float calculateAvgCoExpresionValue()
{
float sumTotal=0;
	for (int count=0; count<countListOfGenesPicked; count++)
	{
		for (int count2=0; count2<countListOfGenesPicked; count2++)
		{
			if (count!=count2)
			{
				sumTotal=sumTotal+coExpresionMatrix[listOfGenesPicked[count]][listOfGenesPicked[count2]];
			}
		}
	}
return (float)sumTotal/(float)(countListOfGenesPicked*(countListOfGenesPicked-1));
}

int traverseBack(int coExpresionLength, int colorList, int nodeId, int controlServeMut, FILE *fp)
{
int edgeCoExpresionLength=0;
	//fprintf(fp,"%i %i\n", colorList, colorBinaryRep[colorAssigned[nodeId]]);
	if (colorList==colorBinaryRep[colorAssigned[nodeId]] && coExpresionLength==numColor-1)
	{
		//DPMatrix[controlW][colorList][nodeId]=listNodes[nodeId].weightCases;
		listOfGenesPicked[countListOfGenesPicked]=nodeId;
		fprintf(fp,"%lf %s\n",listNodes[nodeId].weightCases, listNodes[nodeId].nodeName);
		countListOfGenesPicked++;
		return 1;
	}


	for (int countDegree=0; countDegree<listNodes[nodeId].degree; countDegree++)
	{
			
		if (coExpresionMatrix[nodeId][listNodes[nodeId].neighbours[countDegree]]>coExpressHardCut)
			edgeCoExpresionLength=0;
		else edgeCoExpresionLength=1;
//	if (hasColor(colorList, colorAssigned[listNodes[nodeId].neighbours[countDegree]]))
	//	fprintf(fp,"L269 %i %i %i\n", edgeCoExpresionLength, controlServeMut-listNodes[nodeId].weightControl, controlServeMut);
		if (hasColor(colorList, colorAssigned[listNodes[nodeId].neighbours[countDegree]]) && colorAssigned[listNodes[nodeId].neighbours[countDegree]]!=colorAssigned[nodeId] && coExpresionLength-edgeCoExpresionLength> (numColor-2) &&  controlServeMut-listNodes[nodeId].weightControl>=0)
		{
	//fprintf(fp,"L272\n");
			if (DPMatrix[coExpresionLength][colorList][nodeId][controlServeMut]==DPMatrix[coExpresionLength-edgeCoExpresionLength][removeColor(colorList, colorAssigned[nodeId])][listNodes[nodeId].neighbours[countDegree]][controlServeMut-listNodes[nodeId].weightControl]+listNodes[nodeId].weightCases)
			{
				listOfGenesPicked[countListOfGenesPicked]=nodeId;
				countListOfGenesPicked++;
				//fprintf(fp,"%lf %s\n", listNodes[nodeId].weightCases, listNodes[nodeId].nodeName);
				if (traverseBack(coExpresionLength-edgeCoExpresionLength, removeColor(colorList, colorAssigned[nodeId]), listNodes[nodeId].neighbours[countDegree], controlServeMut-listNodes[nodeId].weightControl, fp))
				{		
					fprintf(fp,"%lf %s\n", listNodes[nodeId].weightCases, listNodes[nodeId].nodeName);
					return 1;
				}
				else countListOfGenesPicked--;
			}
		}

	}
return 0;

		
}





int numberColor(int list)//calculates number of colors in each list
{
int countBinary=0;
	for (int count=0; count<numColor+1; count++)
	{
		if (list%2==1)
			countBinary++;	
		list=list >> 1;	
	}
return countBinary;
}

int reinitializeMatrix()
{
	
	for (int coExpresionLength=0; coExpresionLength<numColor; coExpresionLength++)
	{
		for (int colorList=1; colorList<pow(2,numColor)+1; colorList++)
		{
			for (int nodeId=0; nodeId<20000; nodeId++)
			{
				for (int controlSever=0; controlSever<maxNumServeMutInControl; controlSever++)
				{
					DPMatrix[coExpresionLength][colorList][nodeId][controlSever]=-1;
				}
			}
		}
	}
}


int runColorCodingMethod(FILE *fp)
{
double maxScore=-2000;
//int minControl=0;
int minLength=0;

//int bestControlWeight;
int bestColorList;
int bestNodeId;
int minCoExprLength;
int bestControlSeverMut=0;
int redo=1000;
	while(redo>0)
	{
	randomColor();
	reinitializeMatrix();
	maxScore=-2000;
	//minControl=0;
	minLength=0;
	//bestControlWeight=0;
	bestColorList=0;
	bestNodeId=0;
	minCoExprLength=0;
	//for (int coExpresionLength=0; coExpresionLength<numColor; coExpresionLength++)
	int coExpresionLength=numColor-1;
	{
		for (int colorList=1; colorList<pow(2,numColor)+1; colorList++)
		{
			for (int nodeId=0; nodeId<numNodes; nodeId++)
			{
				for (int controlSeverMut=0; controlSeverMut<controlServMutAllowed; controlSeverMut++)
				{
				
					DPMatrix[coExpresionLength][colorList][nodeId][controlSeverMut]=dynamicProgrammingFillMatrix(coExpresionLength, colorList, nodeId, controlSeverMut);
					
				//	printf("%i %i %i %lf %lf\n", controlSeverMut, colorList, nodeId, DPMatrix[coExpresionLength][colorList][nodeId][controlSeverMut], maxScore);
					if (maxScore<DPMatrix[coExpresionLength][colorList][nodeId][controlSeverMut] && coExpresionLength> 1*(numColor-2) && numberColor(colorList)==numColor)
					{ 
				//		printf("%i %i %i %lf %lf\n", controlSeverMut, colorList, nodeId, DPMatrix[coExpresionLength][colorList][nodeId][controlSeverMut], maxScore);
						maxScore=DPMatrix[coExpresionLength][colorList][nodeId][controlSeverMut];
						//traverseBack(controlWeight, colorList, nodeId);
						//minControl=controlWeight;
						bestControlSeverMut=controlSeverMut;
						minLength=numberColor(colorList);
						minCoExprLength=coExpresionLength;
						//bestControlWeight=controlWeight;
						bestColorList=colorList;
						bestNodeId=nodeId;
					}/*else if (maxScore==DPMatrix[controlWeight][colorList][nodeId] && minControl>controlWeight)
					{
						
						//printf("%i %i %i %f\n", controlWeight, colorList, nodeId, DPMatrix[controlWeight][colorList][nodeId]);
						maxScore=DPMatrix[controlWeight][colorList][nodeId];
						//traverseBack(controlWeight, colorList, nodeId);
						minControl=controlWeight;
						minLength=numberColor(colorList);
						bestControlWeight=controlWeight;
						bestColorList=colorList;
						bestNodeId=nodeId;
					}else if (maxScore==DPMatrix[controlWeight][colorList][nodeId] && minControl==controlWeight && numberColor(colorList)<minLength)
					{
						
						//printf("%i %i %i %f\n", controlWeight, colorList, nodeId, DPMatrix[controlWeight][colorList][nodeId]);
						maxScore=DPMatrix[controlWeight][colorList][nodeId];
						//traverseBack(controlWeight, colorList, nodeId);
						minControl=controlWeight;
						minLength=numberColor(colorList);
						
						bestControlWeight=controlWeight;
						bestColorList=colorList;
						bestNodeId=nodeId;
					}*/
					/*else{
						
						printf("WTF %i %i %i %f\n", controlWeight, colorList, nodeId, DPMatrix[controlWeight][colorList][nodeId]);
						traverseBack(controlWeight, colorList, nodeId);
					}*/

				}
			}
		}
	}
		countListOfGenesPicked=0;
		if (traverseBack(minCoExprLength, bestColorList, bestNodeId, bestControlSeverMut, fp))
		{
			
			fprintf(fp,"recolor %i\n", redo);
			redo--;
			fprintf(fp,"%i %i %i %f\n", numberColor(bestColorList),  minCoExprLength, bestControlSeverMut,  DPMatrix[minCoExprLength][bestColorList][bestNodeId][bestControlSeverMut]);
		//traverseBack(minCoExprLength, bestColorList, bestNodeId, bestControlSeverMut, fp);
			fprintf(fp,"Avg CoExpr:%f %f %f\n", calculateAvgCoExpresionValue(), calculateAvgEdgeDensity(), DPMatrix[minCoExprLength][bestColorList][bestNodeId][bestControlSeverMut]);
		}

	}
}



int main(int argv, char *argc[])
{
	int randomRunId=0;
	char fileName[100];
	int avg_bool = 0; // TODO initialize the bool showing if you should do average or maximum score for multiple seeds
	FILE *fp=fopen(argc[1],"r");//PPI Network
	FILE *fp2=fopen(argc[2],"r");//The name of genes as center (such as FMR1)
	FILE *fp3=fopen(argc[3],"r");//The control
	FILE *fp4=fopen(argc[4],"r");//The gene Length
	//maxControlWeight=atoi(argc[5]);
	FILE *fp5=fopen(argc[5],"r");//The gene coexpresions
	FILE *fp6=fopen(argc[6],"r");
	randomRunId = atoi(argc[7]);
	sprintf(fileName,"RandomGeneList.%i\0", randomRunId);
	FILE *fp7=fopen(fileName,"w");
	// TODO you should accept an argument from the user which specifies if will use the average or maximum score for multiple seeds
	if (strcmp(argc[8], "-avg")==0)  // then do the average
	{
		avg_bool = 1;
	}
	//numColor=atoi(argc[7]);//
	//controlServMutAllowed=atoi(argc[8]);
	createPPI_Graph(fp2, fp);
	fclose(fp2);
	fp2=fopen(argc[2],"r");
	printf("L446\n");
	createCoExpresionGeneHash(fp5);
	printf("L448\n");
	createCoExpresionMatix(fp6);
	printf("L450\n");
	assignScoreToBothControlandCases(fp2, fp3, fp4, fp7, avg_bool);
	//for (int countNumColor=8; countNumColor<9; countNumColor++)
	printf("L453\n");
	for (int countNumColor=5; countNumColor<9; countNumColor++)
	{
		numColor=countNumColor;
		for (int countTrunMut=1; countTrunMut<5; countTrunMut++)
		//for (int countTrunMut=7; countTrunMut<8; countTrunMut++)
		{
			sprintf(fileName,"BestPaths.Length%i.Control%i.Run%i\0", countNumColor, countTrunMut, randomRunId);
			fp7=fopen(fileName, "w");
			controlServMutAllowed=countTrunMut;
			initialize();
			runColorCodingMethod(fp7);
			freeMatrixes();
			fclose(fp7);
		}
	}
}

