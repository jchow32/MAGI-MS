#include <stdio.h>
#include "PPI_Graph.h"
#include "color_coding.h"
#include <math.h>
#include <time.h>
#include <string.h>
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
double maxWeight=-2000;//-2000 is equal to - infinity
double tempWeight=-2000;
int edgeCoExpresionLength;
	if (DPMatrix[coExpresionLength][colorList][nodeId][controlServeMut]!=-1)
	{
		return DPMatrix[coExpresionLength][colorList][nodeId][controlServeMut];
	}

	if (colorList==colorBinaryRep[colorAssigned[nodeId]] && coExpresionLength==numColor-1 && controlServeMut==listNodes[nodeId].weightControl)
	{
		DPMatrix[coExpresionLength][colorList][nodeId][controlServeMut]=listNodes[nodeId].weightCases;
		return listNodes[nodeId].weightCases;
	}		


	if (colorList==colorBinaryRep[colorAssigned[nodeId]] && coExpresionLength!=numColor-1 && controlServeMut!=listNodes[nodeId].weightControl)
	{
		DPMatrix[coExpresionLength][colorList][nodeId][controlServeMut]=-2000;
		return DPMatrix[coExpresionLength][colorList][nodeId][controlServeMut];
	}
	if (oneColorsOrMore(colorList) && colorList!=colorBinaryRep[colorAssigned[nodeId]])
	{
		DPMatrix[coExpresionLength][colorList][nodeId][controlServeMut]=-2000;
		return DPMatrix[coExpresionLength][colorList][nodeId][controlServeMut];
	}

	if (hasColor(colorList, colorAssigned[nodeId])==false || controlServeMut<listNodes[nodeId].weightControl)///THIS MUST BE FIXED
	{	
		DPMatrix[coExpresionLength][colorList][nodeId][controlServeMut]=-2000;
		return DPMatrix[coExpresionLength][colorList][nodeId][controlServeMut];
	}
	for (int countDegree=0; countDegree<listNodes[nodeId].degree; countDegree++)
	{
		

		//HERE CALCULATE THE edgeCoExpresionLength
		if (coExpresionMatrix[nodeId][listNodes[nodeId].neighbours[countDegree]]>minCoExpresThreshold)
			edgeCoExpresionLength=0;
		else edgeCoExpresionLength=1;

		if (hasColor(colorList, colorAssigned[listNodes[nodeId].neighbours[countDegree]]) && colorAssigned[listNodes[nodeId].neighbours[countDegree]]!=colorAssigned[nodeId] && coExpresionLength-edgeCoExpresionLength> (numColor-2) && controlServeMut-listNodes[nodeId].weightControl>=0)
		{
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
	if (colorList==colorBinaryRep[colorAssigned[nodeId]] && coExpresionLength==numColor-1)
	{
		listOfGenesPicked[countListOfGenesPicked]=nodeId;
		fprintf(fp,"%lf %s\n",listNodes[nodeId].weightCases, listNodes[nodeId].nodeName);
		countListOfGenesPicked++;
		return 1;
	}


	for (int countDegree=0; countDegree<listNodes[nodeId].degree; countDegree++)
	{
			
		if (coExpresionMatrix[nodeId][listNodes[nodeId].neighbours[countDegree]]>minCoExpresThreshold)
			edgeCoExpresionLength=0;
		else edgeCoExpresionLength=1;
		if (hasColor(colorList, colorAssigned[listNodes[nodeId].neighbours[countDegree]]) && colorAssigned[listNodes[nodeId].neighbours[countDegree]]!=colorAssigned[nodeId] && coExpresionLength-edgeCoExpresionLength> (numColor-2) &&  controlServeMut-listNodes[nodeId].weightControl>=0)
		{
			if (DPMatrix[coExpresionLength][colorList][nodeId][controlServeMut]==DPMatrix[coExpresionLength-edgeCoExpresionLength][removeColor(colorList, colorAssigned[nodeId])][listNodes[nodeId].neighbours[countDegree]][controlServeMut-listNodes[nodeId].weightControl]+listNodes[nodeId].weightCases)
			{
				listOfGenesPicked[countListOfGenesPicked]=nodeId;
				countListOfGenesPicked++;
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
double maxScore=-2000; //-2000 is equal to -infinity
//int minControl=0;
int minLength=0;

//int bestControlWeight;
int bestColorList;
int bestNodeId;
int minCoExprLength;
int bestControlSeverMut=0;
int redo=1000; //number of times to run color-coding (number of paths outputed)
	while(redo>0)
	{
	randomColor();
	reinitializeMatrix();
	maxScore=-2000;
	minLength=0;
	bestColorList=0;
	bestNodeId=0;
	minCoExprLength=0;
	int coExpresionLength=numColor-1;
	{
		for (int colorList=1; colorList<pow(2,numColor)+1; colorList++)
		{
			for (int nodeId=0; nodeId<numNodes; nodeId++)
			{
				for (int controlSeverMut=0; controlSeverMut<controlServMutAllowed; controlSeverMut++)
				{
				
					DPMatrix[coExpresionLength][colorList][nodeId][controlSeverMut]=dynamicProgrammingFillMatrix(coExpresionLength, colorList, nodeId, controlSeverMut);
					
					if (maxScore<DPMatrix[coExpresionLength][colorList][nodeId][controlSeverMut] && coExpresionLength> 1*(numColor-2) && numberColor(colorList)==numColor)
					{ 
						maxScore=DPMatrix[coExpresionLength][colorList][nodeId][controlSeverMut];
						bestControlSeverMut=controlSeverMut;
						minLength=numberColor(colorList);
						minCoExprLength=coExpresionLength;
						bestColorList=colorList;
						bestNodeId=nodeId;
					}

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
			fprintf(fp,"Avg CoExpr:%f %f %f\n", calculateAvgCoExpresionValue(), calculateAvgEdgeDensity(), DPMatrix[minCoExprLength][bestColorList][bestNodeId][bestControlSeverMut]);
		}

	}
}



int main(int argv, char *argc[])
{

	int randomRunId=0;
	char fileName[100];

	FILE *fp, *fp2, *fp3, *fp4, *fp5, *fp6, *fp7;
	int countNumParamters=0;
	bool filter=false;


	for (int count=0; count<argv; count++)
	{
		if (strcmp(argc[count], "-p")==0) // PPI
			{
				fp=fopen(argc[count+1],"r");
				countNumParamters++;
			}
		if (strcmp(argc[count], "-c")==0) // CASES
			{
				fp2=fopen(argc[count+1],"r");
				countNumParamters++;
			}
		if (strcmp(argc[count], "-d")==0) // CONTROL : ESP
			{
				fp3=fopen(argc[count+1],"r");
				countNumParamters++;
			}
		if (strcmp(argc[count], "-l")==0) // Gene Length
			{
				fp4=fopen(argc[count+1],"r");
				countNumParamters++;
			}
		if (strcmp(argc[count], "-h")==0) // Gene coExpression Hash Table
			{
				fp5=fopen(argc[count+1],"r");
				countNumParamters++;
			}
		if (strcmp(argc[count], "-e")==0) // CoExpression Matrix
			{
				fp6=fopen(argc[count+1],"r");
				countNumParamters++;
			}
		if (strcmp(argc[count], "-i")==0) // ranodm id run
			{
				randomRunId = atoi(argc[count+1]);
				countNumParamters++;
			}
		if (strcmp(argc[count], "-f")==0) // genes to filter out (by assiging a very high ESP/Control mutations to it.
			{
				fp7=fopen(argc[count+1], "r");
				filter=true; 
			}
	}

	if (countNumParamters!=7)
	{
		printf("usage:\n needs seven paramters to run\n -p <PPI Network> \n -c <cases mutaion list> \n -d <control mutation list> \n -l <Length of each genes> \n -h <Gene CoExpression Id> \n -e <CoExpression Matrix> \n -i <random run id> \n -f <filter genes (optional, defualt is no gene is filtered) \n");
		return 0;
	}

//	FILE *fp=fopen(argc[1],"r");//PPI Network
//	FILE *fp2=fopen(argc[2],"r");//The cases
//	FILE *fp3=fopen(argc[3],"r");//The control
//	FILE *fp4=fopen(argc[4],"r");//The gene Length
	//maxControlWeight=atoi(argc[5]);
//	FILE *fp5=fopen(argc[5],"r");//The gene coexpresions
//	FILE *fp6=fopen(argc[6],"r");
//	randomRunId = atoi(argc[7]);
	sprintf(fileName,"RandomGeneList.%i\0", randomRunId);
	FILE *fp8=fopen(fileName,"w");
	//numColor=atoi(argc[7]);//
	//controlServMutAllowed=atoi(argc[8]);
	createPPI_Graph(fp);
	assignScoreToBothControlandCases(fp2, fp3, fp4, fp7, fp8, filter);
	createCoExpresionGeneHash(fp5);
	createCoExpresionMatix(fp6);

	for (int countNumColor=5; countNumColor<9; countNumColor++)
	{
		numColor=countNumColor;
		for (int countTrunMut=1; countTrunMut<5; countTrunMut++)
		{
			sprintf(fileName,"BestPaths.Length%i.Control%i.Run%i\0", countNumColor, countTrunMut, randomRunId);
			fp8=fopen(fileName, "w");
			controlServMutAllowed=countTrunMut;
			initialize();
			runColorCodingMethod(fp8);
			freeMatrixes();
			fclose(fp8);
		}
	}
}

