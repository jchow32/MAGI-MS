# MAGI-S
Merging affected genes into integrate networks - seed (MAGI-S) is a seed-centric method that identifies modules that consist of genes hypothesized to contribute to similar biological functions through the use of protein-protein interaction and co-expression networks, mutations observed in control populations, and a single seed gene selected by the user. A module is constructed combinatorially by first generating seed pathways, which contain genes that are highly connected in the supplied interaction networks to the input seed gene. Seed pathways are then merged via a random walk, and optimal modules are refined by local search. MAGI-S produces a module specific to the interaction profile of the selected seed gene. Supplying MAGI-S with a variety of candidate genes relevant to certain phenotypes can dissect comborbid phenotypes in complex disorders, as seen in https://doi.org/10.1186/s13073-019-0678-y (Chow et al. 2019, Genome Medicine).


# How to run MAGI-S
This package has two different makefiles (magi-s/Gene_Centric/Clusters/makefile_Pathway and magi-s/Gene_Centric/Clusters/makefile_Cluster). Both need to be used to build it. 
makefile -f makefile_Pathway
makefile -f makefile_Cluster

The two executable files are:
./Pathway_GeneCenter
./Cluster2

The first executable (Gene_Centric/Paths/Pathway_GeneCenter) create the seeds of size 5 to 8 genes (the names of the output files are in format of BestPaths.Length*.Control*.Run*).
The second executable (Gene_Centric/Clusters/Cluster2) merges the seeds into clusters/modules.

*Pathway_GeneCenter*

PPI Networks, List of de novo mutations in cases, list of mutations in controls, length of genes (genes which length is not provided are given a default length of 3300bp), hash table of the gene names to exact coexpression values, pair-wise gene coexpressions, and a list if gene to filter/remove (optional)

	Required Parameters 
	$PPI <PPI Network> 
		The input file should be a file of binary interactions. File input format:<Gene_Name_1><\t><Gene_Name_2> (Example: String_HPRD_PPI) 
	$caseGene <seed gene>
		Seed gene. File input format:<Gene_Name> (Example:file containing one line with gene name, like SCN1A) 
	$geneCoexpressionID <Gene CoExpression Id>
		The input gives the order of each gene appearing in the coExpression table. FIle input format: <Gene_row_id><\t>OTHER_ID<\t><Gene_Name> (Example: GeneCoExpresion_ID)
	$coexpressionMatrix <CoExpression Matrix>
		Pairwise gene coexpression values. File input format:<Gene_Name_1><\t><Gene_Name_2><\t><CoExpression - r^6> (example: adj1.csv.Tab.BinaryFormat)
		Note that the genes pairs are sorted based on their <Gene_row_id> provided in the file inputed as -h parameters. 
	$controlMutationList <control mutation list>
		The number of mutations in each gene in controls. File input format:<Gene_Name_1><\t><number of mutations in control> (example: New_ESP_Sereve)
	$geneLength <Length of genes>
		This file includes the length of each gene. File input format:<Gene Name><\t><Length> (example : Gene_Name_Length)
	$ID <run id>
		The id for this run (integer)
	
Output:
1) RandomGeneList.%i: The list of genes with their assigned score (based on mutations in cases) and total number of mutation in controls 
2) BestPaths.Length%i.Control%i.Run%i: The seeds created (a total of 1000 seed per type), for different lengths (default is from 5 to 8) and total mutations in control (ranging from 0 to 4)


*Cluster2*
	Required Parameters

	-p <PPI Network> $PPI
		The input file should be a file of binary interactions. File input format:<Gene_Name_1><\t><Gene_Name_2> (Example: String_HPRD_PPI)
	-c <case/control mutation scores>
		The input file is the RandomGeneList.%i file created in previous step
	-h <Gene CoExpression Id> $geneCoexpressionID
 		The input gives the order of each gene appearing in the coExpression table. FIle input format: <Gene_row_id><\t>OTHER_ID<\t><Gene_Name> (Example: GeneCoExpresion_ID)
	-e <CoExpression Matrix> $coexpressionMatrix
		Pairwise gene coexpression values. File input format:<Gene_Name_1><\t><Gene_Name_2><\t><CoExpression - r^6> (example: adj1.csv.Tab.BinaryFormat)
		Note that the genes pairs are sorted based on their <Gene_row_id> provided in the file inputed as -h parameters.
	-s <Seed File> 
		The file names of different seeds. <Seed Name File><\t><Number of seeds><\t><Length of Seeds> (example: Paths. The Paths file references files (BestPaths.Length%i.Control%i.Run%i) created from running Pathway_GeneCenter in the previous step)
	-m <upper bound on control mutations>
		The total number of mutations in control's allowed. 
	-l <lower bound on the size of the module>
		The minimum number of genes in the module 
	-u <upper bound on the size of the module>
		The maximum number of genes in the module
	-a <minimum ratio of seed score allowed>
		For each seed type the ratio of the score from maximum score of the seed allowed (in the paper 0.5 was used)
	-i <run id>  
		The id for this run.

	Optional Parameters
	
	-minCoExpr <minimum pair-wise coexpression value>
		The minimum pair-wise coexpression value per gene allowed (the default is 0.01, i.e. r^2>0.01, which is the median coexpression value in the input adj1.csv.Tab.BinaryFormat)
	-avgCoExpr <minimum average coexpression of the module>
		The minimum average coexpression of the modules allowed (the default is 0.415)
	-avgDensity <minimum density of PPI>
		The minimum PPI density of the modules allowed (the default is 0.08)

Example: 
Gene_Centric/Clusters/Cluster2 -p $PPI -c ../RandomGeneList.1 -h $geneCoexpressionID -e $coexpressionMatrix -s Paths -m 6 -l 20 -u 100 -a 0.5 -i 1 -minCoExpr 0.01 -avgCoExpr 0.440000 -avgDensity 0.14 > CCluster_Len_20_MinCoExpr0.01_LessCoexpr_New_0.440000

# Tips
In practice, Cluster2 can be run multiple times on the same BestPaths* generated from a single run of Pathway_GeneCenter, and the highest scoring module can be retrieved from all runs of Cluster2. 
36 example input parameters (a single combination of parameters per line) for Cluster2 are displayed in the file magi-s/Inputs, in which the parameters -l, -i, -avgCoExpr, and -avgDensity are varied. 36 example output names are also provided in the file magi-s/Outputs2, which correspond to the suggested input parameters listed in magi-s/Inputs. Note that the first line of magi-s/Inputs and magi-s/Outputs2 are substituted into the Example Cluster2 command shown above. 
