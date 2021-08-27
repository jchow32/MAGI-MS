# MAGI-MS
Merging affected genes into integrate networks - multiple seeds (MAGI-MS) is an extension of the existing method MAGI-S, as seen in Chow et al. 2019, Genome Medicine (https://doi.org/10.1186/s13073-019-0678-y). MAGI-MS is seed-centric method that identifies modules that consist of genes hypothesized to contribute to similar biological functions through the use of protein-protein interaction (PPI) and co-expression networks, truncating mutations observed in control populations, and up to 3 seed genes selected by the user (**Figure 1**). <br><br>

<p align="center"> <img src="https://github.com/jchow32/MAGI-MS/files/7068018/Figure1.pdf"> <br>  </p> 
<b>Figure 1. </b>General methods overview of MAGI-MS. User-selected seed gene(s), a protein-protein interaction (PPI) network, a co-expression network, and loss-of-function mutations observed in a control population are provided as input to construct modules specific to biological pathways associated with the provided seed genes. During <i>Pathway Gene Center</i> ('Pathway_GeneCenter'), scores are assigned to genes to describe their degree of co-expression with seed gene(s), and seed pathways consisting of high scoring genes are formed. During <i>Clustering</i> ('Cluster'), seed pathways are merged and refined to produce candidate modules. <br> <br>

A module is constructed by first assigning scores to every gene (<i>G</i>) in the PPI network according to their degree of co-expression with the seed gene(s) (**Equation 1**). For example, if two seed genes are provided, then a gene in the PPI network is associated with two scores; z-scoring is applied for scores calculated relative to a given seed gene. To assign a final gene score for each gene in the PPI network, the average (-avg) or minimum (-min) score is selected and the unary operation (-) is applied to the final score, such that a larger score indicates a greater degree of co-expression with the seed gene(s). <br>

<p align="center"><img width="332" alt="Equation1" src="https://user-images.githubusercontent.com/7622394/131165256-b9485871-4a1b-42a0-a748-a9c93fbf6d0b.png"> <br> </p>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- <i>H1</i>: number of pairwise comparisons of co-expression values for which co-expression of (seed gene & gene to be scored) &gt; (seed gene & and another gene in the PPI network) <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- <i>H2</i>: number of pairwise comparisons of co-expression values for which co-expression of (seed gene & gene to be scored) &gt; (gene to be scored & another gene in the PPI network) <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- <i>N</i>: total number of genes within the PPI network <br> <br>

Seed pathways are created using high scoring genes, which are genes that are highly connected in the supplied interaction networks to the input seed gene(s). Seed pathways are then merged via a random walk, and optimal modules are refined by local search (<b>Figure 1</b>). <br><br>
MAGI-MS produces a module specific to the interaction profile of the selected seed gene(s), thereby permitting the identification of modules specific to certain biological pathways relevant to the seed gene(s). Supplying MAGI-MS with a variety of candidate genes relevant to certain phenotypes can dissect comborbid phenotypes in complex disorders.

# How to run MAGI-MS
This package has two different makefiles (`magi-ms/Gene_Centric/Clusters/makefile_Pathway` and `magi-ms/Gene_Centric/Clusters/makefile_Cluster`). Compiling each makefile separately yields the executables `./Pathway_GeneCenter` and `./Cluster` .<br>
```
# How to compile makefiles
makefile -f makefile_Pathway
makefile -f makefile_Cluster
```

The first executable (`magi-ms/Gene_Centric/Paths/Pathway_GeneCenter`) create the seeds of size 5 to 8 genes. `Pathway_GeneCenter` produces 16 output files named in the format of `BestPaths.Length*.Control*.Run*`). <br>
The second executable (`magi-ms/Gene_Centric/Clusters/Cluster2`) merges the seeds into clusters to form modules.

## *Pathway_GeneCenter*

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


## *Cluster*
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
