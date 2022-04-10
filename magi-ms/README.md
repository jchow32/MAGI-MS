# MAGI-MS
Merging affected genes into integrate networks - multiple seeds (MAGI-MS) is an extension of the existing method MAGI-S, as seen in Chow et al. 2019, Genome Medicine (https://doi.org/10.1186/s13073-019-0678-y). MAGI-MS is seed-centric method that identifies modules that consist of genes hypothesized to contribute to similar biological functions through the use of protein-protein interaction (PPI) and co-expression networks, truncating (loss-of-function, LOF) mutations observed in control populations, and up to 3 seed genes selected by the user (**Figure 1**). <br><br>

<p align="center"> <img src="https://github.com/jchow32/MAGI-MS/files/7068018/Figure1.pdf"> <br>  </p> 
<b>Figure 1. </b>General methods overview of MAGI-MS. User-selected seed gene(s), a protein-protein interaction (PPI) network, a co-expression network, and loss-of-function mutations observed in a control population are provided as input to construct modules specific to biological pathways associated with the provided seed genes. During <i>Pathway Gene Center</i> ('Pathway_GeneCenter'), scores are assigned to genes to describe their degree of co-expression with seed gene(s), and seed pathways consisting of high scoring genes are formed. During <i>Clustering</i> ('Cluster2'), seed pathways are merged and refined to produce candidate modules. <br> <br>

A module is constructed by first assigning scores to every gene (<i>G</i>) in the PPI network according to their degree of co-expression with the seed gene(s) (**Equation 1**). For example, if two seed genes are provided, then a gene in the PPI network is associated with two scores; z-scoring is applied for scores calculated relative to a given seed gene. To assign a final gene score for each gene in the PPI network, the average (-avg) or minimum (-min) score is selected, such that a larger score indicates a greater degree of co-expression with the seed gene(s). <br>

<p align="center"><i> G<sub>s,i</sub> = ((H<sub>1</sub>)(H<sub>2</sub>)) / N<sup>2</sup> </i><br> </p>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- <i>H<sub>1</sub></i>: number of pairwise comparisons of co-expression values for which co-expression of (seed gene & gene to be scored) &gt; (seed gene & and another gene in the PPI network) <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- <i>H<sub>2</sub></i>: number of pairwise comparisons of co-expression values for which co-expression of (seed gene & gene to be scored) &gt; (gene to be scored & another gene in the PPI network) <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- <i>N</i>: total number of genes within the PPI network <br> <br>

Seed pathways are created using high scoring genes, which are genes that are highly connected in the supplied interaction networks to the input seed gene(s). Seed pathways are then merged via a random walk, and optimal modules are refined by local search (<b>Figure 1</b>). <br><br>
MAGI-MS produces a module specific to the interaction profile of the selected seed gene(s), thereby permitting the identification of modules specific to certain biological pathways relevant to the seed gene(s). Supplying MAGI-MS with a variety of candidate genes relevant to certain phenotypes can dissect comborbid phenotypes in complex disorders.

# How to run MAGI-MS
This package has two different makefiles (`magi-ms/Gene_Centric/Paths/makefile` and `magi-ms/Gene_Centric/Clusters/makefile_Cluster2`). The [Pathway_GeneCenter](https://github.com/jchow32/MAGI-MS/blob/main/magi-ms/Gene_Centric/Paths/Pathway_GeneCenter) and [Cluster2](https://github.com/jchow32/MAGI-MS/blob/main/magi-ms/Gene_Centric/Clusters/Cluster2) programs have been compiled from these makefiles and are available for immediate use. <br>

The first executable (`magi-ms/Gene_Centric/Paths/Pathway_GeneCenter`) creates the seed pathways of size 5 to 8 genes. `Pathway_GeneCenter` produces 16 output files named in the format of `BestPaths.Length*.Control*.Run*`). <br>
The second executable (`magi-ms/Gene_Centric/Clusters/Cluster2`) merges the seed pathways written by `Pathway_GeneCenter` into clusters to form modules.

## *Pathway_GeneCenter*

	Required Parameters 
	$PPI <PPI network> 
		The input file should be a file of binary interactions. File input format:<Gene_Name_1><\t><Gene_Name_2> (Example: StringNew_HPRD) 
	$caseGene <Seed gene>
		Seed gene(s), separated by new line characters in a text file. File input format:<Gene_Name_1>\n<Gene_Name_2> (Example: caseGene) 
	$controlMutationList <Control Mutation List>
		The number of truncating (loss-of-function) mutations in each gene in a control population. File input format:<Gene_Name_1><\t><Number_LOF_in_Controls> (example: New_ESP_Sereve)
	$geneLength <Length of genes>
		This file includes the length of each gene. File input format:<Gene_Name><\t><Length> (example : Gene_Name_Length)
	$geneCoexpressionID <Gene CoExpression ID>
		The input gives the order of each gene appearing in the coExpression table. File input format: <Gene_row_id><\t>OTHER_ID<\t><Gene_Name> (Example: GeneCoExpresion_ID)
	$coexpressionMatrix <CoExpression Matrix>
		Pairwise gene co-expression values. File input format:<Gene_Name_1><\t><Gene_Name_2><\t><CoExpression - r^6> (Example: adj1.csv.Tab.BinaryFormat)
		Note: Genes pairs are sorted based on their <Gene_row_id> provided in $geneCoexpressionID. 
	$ID <run id>
		The id for this run (integer). This integer is used to seed a random number generator. 
	$scoring <Final gene scoring method>
		This string indicates how scores will be assigned to genes. String options: '-avg' or '-min' (without apostrophes). 
		Note: By default, the minimum score (which maximizes the final score assigned to a gene), rather than the average score, is selected. 
	
Output:
1) `RandomGeneList.%i`: A text file listing genes with their assigned score. File output format: <Gene_Name> <Score> <0> <0> <Number_LOF_in_Controls> <0>
2) `BestPaths.Length%i.Control%i.Run%i`: The seed pathways created, for different lengths (default is from 5 to 8) and total mutations in control (ranging from 0 to 4), with 1,000 iterations per combination of length and total mutations allowed in controls, resulting in a total of 16,000 seed pathways. A total of 16 files, each with 1,000 seed pathways, are created.
	
	
Example downloadable input files: 
* PPI network: [StringNew_HPRD](https://eichlerlab.gs.washington.edu/MAGI/Data/StringNew_HPRD)
* Seed gene(s): [caseGene](https://github.com/jchow32/MAGI-MS/blob/main/magi-ms/caseGene)
* Gene CoExpression ID: [GeneCoExpresion_ID](https://eichlerlab.gs.washington.edu/MAGI/Data/GeneCoExpresion_ID)
* CoExpression Matrix: [adj1.csv.Tab.BinaryFormat](https://eichlerlab.gs.washington.edu/MAGI/Data/adj1.csv.Tab.BinaryFormat)
* Control Mutation List: [New_ESP_Sereve](https://eichlerlab.gs.washington.edu/MAGI/Data/New_ESP_Sereve)
* Length of genes: [Gene_Name_Length](https://eichlerlab.gs.washington.edu/MAGI/Data/Gene_Name_Length)

Example command:
```
magi-ms/Gene_Centric/Paths/Pathway_GeneCenter \
StringNew_HPRD \
caseGene \
New_ESP_Sereve \
Gene_Name_Length \
GeneCoExpresion_ID \
adj1.csv.Tab.BinaryFormat \
1 \
-min
```


## *Cluster2*
	Required Parameters
	-p <PPI Network> $PPI
		The input file should be a file of binary interactions. File input format:<Gene_Name_1><\t><Gene_Name_2> (Example: StringNew_HPRD)
	-c <case/control mutation scores>
		The input file is the RandomGeneList.%i file created in previous step.
	-h <Gene CoExpression ID> $geneCoexpressionID
 		The input gives the order of each gene appearing in the coExpression table. File input format: <Gene_row_id><\t>OTHER_ID<\t><Gene_Name> (Example: GeneCoExpresion_ID)
	-e <CoExpression Matrix> $coexpressionMatrix
		Pairwise gene coexpression values. File input format:<Gene_Name_1><\t><Gene_Name_2><\t><CoExpression - r^6> (example: adj1.csv.Tab.BinaryFormat)
		Note that the genes pairs are sorted based on their <Gene_row_id> provided in the file inputed as -h parameter.
	-s <Seed File> 
		The file names of different seed pathways. <Seed Name File><\t><Number of seeds><\t><Length of Seeds> (example: Paths. The Paths file references files created from running Pathway_GeneCenter in the previous step)
	-m <upper bound on control mutations>
		The total number of mutations in controls allowed (recommended value: 6).
	-l <lower bound on the size of the module>
		The minimum number of genes in the module.
	-u <upper bound on the size of the module>
		The maximum number of genes in the module.
	-a <minimum ratio of seed score allowed>
		For each seed type, the ratio of the score from the maximum score of the seed allowed (recommended value: 0.5)
	-i <run id>  
		The id for this run (integer).

	Optional Parameters
	
	-minCoExpr <minimum pair-wise co-expression value>
		The minimum pair-wise co-expression value per gene allowed (the default is 0.01, i.e. r^2>0.01, which is the median co-expression value in the input adj1.csv.Tab.BinaryFormat)
	-avgCoExpr <minimum average co-expression of the module>
		The minimum average co-expression of the modules allowed (the default is 0.415, recommended range: 0.425-0.52)
	-avgDensity <minimum density of PPI>
		The minimum PPI density of the modules allowed (the default is 0.08, recommended range 0.080-0.14)

	
Example downloadable input files:
* Seed pathway file: [Paths](https://github.com/jchow32/MAGI-MS/blob/main/magi-ms/Paths)

Example command: 
```
magi-ms/Gene_Centric/Clusters/Cluster2 \
-p StringNew_HPRD \
-c RandomGeneList.1 \
-h GeneCoExpresion_ID \
-e adj1.csv.Tab.BinaryFormat \
-s Paths \
-m 6 -l 20 -u 100 -a 0.5 \
-i 1 \
-minCoExpr 0.01 -avgCoExpr 0.440000 -avgDensity 0.14 > CCluster_Len_20_MinCoExpr0.01_LessCoexpr_New_0.440000
```
	
## Tips
* For any given output file generated by <i>Cluster2</i> (for example, called `$output_file_name`), the highest scoring module can be retrieved by running the bash script [getOptimalModule.sh](https://github.com/jchow32/MAGI-MS/blob/main/magi-ms/getOptimalModule.sh): 
	`./getOptimalModule.sh $output_file_name`
* In practice, <i>Cluster2</i> can be run multiple times on the same `BestPaths*` generated from a single run of <i>Pathway_GeneCenter</i>, and the highest scoring module can be retrieved from all runs of <i>Cluster2</i>.
* 36 example input parameters (a single combination of parameters per line) for <i>Cluster2</i> are displayed in the [Inputs](https://github.com/jchow32/MAGI-MS/blob/main/magi-ms/Inputs) file, in which the parameters -l, -i, -avgCoExpr, and -avgDensity are varied. 36 example output names are also provided in the file [Outputs](https://github.com/jchow32/MAGI-MS/blob/main/magi-ms/Outputs2), which correspond to the suggested input parameters listed in [Inputs](https://github.com/jchow32/MAGI-MS/blob/main/magi-ms/Inputs). Note that the first line of [Inputs](https://github.com/jchow32/MAGI-MS/blob/main/magi-ms/Inputs) and [Outputs](https://github.com/jchow32/MAGI-MS/blob/main/magi-ms/Outputs2) are substituted into the example <i>Cluster2</i> command shown above. 

