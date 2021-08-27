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
