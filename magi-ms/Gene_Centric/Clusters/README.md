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
