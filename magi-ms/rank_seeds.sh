score_list=$1
gene_list=$2
output_name=$3

join -t$'\t' -1 1 -2 1 <(sed 's/ /\t/g' $score_list | sort -k1,1) <(sort $gene_list) | sort -k2,2rg | cut -f1,2 > $output_name
