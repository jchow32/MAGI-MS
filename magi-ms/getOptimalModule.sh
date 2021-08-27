output_file=$1
highScore=$(cut -d' ' -f7 $output_file | awk '{if ($0 ~ /\./) {print $0}}' | sort -k1,1r | grep -v A | head -n1)
highScore_line=$(grep -n $highScore $output_file | cut -d: -f1 | head -n1) 
moduleSize=$(head -n $highScore_line $output_file | grep '^[0-9]' | tail -n2 | head -n1)
start_index="$(($highScore_line-$moduleSize))"
sed -n "$start_index,$highScore_line p" $output_file | head -n $moduleSize > Optimal_module.txt
