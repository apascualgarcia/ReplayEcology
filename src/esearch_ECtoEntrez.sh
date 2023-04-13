#/bin/bash
# Requires: Entrez direct installed: https://anaconda.org/bioconda/entrez-direct
# Input: the tab-separated table of metagenomic predictions with EC numbers in the first column without "EC", (e.g. 1.1.1.105)
# Output: a list of EC and Gene.
tableIn="pred_metagenome_unstrat_2BioCyc.tsv"
listOut="EC_codes_wGene.list"

##### STOP EDITING
echo "Extract EC numbers"
tableTmp="EC_codes.list"
length=$(wc -l $tableIn | cut -f 1 -d ' ')
length=$((length-1))
#length=2 # debug
cut -f 1  $tableIn | tail -n $length > $tableTmp

#ECS from the table, store into an array
echo "Loading the codes in an array"

arr=()
while IFS= read -r line; do
   arr+=("$line")
done < $tableTmp
# alternative
#arr=$(head -n 2 $tableTmp) # | xargs);

echo "Start searching, check file $listOut for progression"
i=1
for ECS in ${arr[@]}
#head -n 2 $tableTmp | while IFS="" read -r ECS || [ -n "$ECS" ]
do
i=$((i+1)) # debug
# this line do work --> Change to Prot-ref_desc the field to retrieve
GENE=$(esearch -db gene -query ${ECS} | efetch -format xml | xtract -pattern Entrezgene -unless OrgName_lineage -contains Bacteria -element Gene-ref_locus | tr '[:upper:]' '[:lower:]' | awk '{count[$1]++} END {for (word in count) print word, count[word]}' | sort -gr -k 2,2 | head -n 1 | cut -f 1 -d ' ' )


echo -e ${ECS}'\t'${GENE}
sleep 0.5
done > $listOut
echo "Done!"
rm -f $tableTmp

# Debug: this line do the search but returns only one line
#GENE=$(esearch -db gene -query "1.1.1.10" | efetch -format xml | xtract -pattern Entrezgene -unless OrgName_lineage -contains Bacteria -element Gene-ref_locus | tr '[:upper:]' '[:lower:]' | awk '{count[$1]++} END {for (word in count) print word, count[word]}' | sort -gr -k 2,2 | head -n 1 | cut -f 1 -d ' ' )

# Debug: This line simply takes the next element of the table to print it with the current one. It works fine.
#GENE=$(head -n $i $tableTmp | tail -n 1) # debug


#do echo "$(esearch -db gene -query ${ECS} | efetch -format xml | xtract -pattern Entrezgene -unless OrgName_lineage -contains Bacteria -element Gene-ref_locus)"
#do echo "$(esearch -db gene -query ${ECS} | efetch -format xml | xtract -pattern Entrezgene -unless OrgName_lineage -contains Bacteria -element Gene-ref_locus | tr '[:upper:]' '[:lower:]' | awk '{count[$1]++} END {for (word in count) print word, count[word]}' | sort -gr -k 2,2 | head -n 1 | cut -f 1 -d ' ')"

##ECS2="${arr[@]::$i}" # debug
