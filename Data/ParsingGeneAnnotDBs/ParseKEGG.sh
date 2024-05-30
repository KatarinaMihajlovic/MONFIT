###run the script as bash name.sh
#!/usr/bin/env bash

# Get list of human KEGG pathways
wget --no-check-certificate -O input/KEGG_PathwaysList.txt http://rest.kegg.jp/list/pathway/hsa
# Get list of links between human genes and KEGG pathways
wget --no-check-certificate -O input/KEGG_Pathway2Gene.txt http://rest.kegg.jp/link/hsa/pathway/ 
