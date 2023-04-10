# Dependencies

software

        MAFFT

Python Modules
        
        outlier-utils
# Usage

        ls *.prm >prm.file 
        ls *.cds >cds.file

1、Extract variants and extart variant sites associate with traits

use cds file :

        python3 PanMarker.py -i cds.file -p trait.file -e FPKM.file -s cds -o prefix -g T\F【-t num -a person_cor】
                -i INPUTFILE, --inputfile  input file(cds or prm file list)
                -p PHE, --phe              phenotype file
                -e EXP, --exp              expression profile
                -s TYPE --type             file type, cds or prm
                -o OUTPUT, --output        output file prefix
                -g GRU, --gru              T or F,Ture or False, whether to perform phenotype outlier filtering
                -t THREAT, --threat        Number of threads (default=10）
                -a PERVALUE --pervalue     Pearson correlation coefficient(default=0.3)

trait.file :

        sample1 value1

        sample2 value2

        ...

        sampleN valueN

FPKM.file :

                gene1 gene2 gene3 ... geneN
       
        sample1 xx  xx  xx  ... xx

        sample2 xx  xx  xx  ... xx

        ...

        sampleN xx  xx  xx  ... xx

use prm file:
        
        python3 PanMarker.py -i prm.file -p FPKM.file -s prm -o prefix -g T\F【-t num】
                -i INPUTFILE, --inputfile  input file(cds or prm file list)
                -p PHE, --phe              expression profile
                -s TYPE --type             file type, cds or prm
                -o OUTPUT, --output        output file prefix
                -g GRU, --gru              T or F,Ture or False, whether to perform phenotype outlier filtering
                -t THREAT, --threat        Number of threads (int,default=10）
                -a PERVALUE --pervalue     Pearson correlation coefficient(float,default=0.3)
        
2、result

        prefix.result(all variant sites associate with phenotype)
                
        prefix.out(top 5% variant sites associate with phenotype) 
 
                
