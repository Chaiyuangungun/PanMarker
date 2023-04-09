# Usage

        ls *.prm >prm.file 
        ls *.cds >cds.file

1、Use cds

        python3 PanMarker.py -i cds.file -p trait.file -e FPKM.file -s cds -o prefix -g T\F【-t num -a person_cor】

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

2、Use prm

        python3 PanMarker.py -i prm.file -p FPKMsfile -s prm -o prefix -g T\F【-t num】

3、result

        prefix.result(all variant sites associate with traits)
        
        prefix.out(top 5% variant sites associate with traits)
