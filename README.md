# PanMarker
ls *.prm >prm.file 
ls *.cds >cds.file
use cds file
python3 PanMarker.py -i cds.file -p trait.file -e FPKM.file -s cds [-t num ]
