from Bio.Blast import NCBIXML

file = open("BIOPYTHON/result.xml")
records = NCBIXML.parse(file)
for record in records : #records k andar individual seq objects as record
   for alignment in record.alignments : #record.alignments ek attribute hai records ka
      #jo ki ek list ki form mai hai ab jiske andar alignment objs hai
     for hsp in alignment.hsps :
        print(hsp)

"""records (iterator)
    ↓
Record object
    ↓
alignments (list)
    ↓
Alignment object
    ↓
hsps (list)
    ↓
HSP object
    ↓
expect, score, bits..."""