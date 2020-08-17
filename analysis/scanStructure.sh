
#!/bin/bash

#build structure database
for READ in `ls stockholm`;
do
    echo "$READ"
      cmbuild 'cm/'$READ'.cm' 'stockholm/'$READ''
      cmcalibrate --forecast 'cm/'$READ'.cm'
      cmcalibrate 'cm/'$READ'.cm'
done

cat cm/*cm > database/Template42.cm

cmpress database/Template42.cm
#If a cm database already exists, please mask the above code.

#search structure
cmscan -o Sym_ITS2_core.cm --textw 2000 database/Template42.cm "Sym_ITS2_core.fasta" 
