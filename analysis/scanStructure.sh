
#!/bin/bash

#build structure database
for READ in `ls stockholm`;
do
    echo "$READ"
      cmbuild 'cm/'$READ'.cm' 'stockholm/'$READ''
      cmcalibrate --forecast 'cm/'$READ'.cm'
      cmcalibrate 'cm/'$READ'.cm'
done
cd cm
cat .*cm > ../database/Template42.cm
cd ..
cmpress database/Template42.cm
#If a cm database already exists, please mask the above code.

#search structure
cmscan -o Sym_ITS2_core.cm database/Template42.cm "Sym_ITS2_core.fasta"
