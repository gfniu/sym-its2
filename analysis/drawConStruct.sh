
#!/bin/bash

SAVEIFS=$IFS
IFS=$(echo -en "\n\b")

for READ in `ls input`;
do
     if [[ $READ == *.sto ]];
     then
           echo "$READ"
           GSCfilename='GSC'$READ''
           r2r --GSC-weighted-consensus 'input/'$READ'' 'temp/'$GSCfilename'' 3 0.97 0.9 0.75 4 0.97 0.9 0.75 0.5 0.1
           echo 'temp/'$GSCfilename'' > 'temp/'$READ'.r2r_meta'
           r2r 'temp/'$READ'.r2r_meta' 'output/'$READ'.pdf'
     fi
done
IFS=$SAVEIFS
