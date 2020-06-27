
#!/bin/bash


#run raxml by intelmpi
mpirun -np 250 raxmlHPC-MPI-SSE3 -m GTRGAMMA -p 12345 -N 1000 -s alignment-1322.phy -n step1  -o out_JN558106,out_JN558107,out_JN558108,out_LC068842,out_LC068843,out_MG551867,out_LC413946,out_LC413949,out_AY999082,out_FJ024304

mpirun -np 250 raxmlHPC-MPI-SSE3 -m GTRGAMMA -p 12345 -b 12345 -N 1000 -s alignment-1322.phy -n step2  -o out_JN558106,out_JN558107,out_JN558108,out_LC068842,out_LC068843,out_MG551867,out_LC413946,out_LC413949,out_AY999082,out_FJ024304

raxmlHPC-SSE3 -m GTRGAMMA -p 12345 -f b -t RAxML_bestTree.step1 -z RAxML_bootstrap.step2 -n tre


# run mrbayes
BEGIN MRBAYES;
outgroup out_JN558106;
outgroup out_JN558107;
outgroup out_JN558108;
outgroup out_LC068843;
outgroup out_LC413949;
outgroup out_LC413946;
outgroup out_MG551867;
outgroup out_LC068842;
outgroup out_AY999082;
outgroup out_FJ024304;
Lset  nst=6  rates=gamma;
Prset statefreqpr=dirichlet(1,1,1,1);
mcmc Ngen=20000000 Samplefreq=1000 Printfreq=1000 savebrlens=yes nchains=16 Nswaps=5 Mcmcdiagn=yes Diagnfreq=1000 relburnin=yes burninfrac=0.25 checkfreq=1000;
sumt burnin=5000;
end;


#paup*
begin paup;
Execute A.nex;
set criterion=parsimony rootmethod=outgroup storebrlens=yes increase=no notifybeep=no outroot=monophyl storetreewts=yes;
pset gapmode=newstate;
Outgroup out_JN558106 out_JN558107 out_JN558108 out_LC068843 out_LC413949 out_LC413946 out_MG551867 out_LC068842;
bootstrap nreps = 200 search = heuristic keepall = yes grpfreq = yes treefile = bootstrap-A.tre brlens = yes / start = stepwise swap = tbr multrees = yes nreps = 200 addseq = random dstatus = 60 rstatus = yes;
savetrees file=A.tre BrLens=yes savebootp=nodelabels;
pscores all / TL = yes CI = yes RI = yes RC = yes HI = yes Scorefile = pscores-A.scores;
gettrees file = bootstrap-A.tre StoreTreeWts = yes StoreBrLens = yes mode = 3 duptrees = keep;
filter best = yes;
savetrees file=best-A.tre BrLens=yes savebootp=nodelabels;
gettrees file = bootstrap-A.tre StoreTreeWts = yes StoreBrLens = yes mode = 3 duptrees = eliminate;
contree all/strict=no Percent=50 MajRule=yes LE50=no Rootmethod=outgroup usetreewts=yes outroot=monophyl treefile=contree-A.tre;
end;
