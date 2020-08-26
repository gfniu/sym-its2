## Required software
        CD-HIT v4.8.1  
        GraphClust v0.7.6  
        HMMER v3.3  
        Hapstar v0.7  
        Infernal v1.1.2
        MAFFT v7.427
        MrBayes v3.0b4
        MrModeltest v2.4
        PAUP* v4.0a165        
        RAxML v8.2.10
        R2R v1.0.5
        ViennaRNA v2.4.11
        WebLogo v2.8.2

## Bioinformatic pipeline:
+ ### GenBank
1. *Download sequences*  
The script **searchGB.py** downloads Symbiodiniaceae rDNA sequence containing ITS2 from GenBank using keyword search command. The downloaded sequences are written to "~/data/1_GenBank/gb_its2.fasta".
2. *Multiple sequence alignment*  
The gb_its2 sequences having 5.8S or 28S annotations are aligned and written to "\~/data/1_GenBank/gb_its2_anno_aln.fasta".
3. *Extract boundary information*  
The script **extractBounds.py** extracts boundary information based on GenBank annotation, which is used to map the rDNA operon elements onto the gb_its2 alignment. This step generates the stacked bar chart ("\~/data/1_GenBank/gb_its2_anno_stack.xlsx") as shown in Fig. 1 of the manuscript.
* ### HMM
1. The gb_its2 sequences with both 5.8S and 28S annotations ("\~/data/2_HMM/gb_its2_5.8S_28S_anno.gb") are scanned for the ITS2-proximal stem using the ITS2 database web tool (http://its2.bioapps.biozentrum.uni-wuerzburg.de).  
2. The file "\~/tables/Table S1.htm" lists ITS2-proximal stems identified from 1,328 gb_its2 sequences that are used to create the logo in Fig. 2 of the manuscript.  
3. A Symbiodiniaceae-specific profile Hidden Markov Model of the ITS2-proximal stem is generated and saved in "\~/data/2_HMM/hmmdb/".
* ### Core references
1. *Merge reference datasets*  
The Arif et al. 2014 (\~/data/3_Core_references/Arif_derep.fasta) and Cunning et al. 2017 (\~/data/3_Core_references/Cunning_derep.fasta) reference sequences are combined and dereplicated. The result is written to "\~/data/3_Core_references /ref_its2.fasta".
2. *ITS2 delineation*  
The script **scanStem.py** delineates ITS2 sequences from the "\~/data/3_Core_references/ref_its2.fasta" and "\~/data/1_GenBank/gb_its2.fasta" files by scanning ITS2-proximal 5.8S/28S stem with the Symbiodiniaceae-specific HMM.
3. *ITS2 clustering*  
The script **clusterSeq.sh** clusters the delineated ref_its2 sequences using cd-hit-4.8.1 to generate the core reference dataset ("\~/data/3_Core_references/Sym_ITS2_core.fasta"). The script also collapses some gb_its2 sequences into the reference core (see details in Table S2 of the manuscript). GenBank sequences that do not collapse into the reference core are saved as "\~/data/3_Core_references /Sym_ITS2_extra.fasta". Table S3 of the manuscript contains detailed delineation and clustering information for this extra dataset.
* ### Haplotype
1. *Multiple sequence alignment*  
The script **alignSeq.sh** aligns each clade of Symbiodiniaceae ITS2 sequences in the core dataset ("\~/data/4_Haplotype/sequences/*.fasta") using MAFFT, and writes the results to "\~/data/4_Haplotype/alignment/".
2. *Minimum spanning tree*  
The script **mst.R** generates a matrix of absolute pairwise nucleotide difference ("\~/data/4_Haplotype/MST/\*.matrix.csv") with the haplotypes package v1.1 in R, which is subsequently converted to a minimum spanning tree (MST) ("\~/data/4_Haplotype/MST/\*.txt/").
3. *Haplotype network*  
The MST is processed and visualized with hapstar. The generated haplotype networks for clade C and other clades correspond to Fig. 3 and Fig. S2 of the manuscript, respectively.
* ### Phylogeny
1. *Multiple sequence alignment*  
The script **alignSeq.sh** aligns the extended set of 1,322 distinct Symbiodiniaceae ITS2 variants ("\~/data/5_Phylogeny/Sym_ITS2_extended.fasta") using MAFFT. The alignment is saved as "\~/data/5_Phylogeny/Sym_ITS2_extended_aln.fasta".
2. *Phylogenetic tree reconstruction*  
The script **buildTree.sh** reconstructs maximum likelihood trees, inferred from Bayesian posterior probabilities and bootstrap support with MrBayes and RAxML, respectively, using the best-fit model of nucleotide substitution estimated by MrModeltest. The script also reconstructs phylogenetic trees under maximum parsimony criterion with PAUP* 4.0a165, but using within-clade subsets of sequences.
* ### Genetic distance
1. *Multiple sequence alignment*  
The script **alignSeq.sh** aligns 583 ITS2 sequences ("\~/data/6_Genetic_distance/Sym_ITS2_core_583.fasta") in the reference core that span all the nine clades using MAFFT.
2. *Calculation of genetic distances*  
The script **calDist.nex** takes the alignment generated above ("\~/data/6_Genetic_distance/Sym_ITS2_core_583_aln.fasta") as input and calculates GTR model-corrected pair-wise within-clade genetic distances using PAUP* 4.0a165. The result is written to “\~/data/6_Genetic_distance/Sym_ITS2_core_583_distance.csv”.
* ### Secondary structure
1. *Build and scan CM database for similar structures*  
The script **scanStructure.sh** takes 42 published secondary structures as template ("\~/data/7_Secondary_structure/template/") to build a calibrated covariance model (CM) database ("\~/data/7_Secondary_structure/database/Template42.cm.*") using Infernal 1.1. ITS2 sequences with unknown structures from the core set ("\~/data/7_Secondary_structure/Sym_ITS2_core.fasta") are scanned against this CM database for regions with similarly base-paired helix structures.
2. *Format Infernal output*   
The script **cm2xfasta.py** converts the Infernal output ("\~/data/7_Secondary_structure/Sym_ITS2_core.cm") to XFASTA format.
3. *Calculate free energy*  
The script **calEnergy.py** calculates free energy for each secondary structure using ViennaRNA-2.4.9.
4. *View the structures*  
The structures are plotted and visualized in R2R-1.0.5. The predicted secondary structures of 619 core ITS2 sequences, in XFASTA format as well as individual images, can be found in the compressed file "\~/files/File S1.zip" accompanying the manuscript.
* ### Consensus structure
1. *Cluster secondary structures*    
For each clade of Symbiodiniaceae ITS2 sequences in the core dataset ("\~/data/3_Core_references /Sym_ITS2_core.fasta"), clustering of their secondary structures is performed with graphclust-0.7.6, generating a consensus sequence per each clade ("\~/data/8_Consensus_structure/con_struct_seqs.fasta")
2. *Predict clade-specific secondary structures*  
The consensus sequences from step 1 is used as input for Infernal-1.1.2 to search against the template database ("\~/data/7_Secondary_structure/database/Template42.cm.*") for similar secondary structures, which are then used as templates to generate structural alignments for the sequences clustered into each clade.
3. *Depict consensus secondary structures*  
The script **drawConStruct.sh** takes structural alignments generated in step 2 as input to draw consensus secondary structures using R2R-1.0.5. The clade-specific multi-sequence structural alignments in STOCKHOLM format and the corresponding consensus secondary structures deduced from 565 core sequences, as shown in Fig. 5 of the manuscript, are saved in the compressed file "\~/files/File S2.zip" accompanying the manuscript.
