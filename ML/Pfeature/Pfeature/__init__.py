# -*- coding: utf-8 -*-
"""
Pfeature
========

        Pfeature is a comprehensive package which will allow users to compute most of the protein features that have been discovered over the past decades. Using different functional modules of this package, users will be able to evaluate six major categories of protein features: i) Composition-based features, ii) Binary profile of sequences, iii) evolutionary information based features, iv) structural descriptors, v) pattern based descriptors, and vi) model building, for a group of protein/peptide sequences. Additionally, users will also be able to generate these features for sub-parts of protein/peptide sequences. This will be helpful to annotate structure, function and therapeutic properties of proteins.


Important:
==========
1: To know how to import functions, please have a look at file "Functions_Tables.pdf".
2: To know more about the headers of the result files, please have a look at file "Pfeature_descriptors.pdf"

**Boths files are available in the zipped folder.**

Using
-----

        Just write in Python

        >>> from Pfeature.pfeature import *      # To import the all function together
        >>> from Pfeature.pfeature import aac_wp #To import the function to calculate amino acid composition for whole protein/peptide sequences.

Abbreviation:
###########Whole sequence and Sub-sequences################
wp    : Whole Protein
nt    : N-terminal residues
ct    : C-terminal residues
rt    : Rest of the residues after removing N- and C-terminal residues
nct   : N- and C-terminal residues together
st    : Number of splits generated from the input sequence(s).
##########Composition#####################################
aac   : Amino acid composition
dpc   : Dipeptide composition
tpc   : Tripeptide composition
atc   : Atomic composition
btc   : Bond composition
pcp   : Physico-chemical properties composition
aai   : Amino-acid indices composition
rri   : Residue repeat information
pri   : Physico-chemical properties repeat information
ddr   : Distance distribution of residues
sep   : Shannon entropy of protein
ser   : Shannon entropy of residues
spc   : Shannon entropy of physico-chemical properties
acr   : Autocorrelation descriptors
ctc   : Conjoint triad descriptors
ctd   : Composition enhanced transition distribution
paac  : Pseudo amino acid composition
apaac : Amphiphilic pseudo amino acid composition
qso   : Quasi sequence order
soc   : Sequence order coupling number
##########Binary#####################################
aab   : Amino acid based binary profile
dpb   : Dipeptide based binary profile
atb   : Atom based binary profile
btb   : Bond based binary profile
pcb   : Physico-chemical properties based binary profile
aib   : Amino-acid indices based binary profile
##########PSSM#####################################
pssm_comp: Composition of the binary profiles
pssm_n1  : It normalizes pssm profile based on 1/(1+e^-x) formula
pssm_n2  : It normalizes pssm profile based on (x-min)/(max-min) formula
pssm_n3  : It normalizes pssm profile based on ((x-min)/(max-min))*100 formula
pssm_n4  : It normalizes pssm profile based on 1/(1+e^-(x/100) formula
##########Pattern#####################################
pat_bin  : It generates patterns of binary profiles for amino acid sequences with given window length
pat_str  : It generates patterns of strings from amino acid sequences with given window length
pat_csv  : It generates patterns from submitted csv file with given window length
pat_pcp  : It calculates physico-chemical properties composition of patterns generated from amino acid sequences with given window length
pat_aai  : It calculates amino acid indices composition of patterns generated from amino acid sequences with given window length
"""
