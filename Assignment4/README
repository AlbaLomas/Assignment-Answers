Assignment: use of BLAST to Discover putative Orthologues 
Orthologue DEFINITION:  Genes that are found in different species that evolved from a common ancestral gene by speciation.
A common first-step in discovery of Orthologues is to do a “reciprocal best BLAST”. That is, you take protein X in Species A, and BLAST it against all proteins in Species B.  The top (significant!!) hit in Species B is then BLASTed against all proteins in Species A.  If it’s top (significant) hit is the Protein X, then those two proteins are considered to be Orthologue candidates.  (there is more work to do after this, but this is a good start)
Using BioRuby to blast and parse the blast reports, find the orthologue pairs between species Arabidopsis and S. pombe.  I have uploaded their complete proteomes to the Moodle for you.  You do not need to create objects for this task (the existing BioRuby objects are sufficient)
To decide on "sensible" BLAST parameters, do a bit of online reading - when you have decided what parameters to use, please cite the paper or website that provided the information.
The BLAST parameters that I implemented in the script for the best detection of orthologs as reciprocal best hits are:
-	the combination of soft filtering with a Smith–Waterman final alignment (-F ‘‘m S’’ -s T). These options resulted in both the highest number of orthologs and the minimal error rates. 
-	An e-value threshold of 1*10^-6, as it was common in different papers of orthologues searching (The lower e-value, the more significant is the alignment). 
-	A query coverage of at least 50%, some papers also use a coverage of 60% but as I already implemented the harder filtering parameters and there were other papers that use 50% as coverage threshold, I decided to use 50%. 



Gabriel Moreno-Hagelsieb, Kristen Latimer, Choosing BLAST options for better detection of orthologs as reciprocal best hits, Bioinformatics, Volume 24, Issue 3, 1 February 2008, Pages 319–324, https://doi.org/10.1093/bioinformatics/btm585
Hernández-Salmerón, J.E., Moreno-Hagelsieb, G. Progress in quickly finding orthologs as reciprocal best hits: comparing blast, last, diamond and MMseqs2. BMC Genomics 21, 741 (2020). https://doi.org/10.1186/s12864-020-07132-6
Wall DP, Fraser HB, Hirsh AE. Detecting putative orthologs. Bioinformatics. 2003 Sep 1;19(13):1710-1. doi: 10.1093/bioinformatics/btg213. PMID: 15593400.

Pearson WR. An introduction to sequence similarity ("homology") searching. Curr Protoc Bioinformatics. 2013 Jun;Chapter 3:Unit3.1. doi: 10.1002/0471250953.bi0301s42. PMID: 23749753; PMCID: PMC3820096.


Bonus:  1%
Reciprocal-best-BLAST is only the first step in demonstrating that two genes are orthologous.  Write a few sentences describing how you would continue to analyze the putative orthologues you just discovered, to prove that they really are orthologues. You DO NOT need to write the code - just describe in words what that code would do.
“The most conclusive evidence that two similar genes are orthologous is the result of phylogenetic analysis of the lineage of that gene. Genes found within the same clade are orthologous because they are descended from the same ancestor.”
The analysis can be continued by obtaining phylogenetic trees of genes that are thought to be orthologues of the target species. I think that I would be necessary to include in these phylogenetic trees more species that have the target genes so it would be possible to know where the “mother” gene of each pair of orthologues was divided in the species and how was the evolution. 

Focusing on Blast results, it is also possible to change the parameters in a harder way to found a smaller and more accurated number of ortholog genes. 
