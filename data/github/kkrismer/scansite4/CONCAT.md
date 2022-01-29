<p>
  <img src="Scansite4/CoreApplication/src/main/webapp/img/logo_scansite.png" alt="Scansite 4" style="height: 103px; width: 100%" />
</p>

[![license: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/licenses/MIT)

https://scansite4.mit.edu

Scansite searches for motifs within proteins that are likely to be phosphorylated by specific protein kinases or bind to domains such as SH2 domains, 14-3-3 domains or PDZ domains.

Putative protein phosphorylation sites can be further investigated by evaluating evolutionary conservation of the site sequence or subcellular colocalization of protein and kinase.

Scansite 4 was developed with GWT 2.8.2, Java 1.8, and MySQL 5.7.

## About

Scansite searches for motifs within proteins that are likely to be phosphorylated by specific protein kinases or bind to domains such as SH2 domains, 14-3-3 domains or PDZ domains.

Optimal phosphorylation sites for particular protein Ser/Thr kinases or protein-Tyr kinases are predicted using the matrix of selectivity values for amino acids at each position relative to the phosphorylation site as determined from the oriented peptide library technique described by Songyang et al., 1994, Current Biology 4, 973-982 and Songyang et al., 1995, Nature 373, 536-539.

Optimal binding sites for SH2 domains, PDZ domains, 14-3-3 domains and other domains are determined using the matrix of selectivity values for amino acids at each position relative to an orienting residue as determined by the oriented peptide library technique described in Songyang et al., 1993 Cell 72, 767-778, Songyang et al., 1997 Science 275, 73-77 and Yaffe et al., 1997 Cell 91, 961-971.

The same matrices of selectivity values are used in an approach to provide relative scores of candidate motifs in protein sequences evaluated. The Motifscanner program utilizes an entropy approach that assesses the probability of a site matching the motif using the selectivity values and sums the logs of the probability values for each amino acid in the candidate sequence. The program then indicates the percentile ranking of the candidate motif in respect to all potential motifs in proteins of a protein database. When available, percentile scores of some confirmed phosphorylation sites for the kinase of interests or confirmed binding sites of the domain of interest are provided for comparison with the scores of the candidate motifs.

In the graphical output, the candidate motifs are superimposed on the predicted domain structure of the protein. Clicking on the domain information button provides access to a hot link to the domain families via Pfam. Clicking on the motif region provides the primary sequence at the motif and the percentile score. The program also provides information about the surface probability of the region of the protein around the motif of interest. 

## References

[1] Obenauer JC, Cantley LC, Yaffe MB  
[Scansite 2.0: Proteome-wide prediction of cell signaling interactions using short sequence motifs](https://www.ncbi.nlm.nih.gov/pubmed/12824383)  
Nucleic Acids Res. 2003 Jul 1; 31(13): 3635-41.

[2] Yaffe MB, Leparc GG, Lai J, Obata T, Volinia S, Cantley LC  
[A motif-based profile scanning approach for genome-wide prediction of signaling pathways](https://www.ncbi.nlm.nih.gov/pubmed/11283593)  
Nat Biotechnol. 2001 Apr; 19(4): 348-53.

[3] Songyang Z, Blechner S, Hoagland N, Hoekstra MF, Piwnica-Worms H, Cantley LC  
[Use of an oriented peptide library to determine the optimal substrates of protein kinases](https://www.ncbi.nlm.nih.gov/pubmed/7874496)  
Curr Biol. 1994 Nov 1; 4(11): 973-82.

## Citation

If you use Scansite in your research, please cite:

Scansite 2.0: proteome-wide prediction of cell signaling interactions using short sequence motifs  
John C. Obenauer, Lewis C. Cantley, Michael B. Yaffe  
Nucleic Acids Res. 2003 July 1; 31(13): 3635–3641.  
PMID: [12824383](https://www.ncbi.nlm.nih.gov/pubmed/12824383)