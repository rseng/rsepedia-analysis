# ChiRA - Chimeric Read Analyzer

ChiRA is a set of tools to analyze RNA-RNA interactome experimental data such as CLASH, CLEAR-CLIP, PARIS, SPLASH etc. Following are the descriptions of the each tool. Here we provide descriptions about the input and ouptput files. For the detailed description of the other parameters please look at the help texts of tools.

## chira_collapse.py
* Description: This tool deduplicates the reads from the FASTQ file and writes into a fasta each read once with it's read count.
* Inputs
  * Quality and adapter trimmed FASTQ file
* Outputs
  * FASTA file with unique sequences. The headers of the sequence are in the following format: >sequence_id|UMI|read_count

## chira_map.py
* Description: This tool handles the mapping of the reads to reference transcriptome. User can choose between the bwa-mem and CLAN alignment tools.
* Inputs
  * A fasta file containing reads
  * A FASTA file containing reference transcript sequences
  * An optional second reference FASTA file, if split-reference is used
* Outputs
  * BED file containing the alignments
  * Optional unmapped FASTA file

## chira_merge.py
* Description: This tool merges the overlapping aligned positions to define the read concentrated loci. If an annotation GTF file produced, the transcriptomic alignment positions are first converted to their corresponding genomic positions.
* Inputs
  * Alignments in BED format
  * An annotation GFF file contaning reference genomic positions.
  * A FASTA file containing reference transcript sequences
  * An optional second reference FASTA file, if split-reference is used
* Outputs
  * BED file containing the alignments with reads categorized into segments depending on which part of the read is aligned.
  * BED file containing merged alignments. 4th column contains all the alignments merged into that location and 5th column contains number of reads.

## chira_quantify.py
* Description: This tool first creates CRLS from merged BED file and quantifies them based on the mapped reads.
* Inputs
  * A BED file containing alignment segments
  * A tabular  file containing merged alignments
* Outputs
  * Tabular file containing the reads and their CRLs with TPM values.

## chira_extract.py
* Description: This tool extracts the best chimeric alignments for each read. User can optionally hybridize the loci where the chimeric arms are mapping to. 
* Inputs
  * Tabular file containg CRLs information
  * Annotation GFF file
  * A FASTA file containing reference transcript sequences
  * An optional second reference FASTA file, if split-reference is used  
  * If your alignments are merged at genomic level in previous step (chira_merge.py), then provide a reference genomic FASTA fille.
* Outputs
  * Tabular file of chimeric reads information
  * Tabular file of detected interactions
