# safb_paper

This repository accompanies the paper "Autonomous transposons tune their sequences to ensure somatic suppression" (Ilık et al., 202X). It containes the R code and data needed to annotate splice junctions reported by STAR (SJ.out.tab files) and analyse their usage. A subset of our data is provided as part of the repository.

## Clone the repository and uncompress raw data

    git clone https://github.com/aktas-lab/safb_paper.git`
    cd safb_paper
    tar -xvzf data/splice_junctions.tar.gz -C data/

## Running the analysis

    Rscript load_annotate_splice_junctions.R

Output table: `output/safb_annotated_junctions.txt`
