
### 191023 Martin Raden

- copraRNA_cleanupWD.sh
  - dedicated script to cleanup a CopraRNA job's working directory
- CopraRNA2.pl
  - cleanup removed with calling copraRNA_cleanupWD.sh
- adding 'suppressPackageStartupMessages' to all lib loadings in R scripts
- homology_intaRNA.pl
  - DAVID enrichment visualization done only if enrichment data present
- CopraRNA2-deps.yml
  - r-base >= 3.6.0 (was ==)
- obsolete 'edit' comments removed
