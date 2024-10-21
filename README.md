# PRF-Search
Scripts to search human protein coding transcript sequences for motifs that can induce -1 ribosomal frameshifts via transmembrane segment insertion via Sec61 translocon
--------------------------------------------------------------------------------
Prerequisites:

  - python=3.8
  - ipython
  - matplotlib
  - numpy
  - openpyxl
  - pandas
  - scikit-learn
  - scipy
  - seaborn
  - statsmodels
  - xlrd
  - viennarna
  - biopython
  - blast

Install a conda environment with the included yml file:

conda env create -f prf_search_env.yml

This will work in Fedora 40 and probably other RHEL-like operating systems. If you are using other operating systems you can try the platform-independent yml file:

conda env create -f prf_search_env_from-history.yml
--------------------------------------------------------------------------------

In order to adjust TMD coordinates using the von Heijne transfer free energy scale use 
von_heijne_scan_functions_v39_signal_pep_topcons.py

This can use TMD predictions from TMHMM, DeepTMHMM, and TOPCONS2

We used TOPCONS2 predictions produced on human protein coding sequences from Ensembl translated to amino acid residues, using the included script:
translate_cds_v2.py

TOPCONS2 produces a text file of outputs that von_heijne_scan_functions_v39_signal_pep_topcons.py can parse. To run the tool that adjusts TMD coordinates use this command:

python von_heijne_scan_functions_v39_signal_pep_topcons.py human_GRCh38_new_cds_translated.fasta topcons_Homo_sapiens_GRCh38_new_pep.txt topcons -lmin 16 -lmax 25 -o topcons_ensembl_output.csv -d dg_topcons_ensembl_output.txt

Here, human_GRCh38_new_cds_translated.fasta is the fasta containing translated Ensembl sequences and topcons_Homo_sapiens_GRCh38_new_pep.txt is the output of TOPCONS2. topcons_ensembl_output.csv is the output that gives the original TMD coordinates and the adjusted TMD coordinates, as well as the transfer free energy of each TMD (both unadjusted and adjusted).
--------------------------------------------------------------------------------
Score the propensity of a heptamer to facilitate -1 ribosomal frameshifts using frameshift_routines_v10.py
This requires the module turner_energies.py in the working folder. To score a list of heptamers choose option 3 when running
python frameshift_routines_v10.py

Upload a text file containing a list of heptamers. The script will output a CSV containing each heptamer and a score that calculates the difference in codon-anticodon binding in the PTC between the zero reading frame and the -1 reading frame.

Scores of zero indicate no difference in the energy of codon-anticodon binding between the zero and -1 reading frames.
Higher scores indicate a greater unfavorability. For our collection of slip site candidates we chose a cutoff of +2.89 kcal/mol. This is sufficient to capture all 24 strictly-defined heptamers established to promote -1 programmed ribosomal frameshifting in viruses. I have included all_heptads.txt and 24_strict_slipsites.txt as lists of heptamers for users to test with this script. They will produce the results all_heptamers_friction_scan.csv and noncanonical_slipsites_2.89_cutoff_friction_scores.csv when used with this script.
--------------------------------------------------------------------------------
Identify Harrington motifs and run a few statistical tests on them using enhanced_harrington_search_v60_production_version_final.py

This requires a fasta file, an Excel file downloaded from Uniprot that contains data about specific Uniprot sequences, a blastp output TSV that searches a blastp database built from Uniprot sequences for Ensembl protein coding sequences translated to protein sequences (this searches for associations between Ensembl coding sequences and Uniprot sequences), a path to TOPCONS2 predictions, a list of slip site candidates, and a CSV containing the slip site candidates and their associated friction scores:

path_fasta = Homo_sapiens.GRCh38_new.cds.all.fa
path_uniprot = uniprotkb_taxonomy_id_9606_AND_reviewed_2024_09_07.xlsx
blastp_results_path = blastp_results_uniprot_all_hs_rev_ensembl_all.tsv
path_tmd_csv = topcons_ensembl_output.csv
path_slipsites = noncanonical_slipsites_2.89_cutoff.txt
path_friction_csv = noncanonical_slipsites_2.89_cutoff_friction_scores.csv

These files will be made available to the community upon request.

We have included a config file config_file_v60.ini that contains file paths and parameters for the Harrington motif search.
To run the pipeline:
python enhanced_harrington_search_v60_production_version_final.py config_file_v60.ini

The script will produce a background distribution of distances between randomly selected heptamers and upstream TMDs. This background is quite large and is saved as .pkl files but can be sampled using pickle_to_csv_sampled_v4.py

This will produce a CSV containing the same number of distances as was present in the list of distances produced for slippery heptamers to TMDs.





