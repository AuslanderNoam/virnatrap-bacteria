viRNAtrap-bacteria extands virnateap to generate predicted viral and bacterial contigs from unmapped RNAseq reads

## Python scripts, in alphabetical order:
### analysis_proteins_esca_vs_eso_fix.py
	Purpose: Perform analyses of bacterial protein representation in TCGA and GTEx. Generated data for Fig. 3A nd plot for Fig. 3B.
	Dependencies: numpy, pandas, seaborn, scipy, lifelines, Matplolib, statsmodels
	Arguments: Hardcoded - uses a database (FASTA) of bacterial proteins, binary matrices of protein presence in TCGA and GTEx samples, and TCGA clinical data
### analysis_proteins_survival.py
	Purpose: Perform survival analyses for bacterial proteins, human gene expression analyses, and survival analyses of human proteins. Generated Figs. 3C and 4.
	Dependencies: same as above
	Arguments: same as above, plus expression of human genes in TCGA samples
### cnn3class.py
	Purpose: trains a 3-class CNN to distinguish human, viral, and bacterial reas
	Dependencies: TensorFlow/Keras, numpy
	Arguments: training and validation sequences, model hyperparameters. Run python cnn3class.py -h for details.
### encode_seqs.py
	Purpose: Preprocess and compress sequences for use in model training. Training and validation sequences used for model training should be pre-processed by this script. Shuffling is optional but hardcoded. Any order of classes is permitted, but we used Virus=0, Human=1, Bacteria=2 throughout.
	Dependencies: numpy, BioPython
	Arguments: an output file, followed by one or more FASTA files. Fasta files are sequentially assigned classes starting with 0. 
### esoph_bact.py
	Purpose: Perform bacterial genera abundance analysis & binomial tests. Also suitable for viral species analysis with modification to process entire species names (see inline comment).
	Dependencies: scipy
	Arguments: Four mappings of contigs to species, GTEx sample metadata, abundance threshold (float). Mappings are the output of run_blast_os_par_oldv.py; first three are assumed to be for non-overlapping TCGA samples and last for all GTEx samples.
### filter_reads.py
	Purpose: Removes sequences containing non-ACTG characters froma FASTA file.
	Arguments: A FASTA file. Writes to STDOUT. 
### fragment_split.py
	Purpose: Divides a set of sequences into training, validation, and test sets, and then splits each sequence into simulated reads of a fixed length. Sequences will be shuffled, then assigned to the training set until the target size is exceeded, then assigned to the validation set until the target size is exceeded, then assigned to the test set. 
	Dependencies: BioPython
	Arguments: A FASTA file, a target number of training reads, a target number of validation reads, a target "read" length, and a target "step" or offset between consecutive reads.  
### perf_plot.py
	Purpose: Generate PRC and ROC plots, used in Fig. 1B-C
	Dependencies: numpy, TensorFlow/Keras, scikit-learn, cnn3class.py, Matplotlib
	Arguments: A TF/Keras model to load, a file of encoded sequences to use, and a file name for the plot
### postproc_blastx.py
	Purpose: Post-process output of run_blastx.py into a matrix
	Dependencies: pandas, numpy, BioPython, multiprocessing, ast
	Arguments: Run python postproc_blastx.py -h for details. Requires the output of run_blastx.py. Optional project name.
### run_blast_os_par_oldv.py
	Purpose: Run blastn to map contigs to species from 5 databases: bacteria, phages, HERVs, reference viruses, and novel viruses. 
	Dependencies: pandas, numpy, BioPython, subprocess, multiprocessing. Also requires blastn available on path and $BLASTDB defined, with appropriate databases present.
	Arguments: Run python run_blastx.py -h for details. Input directory (output directory of virnatrap.py) and output directory. Number of threads is hardcoded. TCGA (True/False) is hardcoded; if True, tries to extract additional sample information. 
### run_blastx.py
	Purpose: Run blastx to assign contigs to bacterial proteins.
	Dependencies: pandas, numpy, BioPython, subprocess, multiprocessing. Also requires blastx available on path and $BLASTDB defined, with appropriate database present.
	Arguments: Run python run_blastx.py -h for details. Requires input directory (output directory of virnatrap.py). Optional project name. Number of threads is hardcoded.
### scatter.py
	Purpose: Generate Fig. 1B (without labels)
	Dependencies: Matplotlib
	Arguments: CSV of TCGA and ESCA abundances (output of esoph_bact.py), list of cancer-abundant species/genera, list of healthy-abundant species (line-seperated)
### test_n.py
	Purpose: Test the performace of a model on sequences with 1-4 artificial mutations.
	Dependencies: numpy, TensorFlow/Keras, scikit-learn, cnn3class.py
	Arguments: A TF/Keras model to load and a file of encoded sequences to use
### test_plot.py
	Purpose: Display a 3D scatter plot of scores of random sequences, used in Fig. 1D
	Dependencies: numpy, TensorFlow/Keras, scikit-learn, Matplotlib
	Arguments: A TF/Keras model to load and a file of encoded sequences to use. 
### test_seed.py
	Purpose: Test the behavior of different "seed read"-selection thresholds. Used to display Fig. S1.
	Dependencies: numpy, TensorFlow/Keras, scikit-learn, cnn3class.py, Matplotlib
	Arguments: A TF/Keras model to load and a file of encoded sequences to use. 
### virnatrap.py: 
	Purpose: Run the contig-generation pipeline. Processes all FASTQs in a directory and writes contigs for each to a seperate file in an output directory, except that FASTQs with an existing output file in that directory are skipped. 
	Dependencies: multiprocessing, numpy, ctypes, TensorFlow/Keras. Requires src/assemble_read_c.c to be compiled into src/assemble_read_c.so. 
	Arguments: Run python virnatrap.py -h for details. Requires an input directory of FASTQs, an output directory, and a model to load. Optional (with defaults) are an expected read length and a number of threads to use.
Common, built-in dependencies (e.g., sys) are omitted.

R scripts:
plot_fig3.R: Used to generate Fig. 3A
