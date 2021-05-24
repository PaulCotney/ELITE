# ELITE

Pre-requisite for running ELITE:

1. Install msBWT from here: https://github.com/holtjma/msbwt
2. Install bowtie from here: http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
3. Install Levenshtein python package from here: https://pypi.org/project/python-Levenshtein/
4. Build msBWT from short-read  (a sample msbwt is given in the sample_bwt folder)
    a. Follow instructions from here https://github.com/holtjma/msbwt

Prepare input data:

1. Create a .csv file with the following information for each sample: (a sample file sample_list.csv is given in the sample_data folder)
    a. sample name: name of you sample
    b. msBWT directory: directory where the sample's msBWT is located
    c. threshold: minimum number of supporting split-reads with TE
    
2. Create a .csv file containing the TE templates: (a sample file TEseq.csv is given in the sample_data folder)
    a. my_id: name of the TE template
    b. TE: sequence of the TE
    c. start_seed: seed kmer (lengths between 21-35) near the TE's proximal boundary [empty string if you want ELITE to find the optimal seed]
    d. end_seed: seed kmer (lengths between 21-35) near the TE's distal boundary [empty string if you want ELITE to find the optimal seed]
    
3. Reference genome: (a sample genome Reference.fa is given in root folder)
    a. The reference genome of the sample species

Required parameters for running ELITE:
    a. sample_list: link to the .csv file containing the sample's information
    b. te_file: link to the .csv file containing the TE templates information
    c. reference: link to the reference genome file
    d. species: name of your samples' species
    e. C: length of TE's context which is the segment of non-TE sequence right next to a TE sequence [C = 25 recommended]
    f. T: length of proximal and distal TE length which is the segment of TE sequence right next to a context [C + T must be less than half of read length]
    g. K: length of seed [25 recommended]
    
a sample script runELITE.sh is given to run ELITE using the sample data. 
It will automatically create bowtie index if not already built. 
But if you have it built, you must use the same reference genome that will be used here.
