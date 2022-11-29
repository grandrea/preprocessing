# Preprocessing for Xi

Preprocessing script recalibrating MS1 and MS2 spectra based on mass error of a linear proteomics search. It uses Xi to perform the linear search. First step in the xiSEARCH workflow. It first converts Thermo .raw files into peakfiles in .mgf format using [ProteoWizard MSconvert](https://proteowizard.sourceforge.io/index.html). The script is designed to work with the windows version of msconvert.

The recalibrated .mgf files from this script may then be used as input for a crosslinking MS search with [xiSEARCH](https://github.com/Rappsilber-Laboratory/XiSearch).

Requirements:
- [ProteoWizard MSconvert](https://proteowizard.sourceforge.io/index.html)
- Seaborn
- Pandas

Before using, edit the path to MSconvert in config.py

Usage:

create a directory with the following structure:

    Top
    |
    |--rawfiles
    |--processed
    |--myfasta.fasta

Put your raw files in the "rawfiles" directory. myfasta.fasta is the sequence database you wish to recalibrate on. "processed" will contain your results

In command line (in windows, this may be powershell, anaconda prompt, or within an IDE), fron the top of the directory, run

    python /path/to/preprocessing_ms2recal.py  --db ./myfasta.fasta --input ./rawfiles --outpath ./processed --xiconf /path/to/resources/xi_linear_by_tryp.conf --config /path/to/config.py

--input folder containing.raw files or single file to process

--outpath directory for output, default is separate folder in input directory

--config path to config.py file

--xiconf  path to .conf file in resources directory

The .conf file is a xi config file set for a linear search with trypsin digestion. Other files may be chosen with different proteases and they are found in the "resources" directory. Documentation on editing config files with custom settings may be found [here](https://github.com/Rappsilber-Laboratory/XiSearch#setting-up-a-search-in-the-advanced-interface-and-editing-config-files) .



