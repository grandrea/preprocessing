# Preprocessing for xiSEARCH

This set of tools recalibrates MS1 and MS2 spectra based on mass error of a linear proteomics search. It uses [xiSEARCH](https://github.com/Rappsilber-Laboratory/XiSearch) to perform the linear search. This is usually done as the first step in the xiSEARCH workflow, prior to a crosslinking MS search, to improve identifications and understand what MS1 and MS2 error tolerances one should set. It first converts Thermo .raw files into peakfiles in .mgf format using [ProteoWizard MSconvert](https://proteowizard.sourceforge.io/index.html). The script can call either a native `msconvert.exe` installation or a Singularity image that exposes `msconvert`.

The recalibrated .mgf files from this script may then be used as input for a crosslinking MS search with [xiSEARCH](https://github.com/Rappsilber-Laboratory/XiSearch).

If you use this preprocessing script, please cite [Lenz *et al.*, Nat. Comm. 2021](https://www.nature.com/articles/s41467-021-23666-z).

#### Requirements:
- [ProteoWizard MSconvert](https://proteowizard.sourceforge.io/index.html) for native mode
- [Singularity](https://sylabs.io/docs/) and a `pwiz-skyline.sif`-style image for container mode
- Seaborn
- Pandas

#### Usage:

Before usage, edit `config.py` to choose the msconvert backend.

Native mode:

    msconvert_mode = 'native'
    msconvert_exe = r'/path/to/msconvert.exe'

Singularity mode:

    msconvert_mode = 'singularity'
    singularity_exe = 'singularity'
    singularity_image = '/path/to/pwiz-skyline.sif'
    singularity_wine_bind_target = '/path/to/wine-prefix'
    singularity_tmp_root = '/path/to/wine-temp-root'
    msconvert_nthr = 1

In singularity mode, `preprocessing_ms2recal.py` keeps the same CLI, but it runs msconvert through:

    singularity run -B <input/output/db/xiconf dirs> -B <tempdir>:<wine bind target> /path/to/pwiz-skyline.sif msconvert ...

The directories implied by `--input`, `--outpath`, `--db`, and `--xiconf` are bound automatically. These paths must already be valid on the host running Singularity; the script does not translate Windows paths into Linux paths.
If the host `/tmp` is too small for Wine initialization, set `singularity_tmp_root` to a larger scratch location.
When using Singularity/Wine, `msconvert_nthr = 1` is recommended even if `nthr` is higher, because the conversion step is not reliably parallel-safe in this container setup.

Create a directory with the following structure (this directory tree is not required, it's just to make the paths in the command clearer):

    Top
    |
    |--rawfiles
    |--processed
    |--myfasta.fasta

Put your raw files in the "rawfiles" directory. myfasta.fasta is the sequence database you wish to recalibrate on. "processed" will contain your results

In command line, from the top of the directory, run

    python /path/to/preprocessing_ms2recal.py  --db ./myfasta.fasta --input ./rawfiles --outpath ./processed --xiconf /path/to/resources/xi_linear_by_tryp.conf --config /path/to/config.py

--input folder containing.raw files or single file to process

--db the .fasta file containing the sequences to be searced.

--outpath directory for output, default is separate folder in input directory

--config path to config.py file (edited to choose either native or singularity msconvert)

--xiconf  path to .conf file in resources directory

The .conf file is a xi config file set for a linear search with trypsin digestion. Other files may be chosen with different proteases and they are found in the "resources" directory. Documentation on editing config files with custom settings may be found [here](https://github.com/Rappsilber-Laboratory/XiSearch#setting-up-a-search-in-the-advanced-interface-and-editing-config-files) .


The output directory contains several files:
- peakfiles recalibrated according to the ms1 and ms2 errors (recal_*.mgf) **these are the files to be used in a crosslinking MS search by xiSEARCH**
- .csv files with the average ms1 and ms2 errors per raw file
- images of the ms1 and ms2 error distributions - these should be symmetric gaussian shapes. If they are not, something may be wrong with the search or the acquisition. 
- peakfiles without any error recalibration (which retain the original file name)
- .csv file with the xiSEARCH output
- The error distributions may then be used to understand and set the tolerances for ms1 and ms2 matching in a subsequent crosslinking MS search in xiSEARCH.

Depositing into ProteomeXChange repositories:
Typically, the recalibrated .mgf files are included in the deposition of crosslinking MS results in PRIDE, JPost or other ProteomeXChange repositories.
