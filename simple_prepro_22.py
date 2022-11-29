import os
import numpy as np
import subprocess
from multiprocessing import Pool
import sys
import re
import getopt
from pyteomics import mzml
from functools import partial
import ProteoFileReader
import mass_recal_ms2
import zipfile
import glob


def read_cmdline():
    try:
        opts, args = getopt.getopt(sys.argv[1:], '', ['input=', 'config=', 'outpath=', 'db=', 'xiconf=', 'shiftcsv=', 'skip_recal='])
    except getopt.GetoptError:
        print('preprocessing.py --input <folder or single file to process> '
              '--outpath <directory for output, default is separate folder in input directory> '
              '--config <path to config file> '
              '--db <path to database to search for recalibration>'
              '--xiconf <path to xi config to use for recalibration>')
        sys.exit()
    recal = True
    recal_conf = {}
    for opt, arg in opts:
        if opt == '--input':
            input_arg = arg
        elif opt == '--outpath':
            outdir = arg
        elif opt == '--config':
            config = arg
        elif opt == '--db':
            recal_conf['db'] = arg
        elif opt == '--xiconf':
            recal_conf['xiconf'] = arg

    if 'input_arg' not in locals() or 'config' not in locals():
        print('preprocessing.py --input <folder or single file to process> '
              '--outpath <directory for output, default is separate folder in input directory> '
              '--config <path to config file> '
              '--db <path to database to search for recalibration>'
              '--xiconf <path to xi config to use for recalibration>')
        sys.exit()
    # if no outdir defined use location of input
    if 'outdir' not in locals() and os.path.isdir(input_arg):
        outdir = os.path.join(input_arg, 'processed')
    elif 'outdir' not in locals() and not os.path.isdir(input_arg):
        outdir = os.path.join(os.path.split(input_arg)[0], 'processed')

    return input_arg, outdir, config, recal_conf, recal


def process_file(filepath, outdir, mscon_settings, split_acq, detector_filter, mscon_exe, cihcd_ms3=False):

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    conv_cmds = mscon_cmd(filepath=filepath, outdir=outdir, settings=mscon_settings, mgf=not split_acq)

    if len(conv_cmds) > 0:
        msconvert = subprocess.Popen([mscon_exe] + conv_cmds)
        msconvert.communicate()

    filename = os.path.split(filepath)[1]
    mzml_file = os.path.join(outdir, filename[:filename.rfind('.')] + '.mzML')

    if cihcd_ms3:
        cihcd_spectra = generate_cihcd_spectra(mzml_file)
        write_mgf(spectra=cihcd_spectra, outfile=os.path.join(outdir, 'CIhcD_ms3_' + filename[:filename.rfind('.')] + '.mgf'))

    if split_acq:
        splitted_spectra = split_mzml(mzml_file, detector_filter)

        for acq in splitted_spectra:
            write_mgf(spectra=splitted_spectra[acq],
                      outfile=os.path.join(outdir, acq + '_' + filename[:filename.rfind('.')]+'.mgf'))


if __name__ == '__main__':
    # read cmdline arguments / get deafult values
    input_arg, outdir, config_path, recal_conf, recal = read_cmdline()
    try:
        execfile(config_path)
    except NameError:
        exec(open(config_path).read())

    # get files in directory
    if os.path.isdir(input_arg):
        file_list = glob.glob(input_arg)
        full_paths = [x for x in file_list if (not os.path.isdir(x)) and ((x[-4:] == '.raw') or (x[-4:] == '.mgf'))]

        print("""file input:
        {}
        """.format('\n'.join(full_paths)))

        # start msconvert for conversion and peak filtering
        pool = Pool(processes=nthr)
        pool.map(partial(process_file, outdir=outdir, mscon_settings=mscon_settings, split_acq=split_acq,
                         detector_filter=detector_filter, mscon_exe=msconvert_exe), full_paths)
        pool.close()
        pool.join()
