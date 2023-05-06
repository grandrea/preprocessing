import os
import subprocess
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
pd.options.mode.chained_assignment = None
import ProteoFileReader
import sys


def xi_wrapper(arguments):
    xi = subprocess.Popen(arguments)
    xi.communicate()
    return


def run_xi_lin(peakfile, fasta, cnf, outpath, xipath, threads='1'):
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    elif os.path.isfile(outpath + '/xi_' + os.path.split(peakfile)[1].replace('.mgf', '.csv')):
        return

    xi_cmds = ['java', '-cp', os.path.join(os.path.dirname(os.path.realpath(__file__)), xipath),
               'rappsilber.applications.Xi', # + '/fastutil-8.1.0.jar;' + xipath + '/XiSearch.jar'
               '--fasta=' + fasta,
               '--xiconf=UseCPUs:' + threads,
               '--peaks=' + peakfile,
               '--config=' + cnf,
               '--output=' + outpath + '/xi_' + os.path.split(peakfile)[1].replace('.mgf', '.csv'),
               '--peaksout=%s_peaks.csv.gz' % peakfile[:len(peakfile) - 4]]

    print('calling ' + subprocess.list2cmdline(xi_cmds))
    xi = subprocess.Popen(xi_cmds) #, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    xi.communicate()


def get_ppm_error(xi_df, peaks_df, outfile):
    xi_df = xi_df[(xi_df.decoy == 0) & (xi_df['match score'] > 6)]
    median_err = np.median(xi_df['Precoursor Error'])
    try:
        fig, ax = plt.subplots()
        sns.histplot(data=xi_df, x='Precoursor Error')
        ax.axvline(median_err)
        plt.title("MS1 Error distribution \n median: " + str(median_err))
        plt.savefig(outfile)
        plt.close()
    except ZeroDivisionError:
        print(xi_df['Precoursor Error'][:5])

    if len(xi_df) < 50:
        try:
            print(os.path.split(outfile)[1] + ': Only %s PSMs found. No correction.' % (len(xi_df), median_err))
            return 0, 0
        except TypeError:
            print("no PSMs found. No correction")
            return 0, 0

    xi_ms2_df = peaks_df[peaks_df["IsPrimaryMatch"] == 1]
    xi_ms2_df["MS2Error_ppm"] = (xi_ms2_df["MS2Error"] * 10. ** 6) / xi_ms2_df["CalcMZ"]
    xi_ms2_df = xi_ms2_df.merge(xi_df[['Scan', 'Run', 'decoy']],
                                      left_on=['ScanNumber', 'Run'], right_on=['Scan', 'Run'], how='inner')
    xi_ms2_df = xi_ms2_df[(xi_ms2_df["MS2Error_ppm"] <= 30) & (xi_ms2_df["MS2Error_ppm"] >= -30)]
    median_err_ms2 = np.median(xi_ms2_df["MS2Error_ppm"])

    fig, ax = plt.subplots()
    sns.histplot(data=xi_ms2_df, x="MS2Error_ppm")
    ax.axvline(median_err_ms2)
    plt.xlabel("mass error")
    plt.title("MS2 Error distribution \n median: " + str(median_err_ms2))
    plt.ylabel("# of identifications")
    plt.xlim(-20, 20)
    plt.savefig(os.path.join(outfile.replace('MS1', "MS2")))
    plt.close()

    return median_err, median_err_ms2


def adjust_prec_mz(mgf_file, ms1_error, ms2_error, outpath):
    outfile = os.path.join(outpath, 'recal_' + os.path.split(mgf_file)[1])
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    elif os.path.isfile(outfile):
        return
    exp = ProteoFileReader.MGF_Reader()
    exp.load(mgf_file)

    out_writer = open(os.path.join(outfile), "w")
    for spectrum in exp:
        prec_mz_new = spectrum.getPrecursorMass()/(1 + ms1_error / 10. ** 6) # TODO wrong sign if newer version
        ms2_mass_new = spectrum.getPeaks()
        for i in range(0, len(ms2_mass_new[0])):
            ms2_mass_new[0][i] = ms2_mass_new[0][i] / (1 + ms2_error / 10. ** 6)

        if sys.version_info.major < 3:
            stavrox_mgf = """
MASS=Monoisotopic
BEGIN IONS
TITLE={}
PEPMASS={} {}
CHARGE={}+
RTINSECONDS={}
{}
END IONS     """.format(spectrum.getTitle(),
                            prec_mz_new, spectrum.getPrecursorIntensity() if spectrum.getPrecursorIntensity() > 0 else 0,
                                int(spectrum.charge), spectrum.getRT(),
                                "\n".join(["%s %s" % (i[0], i[1]) for i in ms2_mass_new if i[1] > 0]))
        else:
            stavrox_mgf = """
MASS=Monoisotopic
BEGIN IONS
TITLE={}
PEPMASS={} {}
CHARGE={}+
RTINSECONDS={}
{}
END IONS     """.format(spectrum.getTitle(),
                        prec_mz_new,
                        spectrum.getPrecursorIntensity() if spectrum.getPrecursorIntensity() > 0 else 0,
                        int(spectrum.charge), spectrum.getRT(),
                        "\n".join(["%s %s" % (mz, ms2_mass_new[1][i]) for i, mz in enumerate(ms2_mass_new[0]) if
                             ms2_mass_new[1][i] > 0]))
        out_writer.write(stavrox_mgf)


def main(mgf, fasta, xi_cnf, outpath, threads, xi_jar='./resources/XiSearch_1.6.745.jar', val_input=None):
    if not os.path.exists(outpath):
        os.makedirs(outpath)

    filename = os.path.split(mgf)[1]
    if val_input is None:
        # linear small search in Xi
        run_xi_lin(peakfile=mgf, fasta=fasta, cnf=xi_cnf, outpath=os.path.join(outpath), xipath=xi_jar, threads=threads)

        # evaluate results, get median ms1 error
        ms1_err, ms2_err = get_ppm_error(xi_df=pd.read_csv(os.path.join(outpath, 'xi_' + filename.replace('.mgf', '.csv'))),
                                         peaks_df=pd.read_csv(os.path.join(outpath, filename.replace('.mgf', '_peaks.csv.gz')), sep='\t', index_col=False, thousands=','),
                                         outfile=os.path.join(outpath, 'MS1_err_' + filename + '.png'))

        error_file = open(outpath + '/ms1_err.csv', 'a')
        error_file.write(filename + ',' + str(ms1_err) + '\n')
        error_file.close()

        error_file_ms2 = open(outpath + '/ms2_err.csv', 'a')
        error_file_ms2.write(filename + ',' + str(ms2_err) + '\n')
        error_file_ms2.close()
    else:
        ms1_input = pd.read_csv(val_input, header=None, index_col=0)
        ms1_err = ms1_input[ms1_input.index.str.contains('_'.join(filename.split('_')[1:]))].values[0][0]

    if ms1_err is not None: # shift all old m/z by value
        adjust_prec_mz(mgf_file=mgf, ms1_error=ms1_err, ms2_error=ms2_err, outpath=os.path.join(outpath))
