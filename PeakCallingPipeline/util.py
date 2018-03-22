############################################
# File Name : util
#
# Purpose : Utility methods for 5prime analysis code.
#
# Creation Date : 02-01-2015
#
# Last Modified : Sat Jan 23 11:25:08 2016
#
# Created By : ahaber
#
# Usage : 
############################################

from termcolor import colored
from time import strftime
import sys
import os
import errno
import subprocess
import threading
import csv
import pandas as pd
import numpy as np
from sklearn import preprocessing
from scipy import stats
from itertools import islice
from pybedtools import Interval, BedTool
import pysam
from shutil import copyfile
from distutils.version import LooseVersion

# def p_flatten(foo):
#     for x in foo:
#         if hasattr(x, '__iter__'):
#             for y in flatten(x):
#                 yield y
#         else:
#             yield x
#
#
# def flatten_list(list_of_lists):
#     return list(p_flatten(list_of_lists))


class Logger:
    def __init__(self):
        self.name = "logger"
    show_threads = False


def check_bedtools_version(min_version="2.0.0", max_version="2.22.1"):
    sys.stdout.write("Checking bedtools version...")
    v_string = subprocess.check_output(["bedtools", "--version"], stderr=subprocess.STDOUT)
    pieces = v_string.rstrip().split("v")
    version = pieces[1]
    if LooseVersion(version) > LooseVersion(max_version) or LooseVersion(version) < LooseVersion(min_version):
        sys.stdout.write("  Incompatible\n")
        error("Incompatible version of bedtools! [" + version + "] The pipeline requires"
              " BedTools " + str(min_version) + "-" + str(max_version))
        sys.exit(1)
    sys.stdout.write("  OK\n")
    return True

    #
    # if not found:
    #     error("Not found! Make sure Java " + str(version) + " is available!")
    #     sys.exit(1)
    # else:
    #     sys.stdout.write("  OK\n")


def check_java_version(version):
    sys.stdout.write("Checking for Java " + str(version) + "...") 
    found = ("1." + str(version)) in subprocess.check_output(["java", "-version"], stderr=subprocess.STDOUT)
    if not found:
        error("Not found! Make sure Java " + str(version) + " is available!")
        sys.exit(1)
    else:
        sys.stdout.write("  OK\n")



def add_columns(in_df, out_df, columns={"Condition": "StandardRNASeq", "Protocol": "Cufflinks"}):
    df = pd.read_csv(in_df, sep="\t")
    for colname, colvalue in columns.iteritems():
        df[colname] = colvalue
    info("Writing annotated df to: " + out_df)
    df.to_csv(out_df, sep="\t", index=False)


def combined_roc_curve(dataset_dir_1, dataset_dir_2, condition_1=None, protocol_1=None):
    """
    Utility only, not used in the pipeline. Produce a combined ROC curve from two seperate
    runs of the pipeline by loading 3 required files from the dataset directory of each
    run. If necessary, annotating which dataset the files came from by adding condition/
    protocol columns.
    """
    info("Loading files from " + dataset_dir_1)
    can_annotated = "can_annotated.txt"
    must_annotated = "must_annotated.txt"
    peaks_annotated = "peaks_annotated.txt"
    all_can = "can_all.txt"
    all_must = "must_all.txt"
    all_peaks = "peaks_all.txt"

    if condition_1 is not None or protocol_1 is not None:
        info("Adding condition/protocol data..")
        add_columns(dataset_dir_1 + "/merged_can_call.txt", can_annotated)
        add_columns(dataset_dir_1 + "/merged_must_call.txt", must_annotated)
        add_columns(dataset_dir_1 + "/merged_peaks.txt", peaks_annotated)

    merge_dbs([dataset_dir_2 + "/merged_can_call_combined.txt", can_annotated], merged_output=all_can)
    merge_dbs([dataset_dir_2 + "/merged_must_call_combined.txt", must_annotated], merged_output=all_must)
    merge_dbs([dataset_dir_2 + "/merged_can_call_combined.txt", peaks_annotated], merged_output=all_peaks)
    info("Generating combined ROC curve..")
    rscript = "/home/unix/ahaber/dev/adam/fiveprime/r/roc_curves.r"
    run_args = ["Rscript", rscript, all_must, all_can, all_peaks, "ROC_curve.pdf"]
    run_external(run_args)


def rename_col(df, index, newname):
    col_names = df.columns.values
    info("Renaming column " + str(index) +".. " + col_names[index] + " -> " + newname)
    new_names = [newname if index == i else x for i, x in enumerate(col_names)]
    df.columns = new_names


def is_iterable(obj, strings_return_true=False):
    if strings_return_true:
        try:
            some_object_iterator = iter(obj)
            return True
        except TypeError, te:
            # print obj, 'is not iterable'
            return False
    else:
        return hasattr(obj, '__iter__')


def head(strlist, n=10):
    for i, ln in enumerate(strlist):
        info("["+str(i)+"]: "+str(ln))
        if i > n:
            return


def line_prepender(filename, line):
    with open(filename, 'r+') as f:
        content = f.read()
        f.seek(0, 0)
        f.write(line.rstrip('\r\n') + '\n' + content)


def merge_dbs(input_files, merged_output=None, n=0, remove_headers=True):
    """
    Concatenate files and don't write the header at the beginning of each file.
    """
    # info("Merging dbs: " + str(input_files))
    if merged_output is None:
        merged_output = "merged_output.txt"
    if remove_headers:
        if n == 0:
            range_string = "FNR > 1"
        else:
            range_string = "FNR > 1 && FNR < " + str(n + 2)
        run_args = ["awk", range_string] + input_files
    else:
        run_args = ["cat"] + input_files
    with open(merged_output, 'w') as outfile:
        run_external(run_args, stdout=outfile, verbose=False)
    if remove_headers:
        # stick the header back into the merged file:
        with open(input_files[0], 'r') as f:
            header_line = f.readline()
        line_prepender(merged_output, header_line)
    # info("Merged dataset written to" + shortname(merged_output))


def shortname(absolutepath):
    base = os.path.basename(absolutepath)
    rval = os.path.splitext(base)[0]
    return rval


def sort_and_test(in_bed_file, out_bed_file=None, sort_with_unix=True):
    # info("Sorting and testing " + shortname(in_bed_file) + " -> " + shortname(out_bed_file)
    #      + " [sort_with_unix= " + str(sort_with_unix) + "]")
    if out_bed_file is None:
        dl = True
        out_bed_file = in_bed_file + ".tmp"
    else:
        dl = False
    if sort_with_unix:
        run_args = ["sort", "-k1,1", "-k2,2n", in_bed_file]
    else:
        run_args = ["sort-bed", "-e", in_bed_file]
    with open(out_bed_file, 'w') as outfile:
        run_external(run_args, stdout=outfile)
    if dl:
        os.remove(in_bed_file)
        copyfile(out_bed_file, in_bed_file)


def num_cols(input_file):

    with open(input_file) as f:
        reader = csv.reader(f, delimiter='\t', skipinitialspace=True)
        first_row = next(reader)
        num_cols = len(first_row)
        info("Number of columns in " + shortname(input_file) + ": "+str(num_cols))
        return num_cols


def end_chop(thestring, ending):
    if thestring.endswith(ending):
        return thestring[:-len(ending)]
    return thestring


def start_chop(thestring, start):
    if thestring.startswith(start):
        return thestring[len(start):]
    return thestring


def blank_lines(n):
    r = range(0, n)
    for i in r:
        print("")


def remove_lines_from_file(test_string, in_file, out_file="none"):
    f = open(in_file, "r")
    lines = f.readlines()
    f.close()
    if out_file == "none":
        f = open(in_file, "w")
    else:
        f = open(out_file, "w")
    for line in lines:
        if test_string not in line:
            f.write(line)
    f.close()


def remove_comments(in_file, out_file="none"):
    f = open(in_file, "r")
    lines = f.readlines()
    f.close()
    if out_file == "none":
        f = open(in_file, "w")
    else:
        f = open(out_file, "w")
    for line in lines:
        if not line.startswith("#"):
            f.write(line)
    f.close()


def get_subdirs(directory, ignore=None, ignore_system=True, verbose=False):
    """
    Count the number of subdirectories in a given directory. If any
    ignore argument is given, then directories with names
    containing that string (case-insensitive) are ignored
    :param dir:
    :param ignore:
    :return:
    """
    # sub_dirs = [x[0] for x in os.walk(directory)] ## this is ALL subdirectories, recursing
    sub_dirs = [name for name in os.listdir(directory)
                if os.path.isdir(os.path.join(directory, name))]
    rval = []
    for sd in sub_dirs:
        should_ignore = ignore is not None and ignore.upper() in sd.upper()
        sys_dir = sd.startswith(".") and ignore_system
        if not (sys_dir or should_ignore):
            rval.append(sd)
        if should_ignore and verbose:
            warn("Ignoring \'" + sd + "\' directory ..")
    return rval


def files_with_extension_in(direc, extension, verbose=False, full_path=False):
    f_list = os.listdir(direc)
    if extension.startswith("."):
        extension = extension[1:]
    if verbose:
        info("Searching for ." + extension + " files in /" + shortname(direc))
    if verbose:
        info("Excluding system files")
    f_list[:] = [tup for tup in f_list if not is_system_file(tup)]

    f_list[:] = [tup for tup in f_list if has_file_extension(tup, extension)]

    if full_path:
        f_list = [os.path.abspath(direc+"/"+f) for f in f_list]
    if verbose:
        info("Returning " + str(len(f_list)) + " files")
    return f_list


def write_ok_file(task_name, dir=os.getcwd()):
    open(dir+"/"+task_name+".ok", 'a').close()


def check_ok_file(task_name, dir=os.getcwd()):
    return os.path.isfile(dir+"/"+task_name+".ok")


def file_contains_string(test_string, in_file, out_file="none"):
    # info("Checking "+in_file)
    f = open(in_file, "r")
    lines = f.readlines()
    f.close()
    rval = False
    for line in lines:
        if test_string in line:
            rval = True
    return rval


def is_system_file(filename):
    base = os.path.basename(filename)
    short_name = os.path.splitext(base)[0]
    if short_name.startswith("."):
        return True
    else:
        return False


def extension_of(filename):
    base = os.path.basename(filename)
    extension=os.path.splitext(base)[1]
    return extension


def has_file_extension(filename, ext):
    base = os.path.basename(filename)
    periods = [pos for pos, char in enumerate(base) if char == "."]
    for p in periods:
        suf = base[p+1:]
        if suf == ext:
            return True
    return False


def is_bed_file(filename):
    base = os.path.basename(filename)
    extension=os.path.splitext(base)[1]
    if extension == ".bed":
        return True
    else:
        return False


# def count_reads(bam_file, qc_passed_only=False, idx=False):
#     """
#     Count the number of reads in a bam file.
#     :param filename:
#     :return:
#     """
#     info("Counting reads in " + shortname(bam_file) + "..")
#     if idx:
#         s = pysam.flagstat(bam_file)
#         if qc_passed_only:
#             rval = int(s[0].split()[0])
#         else:
#             rval = int(s[0].split()[0]) + int(s[0].split()[2])
#     else:
#
#         rval = reduce(lambda x, y: x + y, [eval('+'.join(l.rstrip('\n').split('\t')[2:]))
#                                            for l in pysam.idxstats(bam_file)])
#     return rval
#
#

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def lines_in_file(f):
    return float(sum(1 for line in open(f)))


def run_external(args,
                 stdout=None,
                 test=False,
                 process_output=False,
                 verbose=True):
    try:
        if verbose:
            print "\n                   >> " + flatten(args) + "\n"
        if not test:
            if process_output:
                child = subprocess.Popen(args, stdout=subprocess.PIPE)
                # child.communicate()
                output = child.stdout.read()
                child.wait()
                return output
            else:
                child = subprocess.Popen(args, stdout=stdout)
                exit_code = child.wait()
                return exit_code
        else:
            warn("Test only. Above job not submitted.")
    except KeyboardInterrupt:
        child.terminate()
        raise


def pad_interval(i, pad=100, verbose=False):
    """
    Given a pybedtools interval return a new interval padded *on both sides*
    with the given pad argument (bp).
    :param interval:
    :param pad:
    :return padded:
    """
    pad = int(float(pad))
    if verbose:
        info("Padding an interval of type: " + str(type(i)))
    padded = Interval(i[0], int(float(i.start)) - pad, int(float(i.end)) + pad, strand=i.strand)
    return padded


def load_bt(bedfile, verbose=False):
    if bedfile is None:
        error("Must provide a bed file to load!")
        return False
    if os.path.isfile(bedfile):
        if verbose:
            info("Loading bed tool from " + shortname(bedfile) + ".bed ..")
        bt = BedTool(bedfile)
        return bt
    else:
        error("File " + bedfile + " does not exist!")


def forscreen(time_in_sec):
    m, s = divmod(time_in_sec, 60)
    h, m = divmod(m, 60)
    if h != 0:
        return "{}hr, {}min, {}s".format(round(h,2),round(m,2),round(s,2))
    elif m !=0:
        return "{}min, {}s".format(round(m,2),round(s,2))
    return "{}s".format(round(s,2))


def flatten_list(l):
    return [item for sublist in l for item in sublist]


def flatten(in_array):
    rval = ""
    for p in in_array:
        rval = rval + " " + str(p)
    return rval


def highlight(string):
    print "\n\n\n       |==================================================================|"
    print colored("            "+string+"     ", 'green')
    print "       |==================================================================|\n\n"

def info_green(string):
    print colored("   ["+strftime("%c")+ "] INFO:  "+string, 'green')


def info_cyan(string):
    print colored("   ["+strftime("%c")+ "] INFO:  "+string, 'cyan')


def info_blue(string):
    print colored("   ["+strftime("%c")+ "] INFO:  "+string, 'blue')


def info_red(string):
    print colored("   ["+strftime("%c")+ "] INFO:  "+string, 'red')


def info(string, show_thread_name=False):
    if Logger.show_threads | show_thread_name:
        print "   ["+strftime("%c") + "]-[THREAD-"+threading.current_thread().name+"] INFO:  "+string
    else:
        print "   ["+strftime("%c") + "] INFO:  "+string


def bar():
    print "------------------------------------------------------------"


def warn(string,  show_thread_name=False):
    if Logger.show_threads | show_thread_name:
        print colored("   ["+strftime("%c")+ "]-[THREAD-"+threading.current_thread().name+"] WARN:  "+string, 'yellow')
    else:
        print colored("   ["+strftime("%c")+ "] WARN:  "+string, 'yellow')


def error(string, show_thread_name=False):
    if Logger.show_threads | show_thread_name:
        print colored("   ["+strftime("%c")+ "]-[THREAD-"+threading.current_thread().name+"] ERROR:  "+string, 'red')
    else:
        print colored("   ["+strftime("%c")+ "] ERROR:  "+string, 'red')


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


def convert_strings_to_num(list_str):
    rval = []
    for s in list_str:
        if is_number(s):
            rval.append(float(s))
        else:
            rval.append(s)
    return rval


def normalise_scores(scores,
                     linear=False,
                     n_decimal=8):
    """
    Scale the confidence of a peak caller to a value between 0 and 1. This is
    useful for generating ROC curves, and using the weighted sum consensus
    caller.
    :param scores:
    :return:
    """

    # print("Normalisng values: ")
    # print(scores)

    if not linear:
        warn("Using non-linear normalisation..")
        # scores = np.power(scores, 0.1)  # increase seperation between lower scoring peaks.

        # replaced the above ad hoc method with an attempt to transform the scores onto
        # a normal distribution, using a power transformation (box-cox)
        scores = stats.boxcox(scores)[0]

    min_max_scaler = preprocessing.MinMaxScaler()
    ##UPDATED## minmax = np.round(min_max_scaler.fit_transform(scores), decimals=8)
    minmax = np.round(min_max_scaler.fit_transform(scores.reshape(1,-1)), decimals=8)
    normed_scores = np.round(minmax, decimals=n_decimal)
    # this works but seemed super slow (only tested once)
    #minmax = [(s_i - min(scores)) / (min(scores) - max(scores)) for s_i in scores]
    info("Max score: "
         + str(np.max(normed_scores))+" \n Min score:"+ str(np.min(normed_scores))
         + " \n Mean: " + str(np.mean(normed_scores)) + "\n StdDev: " + str(np.std(normed_scores)))
    return normed_scores


def normalise_peak_scores(bed_file, output=None, n_decimal = 8):
    ##UPDATED!##df = pd.read_csv(bed_file, sep="\t", header=False)
    df = pd.read_csv(bed_file, sep="\t", header=None)
    scores = df.ix[:, 4].values.astype(float)
    # head(scores, n = 100)
    normed_scores = normalise_scores(scores)
    normed_scores=np.ravel(normed_scores)
    print(normed_scores[1:10])
    df.ix[:, 4] = pd.Series(normed_scores)
    if output is None:
        df.to_csv(bed_file, sep="\t", index=False, float_format='%.'+str(n_decimal)+'f', header=False)
    else:
        df.to_csv(output, sep="\t", index=False, float_format='%.'+str(n_decimal)+'f', header=False)


def make_bed_file(tab_delimited_f,
                  output_bed=None,
                  name_id=None,
                  column_mapping={"chr":1,
                                  "start":2,
                                  "end":3,
                                  "name":4,
                                  "score":5,
                                  "strand":6},
                  normalise_scores=False):

    if output_bed is None:
        output_bed = shortname(tab_delimited_f) + ".bed"
    info("Making bed file from " + shortname(tab_delimited_f) + " --> " + shortname(output_bed))
    bed_columns = ["chr","start","end","name","score","strand"]
    with open(output_bed, 'wb') as output_f:
        writer = csv.DictWriter(output_f, bed_columns, delimiter="\t")
        with open(tab_delimited_f, 'r') as input_f:
            reader = csv.reader(input_f, delimiter='\t')
            for row_index, row in enumerate(reader):
                if row[0].startswith("#"):
                    continue
                new_row = {}
                for i, col in enumerate(bed_columns):
                    index_in_input_file = column_mapping[col]
                    if index_in_input_file == -1:
                        # the input file doesn't have a value for this column.
                        if col == "name" and name_id is not None:
                            new_row[col] = name_id + "_" + str(row_index)
                        else:
                            new_row[col] = ""
                    else:
                        new_row[col] = row[index_in_input_file]
                    #print "For bed column " + str(i) + " ("+col+"), use column " + str(index_in_input_file) + " from input"
                writer.writerow(new_row)
    if normalise_scores:
        normalise_peak_scores(output_bed)

#
# pybedtools .each() method appears buggy.
# def remove_chr(feature):
#     """
#     Called on each feature in a bedtool to remove 'chr'
#     """
#     if "chr" in feature[0]:
#         feature[0] = feature[0].replace('chr', '')
#     return feature


def take(n, iterable):
    "Return first n items of the iterable as a list"
    return list(islice(iterable, n))


def remove_chrs(bt):
    """
    Given a bedtool with 'chr' denoting the chromosomes, remove them
    :param bt:
    :return:
    """
    intervals = []
    for i in bt:
        if "chr" in i.chrom:
            i.chrom = i.chrom.replace('chr', '')
        intervals.append(i)
    return BedTool(fn=intervals)


def show_value(s):
    """
    Convert unicode to str under Python 2;
    all other values pass through unchanged
    """
    if sys.version_info.major == 2:
     if isinstance(s, unicode):
         return str(s)
    return s


def query_yes_no(question, default="yes"):
    """Ask a yes/no question via raw_input() and return their answer.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).

    The "answer" return value is True for "yes" or False for "no".
    """
    valid = {"yes": True, "y": True, "ye": True,
             "no": False, "n": False}
    if default is None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        choice = raw_input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "
                             "(or 'y' or 'n').\n")


def is_stranded(bt):
    """
    Return true if the bedtool contains strand info, false
    otherwise.
    :param bt:
    :return:
    """
    for interval in bt:
        strand = interval.strand
        if strand == "+" or strand == "-":
            return True
        else:
            warn("Strand field is populated but not with +/- (it has: " + strand + ")")
            return False
