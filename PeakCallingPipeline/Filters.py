############################################
# File Name : callers
#
# Purpose :
#
# Creation Date : 02-01-2015
#
# Last Modified : Mon May 11 19:25:33 2015
#
# Created By : ahaber
#
# Usage : 
############################################

#from util import info_green, info, mkdir_p, run_external, forscreen, \
    warn, error, sort_and_test, files_with_extension_in, make_bed_file, \
    normalise_peak_scores, has_file_extension
from filter import filter_bed_new
import os
import traceback
from shutil import copyfile
from pybedtools import BedTool
import sys
from pybedtools import create_interval_from_list


##
##paraclu_raw: the input bed file
def run_paraclu(paraclu_raw,paraclu_fin,min_peak_width,max_peak_width,min_density_rise,min_pos_with_data,min_sum,version="paraclu", settings=None, remove_temps=False):
    """
    Run paraclu/reclu peak calling. note, both reclu and paraclu use the same
    algorithm, but the reclu version is modified somehow (see Ohmiya et al. BMC Genomics 2014, 15:269)
    """
    if version == "paraclu":

        info("ParaClu filtering..")
        peaks_bt = BedTool(paraclu_raw)
        # too narrow or wide
        peaks_bt = filter_bed_new(peaks_bt,
                                  lambda x: (len(x) > min_peak_width and len(x) < max_peak_width),
                                  filter_name="Width")
        peaks_bt = peaks_bt.sort()
        initial_count = peaks_bt.count()

        # make sure s=True to only merge peaks on the same strand.
        peaks_bt = peaks_bt.merge(c=[4, 5, 6, 7, 8], o=["max", "max", "distinct", "min", "max"], s=True)  # d=min_resolution) #o=max -- keep the score of the
                                                       #  scoring peak, and 'collapse' keep a list of the names o
        info("[Merge] Initial: " + str(initial_count) + ", removed " + \
                      str(initial_count-peaks_bt.count()) + ", " + str(peaks_bt.count()) + " left.")

        # remove peaks that have less than the minimum density rise
        peaks_bt = filter_bed_new(peaks_bt, lambda x: float(x[7]) / float(x[6]) > min_density_rise,
                                  filter_name="DensityRise")

        peaks_bt = filter_bed_new(peaks_bt, lambda x: max(float(y) for y in x[3].split(",")) > min_pos_with_data,
                                  filter_name="MinPosWData")

        peaks_bt = filter_bed_new(peaks_bt, lambda x: float(x[4]) > min_sum,
                                  filter_name="SumData")

        # name the peaks (and overwrite the names of all the merged peaks.)
        # also, drop the two density columns as they are not standard BED format.
        name_id = "paraclu_peak"
        named_peaks = []
        for i, peak in enumerate(peaks_bt):
            if i < 2:
                print(peak)
            fields = [str(peak.chrom), str(peak.start), str(peak.end), name_id + str(i), str(peak.score), str(peak.strand)]
            named_peak = create_interval_from_list(fields)
            named_peaks.append(named_peak)
        named_peaks_bt = BedTool(fn=named_peaks)
        named_peaks_bt.saveas(paraclu_fin)


if __name__=="__main__":
	print("Read arguments!")
	output=sys.argv[1]
	max_peak_width=int(sys.argv[2])
	min_peak_width=int(sys.argv[3])
	min_density_rise=int(sys.argv[4])
	min_pos_with_data=int(sys.argv[5])
	min_sum=int(sys.argv[6])
	paraclu_raw=output+".raw.bed"
	paraclu_fin=output+".bed"
	print("Filter data!")
	run_paraclu(paraclu_raw,paraclu_fin,min_peak_width,max_peak_width,min_density_rise,min_pos_with_data,min_sum,version="paraclu")
	print("Filtering done!")
