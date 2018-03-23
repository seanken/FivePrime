
from pybedtools import BedTool
import sys
import traceback


__author__ = 'ahaber'


def filter_bed_new(bt, filter_function, verbose=True, filter_name=None, max_error=5):
    initial_count = bt.count()
    intervals = []
    line_count = 0
    for interval in bt:
        try:
            if filter_function(interval):
                intervals.append(interval)
            # if line_count % 10000 == 0 and verbose:
            #     print("Processed " + str(line_count) + " intervals ["+str(round(100*float(line_count)/float(initial_count), 2))+"%]")
        except ValueError, e:
            error(filter_name + " filtering failed for line #" + str(line_count))
            error(str(e))
            print(traceback.format_exc())
            print(sys.exc_info()[0])
            error(str(interval))
            for i in range(1, 8):
                error("Field[" + str(i) + "]=" + str(interval[i]))
            sys.exit(1)


        line_count += 1

    # print("Building bedtool from " + str(len(intervals)) + " intervals..")
    rval = BedTool(fn=intervals)
    # print("Built bedtool.")
    if verbose:
        info_string = "Initial: " + str(initial_count) + ", removed " + \
                      str(initial_count-rval.count()) + ", " + str(rval.count()) + " left."
        if filter_name is not None:
            info_string = "[" + filter_name + "] " + info_string
    return rval


