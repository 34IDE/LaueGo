#!/usr/bin/env python
##!/usr/bin/env /APSshare/bin/python
##/usr/bin/env python

import os
import sys
from math import floor
from time import asctime, strptime

fraction = (12, 27, 42, 57)
#
# /net/ord/share1/burt/2005/12/07/01-42-00.gz
# /net/sec33.xor.aps.anl.gov/export/sector33_34/burt/2009/03/05/11-27-01.gz
burt_path = '/net/sec33.xor.aps.anl.gov/export/sector33_34/burt'


def main(pv, start_date, period=1.0, end_date=-1):
    if end_date == -1:
        end_date = start_date

    start = [int(x) for x in start_date.split("-")]
    end = [int(x) for x in end_date.split("-")]
    current = start[:]
    current.append(0)
    current.append(0)

    f_period = float(period)
    
    period_array = [int(floor(f_period)), int((f_period - floor(f_period)) / 0.25)]
    
    while current[:3] <= end:
	for lastDigit in range(3):
            ## Construct the path and search the file if it exists
            #                                                             year        month       day         hour        minute                second
            current_path = "%s/%i/%02i/%02i/%02i-%i-0%i.gz" % (burt_path, current[0], current[1], current[2], current[3], fraction[current[4]], lastDigit)
            if os.path.isfile(current_path):
                temp = os.popen("zgrep %s %s" % (pv, current_path))
                cmd_output = temp.readline()
                temp.close
                del temp
                if cmd_output != "":
    		## originally did the following, but it doesn't work for eps event details
                    #pv_value_string = cmd_output.split(" ")[7]
                    pv_value_string = " ".join(cmd_output.split(" ")[7:])
                    ordinary_time_string = "%i-%02i-%02i %02i:%i:00" % (current[0], current[1], current[2], current[3], fraction[current[4]])
                    reformatted_time_string = asctime(strptime(ordinary_time_string, "%Y-%m-%d %H:%M:%S"))
                    try:
                        float(pv_value_string)
                    except ValueError:
                        print "%s\t%s\t%s\t%s" % (current_path, reformatted_time_string, pv, pv_value_string)
                    else:
                        print "%s\t%s\t%s\t%f" % (current_path, reformatted_time_string, pv, float(pv_value_string))
                break

        ## Increment the current array by the desired period
        current[3] = current[3] + period_array[0]
        current[4] = current[4] + period_array[1]
        
        ## Handle the overflow of each of the time units
        # Hours
        if current[4] > 3:
            current[3] = current[3] + 1
            current[4] = current[4] - 4
        # Days
        if current[3] > 23:
            days, hours = divmod(current[3],24)
            current[3] = hours
            current[2] = current[2] + days
        # Months
        if current[2] > 31:
            months, days = divmod(current[2],31)
            current[2] = days
            current[1] = current[1] + months
        # Years
        if current[1] > 12:
            years, months = divmod(current[1],12)
            current[1] = months
            current[0] = current[0] + years

### Command line arguments
# sys.argv[0] - burt_tool.py
# sys.argv[1] - 438d:m58:c0:m1.DVAL
# sys.argv[2] - 2005-10-02
### arguments below are optional
# sys.argv[3] - 2005-11-06
# sys.argv[4] - 0.25

count = 0

for x in sys.argv:
    count += 1

if count == 3:
    main(sys.argv[1], sys.argv[2])
elif count == 4:
    ## Allow entering either period or end_date
    try:
        float(sys.argv[3])
    except ValueError:
        main(sys.argv[1], sys.argv[2], end_date=sys.argv[3])
    else:
        main(sys.argv[1], sys.argv[2], sys.argv[3])
elif count == 5:
    ## If both period AND end_date are specified, period must come first
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
else:
    print ""
    print "Usage: burt_tool.py PV START_DATE [PERIOD] [END_DATE]"
    print ""
    print "PV         = process variable name"
    print "START_DATE = starting date in YYYY-MM-DD format"
    print "PERIOD     = time (in hours) between displayed values (default 1, minimum 0.25)"
    print "END_DATE   = ending date in YYYY-MM-DD format"
    print ""
    print "Note: period and end_date are optional arguments"
    print ""

