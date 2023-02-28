import getopt, sys

import datetime

#sys.path.append('/home/users/cjwalden/my-packages')
import kepler_utils as kepler

try:
    opts, args = getopt.getopt(sys.argv[1:], "d:i:o:", ["date=","inpath=","outpath="])
except getopt.GetoptError as err:
        # print help information and exit:
        print(err) # will print something like "option -a not recognized
        sys.exit(2)

data_date = datetime.datetime.now()
datestr = data_date.strftime('%Y%m%d')

inpath = '/Users/cw66/Data/kepler/'
outpath = '/Users/cw66/Data/ncas-mobile-radar-ka-band-1/'

yaml_project_file = '/Users/cw66/git/kepler-radar-utils/amof_projects.yml'
yaml_instrument_file = '/Users/cw66/git/kepler-radar-utils/amof_radars.yml'

tracking_tag = 'AMOF_20220922221548';

for o, a in opts:
    if o == "-d":
        datestr = a;
    elif o == "-i":
        inpath = a;
    elif o == "-o":
        outpath = a;
    else:
        assert False, "unhandled option"

kepler.process_kepler(datestr,inpath,outpath,yaml_project_file,yaml_instrument_file,tracking_tag);

