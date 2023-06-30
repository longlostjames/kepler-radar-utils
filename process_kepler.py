import getopt, sys, os

import datetime

from pathlib import Path
homepath = Path.home()

kepler_rawpath = '/gws/pw/j07/ncas_obs_vol2/cao/raw_data/ncas-radar-mobile-ka-band-1/data'

#sys.path.append('/home/users/cjwalden/my-packages')
import kepler_utils as kepler

try:
    opts, args = getopt.getopt(sys.argv[1:], "d:i:o:c:", ["date=","inpath=","outpath="])
except getopt.GetoptError as err:
        # print help information and exit:
        print(err) # will print something like "option -a not recognized
        sys.exit(2)

data_date = datetime.datetime.now()
datestr = data_date.strftime('%Y%m%d')

tracking_tag = 'AMOF_20220922221548';

campaign = 'woest';


yaml_project_file = os.path.join(homepath,'amof_campaigns','{campaign}_project.yml')
yaml_instrument_file = os.path.join(homepath,'amof_instruments','amof_radars.yml')


for o, a in opts:
    if o == "-d":
        datestr = a;
    elif o == "-i":
        inpath = a;
    elif o == "-o":
        outpath = a;
    elif o == "-c":
        campaign = a;
    else:
        assert False, "unhandled option"

inpath = os.path.join(kepler_rawpath,'campaign',campaign,'mom');
outpath = os.path.join('/gws/pw/j07/ncas_obs_vol2/cao/processing/ncas-radar-mobile-ka-band-1',campaign);

print(tracking_tag);

kepler.process_kepler(datestr,inpath,outpath,yaml_project_file,yaml_instrument_file,tracking_tag);

