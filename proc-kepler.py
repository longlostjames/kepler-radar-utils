import getopt, sys

import datetime

sys.path.append('/home/users/cjwalden/my-packages')
import kepler_utils as kepler

try:
    opts, args = getopt.getopt(sys.argv[1:], "d:p:c:", ["date=","path=","campaign="])
except getopt.GetoptError as err:
        # print help information and exit:
        print(err) # will print something like "option -a not recognized
        sys.exit(2)

data_date = datetime.now()
datestr = data_date.strftime('%Y%m%d')
data_path = '/gws/pw/j07/ncas_obs_vol2/cao/raw_data/ncas-radar-mobile-ka-band-1/data/'


for o, a in opts:
    if o == "-d":
        datestr = a;
    elif o == "-p":
        data_path = a
    elif o == "-c":
        campaign = a.lower()
    else:
        assert False, "unhandled option"

dateyr = datestr[0:4]

if campaign == "long-term":
    mmclxpath = os.path.join(data_path,"long-term",dateyr,datestr);
else:
    mmclxpath = os.path.join(data_path,"campaign",campaign,datestr);

os.chdir(mmclxpath);
files = [os.path.join(mmclxpath,f) for f in glob.glob('*.mmclx')]

print(files);