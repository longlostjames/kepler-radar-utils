import getopt, sys, os
import glob


from datetime import datetime

#sys.path.append('/home/users/cjwalden/my-packages')
import kepler_utils as kepler

try:
    opts, args = getopt.getopt(sys.argv[1:], "d:i:c:o:", ["date=","inpath=","campaign=","outpath="])
except getopt.GetoptError as err:
        # print help information and exit:
        print(err) # will print something like "option -a not recognized
        sys.exit(2)

data_date = datetime.now()
datestr = data_date.strftime('%Y%m%d')
inpath = '/gws/pw/j07/ncas_obs_vol2/cao/raw_data/ncas-radar-mobile-ka-band-1/data/'
outpath = '/gws/pw/j07/ncas_obs_vol2/cao/processing/ncas-radar-mobile-ka-band-1/'


for o, a in opts:
    if o == "-d":
        datestr = a;
    elif o == "-i":
        inpath = a;
    elif o == "-c":
        campaign = a.lower();
    elif o == "-o":
        outpath = a;
    else:
        assert False, "unhandled option"

dateyr = datestr[0:4]

if campaign == "long-term":
    mmclxpath = os.path.join(inpath,"long-term","mom",dateyr,datestr);
    outpath = os.path.join(outpath,"long-term",dateyr,datestr);
else:
    mmclxpath = os.path.join(inpath,"campaign",campaign,"mom",datestr);
    outpath = os.path.join(outpath,campaign,datestr);
os.chdir(mmclxpath);
files = [os.path.join(mmclxpath,f) for f in glob.glob('*.mmclx')]

print(files);

print(outpath);


if not os.path.exists(outpath):
    os.makedirs(outpath)