# utility module for working with Apollo HFE data from various sources

import numpy as np
import pandas as pd
import datetime as dt
import astropy.time as atime
import os
import copy
import re

def load_hfe(pathname='.'):
    data={}

    # uncorrected 2005 NSSDC data
    for m in ['a15','a17']:
        if not m in data.keys():
            data[m]={}
        for p in ['p1','p2']:
            data[m][p]={}
            for s in [1,2,3,4,5]:
                filename = '{pathname}/{m}/{m}_hfe_{p}_{s}.tab'.\
                                    format(pathname=pathname,m=m,p=p,s=s)
                columns = (['Time','HTR','TREF','TC1','TC2','TC3','TC4']
                                                    if s==3 else ['Time','dT','T'])
                data[m][p][s] = pd.read_csv(filename,skiprows=2,
                    names=columns,delim_whitespace=True,skipinitialspace=True)
    return data

def flag_missing_hfe(data):
    for m in data.keys():
        for p in data[m].keys():
            for s in data[m][p].keys():
                name = '{m}:{p}:{s}'.format(m=m,p=p,s=s)
                n = np.array(data[m][p][s]['Time']).shape[0]
                data[m][p][s]['flags']=pd.Series(
                    np.zeros(n,dtype=np.int16),index=data[m][p][s].index)

                # Flag rows with any missing data values of -9999. Nagihara data don't use 
                # this convention (they simply don't include any missing points), but we 
                # create the flags column while we're at it.
                
                # we do *not* exclude values where HTR data is missing; HTR data
                # is essentially unusable, but thermocouple/reference bridge data isn't. 

                if s==3: 
                    data[m][p][s].loc[
                        (data[m][p][s]['Time']==-9999) |
                        (data[m][p][s]['TC1']==-9999) |
                        (data[m][p][s]['TC2']==-9999) |
                        (data[m][p][s]['TC3']==-9999) |
                        (data[m][p][s]['TC4']==-9999),'flags']+=0b1
                else:
                    data[m][p][s].loc[
                        (data[m][p][s]['dT']==-9999) |
                        (data[m][p][s]['T']==-9999) |
                        (data[m][p][s]['Time']==-9999),'flags']+=0b1

# flags time-inverted pairs of points in official data, then sorts all sets by time, which
# also fixes some misplaced time ranges in the Nagihara paper data. checks all sets, 
# although this is a bit wasteful.

def manage_disordered_hfe(data):
    for m in data.keys():
        for p in data[m].keys():
            for s in data[m][p].keys():
                if m=='a15' or m=='a17':
                    f=np.array(data[m][p][s]['flags'])
                    for i in np.arange(data[m][p][s]['Time'].size-1):
                        if data[m][p][s]['Time'][i+1]-data[m][p][s]['Time'][i] <= 0 and\
                        data[m][p][s]['Time'][i] != -9999 and\
                        data[m][p][s]['Time'][i+1] != -9999:
                            f[i]+=0b10000000000
                            f[i+1]+=0b10000000000
                    data[m][p][s]['flags']=f
                data[m][p][s]=data[m][p][s].sort_values(by='Time').reset_index(drop=True)

# functions for writing out reduced sets.

def write_clean_hfe(data,outpath='.',version=''):
    if not os.path.exists(outpath):
        os.mkdir(outpath)
        os.mkdir(outpath+'/a15')
        os.mkdir(outpath+'/a17')
    data_clean_out=copy.deepcopy(data)
    for m in data_clean_out.keys():
        for p in data_clean_out[m].keys():
            for s in data_clean_out[m][p].keys():
                if 'dT_corr' in data_clean_out[m][p][s].columns:
                    data_clean_out[m][p][s]=data_clean_out[m][p][s].reindex(columns=\
                                    ['Time','T','dT','dT_corr','flags'])
                elif s==3:
                    data_clean_out[m][p][s]=data_clean_out[m][p][s].reindex(columns=\
                                    ['Time','HTR','TREF','TC1','TC2','TC3','TC4','flags'])
                else:
                    data_clean_out[m][p][s]=data_clean_out[m][p][s].reindex(columns=\
                                    ['Time','T','dT','flags'])

                # convert to the IBM 1130-equivalent number format, but retain millisecond precision
                # for Nagihara data

                if m[3:4]=='_':
                    for column in data_clean_out[m][p][s].columns:
                        if column=='Time':
                            data_clean_out[m][p][s][column]=data_clean_out[m][p][s][column].\
                                                apply("{:.11E}".format).str.pad(17,'right')
                        else:
                            data_clean_out[m][p][s][column]=data_clean_out[m][p][s][column].\
                                                apply("{:.7E}".format).str.pad(14,'right')
                    # these profoundly ugly regex substitutions are required because
                    # of deficiencies in fixed-width table formatting in pandas, exacerbated
                    # by a regression in pandas.to_string 0.25 that inserts leading 
                    # spaces in non-indexed output.
                    table=re.sub(r'\n\s(\d|-)',r'\n\1',data_clean_out[m][p][s].to_string\
                            (index=False,justify='left'))
                    table=re.sub(r' (?= (\d|-))', r'',table)
                    # create correct line endings when script run from any major OS
                    table=re.sub(r'(\n)|(\r\n)|(\r)',r'\r\n',table)
                    with open('{outpath}/{m}{p}f{s}{v}.tab'.format\
                        (outpath=outpath+'/'+m[0:3],m=m,p=p,s=s,v=version), "w")\
                        as output_file: 
                        print(table, file=output_file)
                else:
                    for column in data_clean_out[m][p][s].columns:
                        data_clean_out[m][p][s][column]=data_clean_out[m][p][s][column].\
                                                apply("{:.7E}".format).str.pad(14,'right')
                    table=re.sub(r'\n\s(\d|-)',r'\n\1',data_clean_out[m][p][s].to_string\
                            (index=False,justify='left'))
                    table=re.sub(r' (?= (\d|-))', r'',table)
                    table=re.sub(r'\n',r'\r\n',table)
                    with open('{outpath}/{m}{p}f{s}{v}.tab'.format(outpath=outpath+'/'+\
                        m[0:3],m=m,p=p,s=s,v=version), "w") as output_file:
                        print(table, file=output_file)

def write_split_hfe(data,outpath='.',version=''):
    data_split_out=copy.deepcopy(data)
    if not os.path.exists(outpath):
        os.mkdir(outpath)
    for m in data_split_out.keys():
        for p in data_split_out[m].keys():
            for s in data_split_out[m][p].keys():
                data_split_out[m][p][s]['Time']=data_split_out[m][p][s]['Time'].\
                                                dt.strftime("%Y-%m-%dT%H:%M:%S.%f").str.slice(0,-3)+"Z"
                data_split_out[m][p][s].to_csv('{outpath}/{m}{p}f{s}{v}_split.tab'.format(
                    outpath=outpath,m=m,p=p,s=s,v=version),index=False,line_terminator='\r\n')

def write_deep_hfe(data,outpath='.',version=''):
    data_deep_out=copy.deepcopy(data)
    if not os.path.exists(outpath):
        os.mkdir(outpath)
    for m in data_deep_out.keys():
        for p in data_deep_out[m].keys():
            data_deep_out[m][p]['Time']=data_deep_out[m][p]['Time'].\
                                    dt.strftime("%Y-%m-%dT%H:%M:%S.%f").str.slice(0,-3)+"Z"
            data_deep_out[m][p].to_csv('{outpath}/{m}{p}{v}_depth.tab'.format(
                outpath=outpath,m=m,p=p,v=version),index=False,line_terminator='\r\n')

# Functions for interpreting data released by Nagihara et. al along with their 2018 paper "Examination of 
# the Long-Term Subsurface Warming Observed at the Apollo 15 and 17 Sites Utilizing the Newly Restored 
# Heat Flow Experiment Data From 1975 to 1977." Output intended primarily as intermediate data for 
# further correction and reduction by other utilities in this module.   
    
# The original reduced Apollo HFE data uses an epoch time format: milliseconds from 24 hours before the beginning of the 
# mission's start year. This is December 31, 1970 for Apollo 15, and December 31, 1971 for Apollo 17. 
# These functions convert between Gregorian time and mission epoch time.

# epoch-to-Gregorian functions. placing in TAI to avoid leap second weirdness. excessive digits of precision are for parity with
# internal astropy values.

def a15_to_greg(x):
    epoch=dt.datetime(1970,12,31,0,0,8,943570)
    return(epoch+dt.timedelta(milliseconds=x))

def a17_to_greg(x):
    epoch=dt.datetime(1971,12,31,0,0,9,889650)
    return(epoch+dt.timedelta(milliseconds=x))

# Gregorian-to-epoch functions (assume TAI input)

def a15_time(x):
    if not x==None:
        epoch=dt.datetime(1970,12,31,0,0,8,943570)
        return (x-epoch).total_seconds()*1000

def a17_time(x):
    if not x==None:
        epoch=dt.datetime(1971,12,31,0,0,9,889650)
        return (x-epoch).total_seconds()*1000

# utility functions for converting between TAI and UTC. this allows us to use astropy to deal with leap 
# seconds without later doing large table comparisons between astropy Time objects (much slower than 
# using datetime).

# 9999 flags empty rows in nagihara 2018 spreadsheet.

def tai_to_utc(x):
    return atime.Time(x,scale="tai").utc.datetime

def utc_to_tai(x):
    if not x==dt.datetime(9999,1,1,0,0,0):
        return atime.Time(x,scale="utc").tai.datetime
    
# silly utility functions for vectorizing over Nagihara datelists. 9999 flags empty rows.

def seconds_interval(x):
        if not np.isnan(x):
            return dt.timedelta(seconds=x)
        else:
            return dt.timedelta(0)
        
def days_since_year(day,year):
        if not np.isnan(day):
            return dt.datetime(year,1,1) + dt.timedelta(days=(int(day)-1))
        else:
            return dt.datetime(9999,1,1,0,0,0)

# Nagihara PDS release uses DOY format; this simply breaks it up to datetime equivalent. 
# TODO: rewrite to use strptime and maybe deal with leap seconds.

def nagihara_doy_to_dt(nagihara_time):
    year=(dt.datetime(int(nagihara_time[0:4]),1,1,0,0))
    day_of_year=dt.timedelta(days=int(nagihara_time[5:8])-1)
    hour=dt.timedelta(hours=int(nagihara_time[9:11]))
    minute=dt.timedelta(minutes=int(nagihara_time[12:14]))
    second=dt.timedelta(seconds=float(nagihara_time[15:21]))
    return year+day_of_year+hour+minute+second

# main function for ingesting data from Nagihara et al. 2018

def ingest_nagihara_2018(nagihara_data={},spreadsheet_path='./source/nagihara/jgre.xlsx'):
    nagihara_datafiles = ['a15_1975','a17_1975','a17_1976','a17_1977']
    nagihara_spreadsheet={}
    for f in enumerate(nagihara_datafiles):
        nagihara_spreadsheet[f[1]]=pd.read_excel(spreadsheet_path,f[0],header=0,
            delim_whitespace=True)  

    # The 1975 data from this paper has been superseded by the related July 2019 PDS release.
    # We read in the whole spreadsheet for convenience, but only actually ingest the 1976 and 1977 data.
    # as elsewhere, m is mission (here mission by year), p is probe, 
    # s is 'file,' i.e. bridge identifier, which conveniently here is only ever '1.'

    nagihara_data = {}
    for m in nagihara_spreadsheet.keys():
        if not m[4:]=='1975':
            if not m in nagihara_data.keys():
                nagihara_data[m] = {}
            for p in ['p1','p2']:
                if not p in nagihara_data[m].keys():
                    nagihara_data[m][p] = {}
                for s in [1]:
                    if not s in nagihara_data[m][p].keys():
                        nagihara_data[m][p][s] = {}

                    # The non-deprecated data from this paper is from the upper gradient bridges 
                    # of both Apollo 17 probes.  

                    # we initially read each 'file' (bridge) as a dict of arrays for convenience.

                    # Time

                    # the original HFE dataset reduced time data into an epoch time format as 
                    # described above. Nagihara et al. reduced time into days from the beginning 
                    # of the then-current calendar year and seconds from the beginning of that day. 
                    # This converts Nagihara et al.'s time data into the format given in the original HFE 
                    # dataset.
                    
                    column_indicator='.1' if p=='p2' else ''
                    nagihara_data[m][p][s]['UTC_Time']=\
                       np.vectorize(days_since_year)\
                        (nagihara_spreadsheet[m]['Day'+column_indicator],\
                        int(m[-4:]))+np.vectorize(seconds_interval)\
                        (nagihara_spreadsheet[m]['Seconds'+column_indicator])
                    nagihara_data[m][p][s]['Time']=\
                        np.vectorize(a17_time)\
                        (np.vectorize(utc_to_tai)\
                        (nagihara_data[m][p][s]['UTC_Time']))
                    # The original HFE dataset reduced temperature data into T and dT, where T is the
                    # average temperature across the bridge and dT is the difference between the upper and lower parts
                    # of the bridge. Nagihara et al. reduced the ALSEP data differently, explicitly giving temperature 
                    # values for 'A' (upper) and 'B' (lower) parts of the bridge. This converts Nagihara et al.'s 
                    # temperature data into the format given in the original HFE dataset.
                    
                    TGA = nagihara_spreadsheet[m]['TG'+p[1]+str(s)+'A']
                    TGB = nagihara_spreadsheet[m]['TG'+p[1]+str(s)+'B']
                    
                    # T

                    nagihara_data[m][p][s]['T'] = (TGA + TGB) / 2
                    
                    # dT

                    nagihara_data[m][p][s]['dT'] = TGB - TGA
                        
                    # We also retain Nagihara et al.'s explicitly-computed bridge values 
                    # to avoid rounding errors later.

                    nagihara_data[m][p][s]['TG'+p[1]+str(s)+'A'] = TGA
                    nagihara_data[m][p][s]['TG'+p[1]+str(s)+'B'] = TGB
                
            # converts 'files' to pandas dataframes and drops empty rows.
            # we don't flag empty rows because they're just formatting artifacts.
            # previous versions wrote intermediate csv files, but we don't bother here;
            # may add later if efficiency becomes an issue or they prove desirable
            # for some other reason.

    for m in nagihara_data:
        for p in nagihara_data[m]:
            for s in nagihara_data[m][p]:
                nagihara_data[m][p][s]=pd.DataFrame.from_dict\
                                           (nagihara_data[m][p][s]).dropna()
    return nagihara_data

def ingest_nagihara_2019(nagihara_data={},pathname='.'):

    # grab files in directory according to Nagihara et al.'s naming convention
    # and read them as pandas dataframes

    # These files are parsable as csv with fields separated by a variable number of spaces.

    # as elsewhere, m is mission (here mission by year), p is probe number, 
    # s is 'file,' i.e. bridge identifier.

    for m in ['a15_1975','a17_1975']:
        nagihara_data[m]={}
        for p in ['p1','p2']:
            nagihara_data[m][p]={}
            for s in [1,2]:
                filename='{pathname}/{m}_hfe_1975_l2_arcsav_tg{p}{s}.tab'\
                         .format(pathname=pathname,m=m[0:3],p=p[1],s=s)    
                if os.path.isfile(filename):
                    nagihara_data[m][p][s]=pd.read_csv(filename,engine='python',\
                                                   sep=r'\ +')
                    
                    # Time
                    
                    # generate a separate column for mission epoch time.
                    # convert native DOY format (e.g. '1975-092T00:04:00.817')
                    # to datetime as an intermediate step.
                    
                    # to_pydatetime is necessary because pandas otherwise gets mad
                    # about loss of (here nonexistent) nanosecond precision.
                    
                    # then remove leap seconds for parity with NSSDC data and convert to mission epoch time. 

                    nagihara_data[m][p][s]['UTC_Time']=\
                        np.vectorize(nagihara_doy_to_dt)(nagihara_data[m][p][s]['time'])
                                                
                    if m[0:3] == 'a17':
                        nagihara_data[m][p][s]['Time']=\
                            np.vectorize(a17_time)\
                            (np.vectorize(tai_to_utc)\
                            (nagihara_data[m][p][s]['UTC_Time'].dt.to_pydatetime()))

                    elif m[0:3] == 'a15':
                        nagihara_data[m][p][s]['Time']=\
                            np.vectorize(a17_time)\
                            (np.vectorize(tai_to_utc)\
                            (nagihara_data[m][p][s]['UTC_Time'].dt.to_pydatetime()))
                    nagihara_data[m][p][s].drop(columns=['time'],inplace=True)
                    
                    # dT 
                    
                    # Nagihara et al. include both the high- and low-sensitivity dT 
                    # measurements for every data point. To construct our reduced sets, 
                    # we choose one of these measurements per point as a canonical dT.
                    # For each point, we simply select the dT measurement that Nagihara 
                    # et al. used to calculate their canonical bridge temperatures.
                    
                    n=nagihara_data[m][p][s]['Time'].size
                    nagihara_data[m][p][s]['dT']=pd.Series(np.zeros(n))
                    
                    dT = nagihara_data[m][p][s]['dT']
                    TA = nagihara_data[m][p][s]['TG{p}{s}A'.format(p=p[1],s=s)]
                    TB = nagihara_data[m][p][s]['TG{p}{s}B'.format(p=p[1],s=s)]
                    DTH = nagihara_data[m][p][s]['DTH{p}{s}'.format(p=p[1],s=s)]
                    DTL = nagihara_data[m][p][s]['DTL{p}{s}'.format(p=p[1],s=s)]
                    
                    # is TA - TB within a rounding error of a given measurement?
                    # then select that measurement as canonical dT, preferring DTH
                    # if both measurements qualify.
                    
                    DTH_index=round(abs(TA-TB-DTH),3) <= 0.01
                    DTL_index=round(abs(TA-TB-DTL),3) <= 0.01
                    
                    nagihara_data[m][p][s].loc[DTL_index,'dT']=DTL.loc[DTL_index]
                    nagihara_data[m][p][s].loc[DTH_index,'dT']=DTH.loc[DTH_index]
                    
                    # are there points where neither measurement qualifies? something has
                    # gone wrong with the analysis.
    
                    if not np.all(np.bitwise_or(DTL_index,DTH_index)):
                        raise ValueError('Margin of error for PDS data dT selection appears to be off.')
                    
                    # then flip the sign for parity with the NSSDC convention

                    nagihara_data[m][p][s]['dT']=nagihara_data[m][p][s]['dT']*-1

                    # rename and reindex columns for parity with NSSDC data and later convenience
                    
                    nagihara_data[m][p][s].\
                        rename(columns={'TG{p}{s}avg'.format(p=p[1],s=s): 'T'},inplace=True)
                    nagihara_data[m][p][s]=nagihara_data[m][p][s].\
                        reindex(columns=['Time','T','dT','UTC_Time','TG{p}{s}A'.\
                        format(p=p[1],s=s),'TG{p}{s}B'.format(p=p[1],s=s)])
        
    return nagihara_data

# functions for polishing dataset and splitting to thermometers.

def discard_flagged(data):
    data_clean=pd.DataFrame(columns=data.columns)
    # retain ambiguous dT corrections and time inversions
    forbidden_flags=sum([0b1,0b100,0b1000,0b10000,0b100000,0b1000000,0b10000000,
                         0b100000000,0b1000000000])
    ix=(np.bitwise_and(data['flags'],forbidden_flags).values==0)
    for c in data.columns:
        data_clean[c]=data[c].loc[ix].values
    return data_clean

def substitute_corrections(data):
    if 'dT_corr' in data.columns:
        data['dT']=data['dT_corr']
        data.drop(columns=['dT_corr'],inplace=True)

def standardize_time(data,mission):

        # use directly-reformatted time from Nagihara; it's given to millisecond precision and 
        # meaningful floating-point errors could plausibly be introduced.

        # otherwise go ahead and do the full conversion from epoch time.

    if 'UTC_Time' in data.columns:
        data['Time']=data['UTC_Time']
        data.drop(columns=['UTC_Time'],inplace=True)
    elif mission[0:3]=='a17':
        data['Time']=np.vectorize(a17_to_greg)(data['Time'])
        data['Time']=np.vectorize(tai_to_utc)(data['Time'])
    elif mission[0:3]=='a15':
        data['Time']=np.vectorize(a15_to_greg)(data['Time'])
        data['Time']=np.vectorize(tai_to_utc)(data['Time'])

def polish_hfe (data):
    for m in data.keys():
        for p in data[m].keys():
            for s in data[m][p].keys():
                data[m][p][s]=discard_flagged(data[m][p][s])
                substitute_corrections(data[m][p][s])
                standardize_time(data[m][p][s],m)

def join_hfe_years(data):
    for m in ['a15_1975','a17_1975','a17_1976','a17_1977']:
        for p in data[m]:
            for s in data[m][p]:
                data[m[0:3]][p][s]=data[m[0:3]][p][s].append\
                    (data[m][p][s]).reset_index(drop=True)
        del data[m]

def split_to_thermometers(data):
    data_split={}
    for m in data.keys():
        data_split[m]={}
        for p in data[m].keys():
            data_split[m][p]={}
            for s in data[m][p].keys():
                if m[3:]=='':
                    # construct per-thermometer temperature fields for NSSDC data
                    probe=p[1]
                    if s==1:
                        sensor=['G','1']
                    elif s==2:
                        sensor=['G','2']
                    elif s==4:
                        sensor=['R','1']
                    elif s==5:
                        sensor=['R','2']
                    if s!=3:
                        # calculate per-thermometer temperatures from NSSDC data
                        A='T'+sensor[0]+p[1]+sensor[1]+'A'
                        B='T'+sensor[0]+p[1]+sensor[1]+'B'
                        data_split[m][p][s]=pd.DataFrame(columns=['Time',A,B])
                        data_split[m][p][s]['Time']=data[m][p][s]['Time']
                        data_split[m][p][s][A]=data[m][p][s]['T']-0.5*data[m][p][s]['dT']
                        data_split[m][p][s][B]=data[m][p][s]['T']+0.5*data[m][p][s]['dT']
                        data_split[m][p][s]['flags']=data[m][p][s]['flags']
                    elif s==3:
                        TC='TC{p}'.format(p=probe)
                        data_split[m][p][s]=data[m][p][s].drop('HTR',axis=1)
                        data_split[m][p][s]=data_split[m][p][s].rename(columns={'TC1': TC+'1',
                        'TC2': TC+'2', 'TC3': TC+'3', 'TC4': TC+'4', 'TREF': 'TREF', 
                        'flags': 'flags'})
                    for column in data_split[m][p][s]:
                    # round time back to second precision rather than the microsecond precision introduced
                    # by astropy time scale conversion. 
                        if column == 'Time':
                            data_split[m][p][s][column]=data[m][p][s][column].dt.round('1s')
                    # retains the (probably spurious) 5 digits after decimal given by the Lamont data, rather than the
                    # (definitely spurious) additional digits of precision introduced by Python's floating point representation.
                        elif column[0] == 'T':
                            data_split[m][p][s][column]=round(data_split[m][p][s][column],5)
                else:
                    # simply use temperatures calculated by Nagihara et al. to avoid introducing errors
                    for s in data[m][p].keys():
                        data_split[m][p][s]=data[m][p][s].drop(columns=['T','dT'])
                        for column in data_split[m][p][s]:
                            # retain millikelvin precision of Nagihara et al. sets rather than
                            # spurious precision introduced by numpy floating point representation
                            if column[0:2]=='TG':
                                data_split[m][p][s][column]=round(data_split[m][p][s][column],3)
    return data_split

# functions & depth dictionary for writing combined 'depth' set

# cm below surface for each sensor. 'S' marks sensors lying somewhere on the surface, 
# 'N' marks the fourth thermocouple at the same position as the top gradient thermometer,
# 'X' marks off-scale bridges. the assembly functions for the third / depth set
# exclude these ranges from consideration; they are included for possible future work.

depthdict = {
    'a17':{
     'p1':
      {1:
       {'TG11A':130,'TG11B':177},
       2:
       {'TG12A':185,'TG12B':233},
       3:
       {'TC11':'S','TC12':14,'TC13':66,'TC14':'N','TREF':'S'},
       4:
       {'TR11A':140,'TR11B':167},
       5:
       {'TR12A':195,'TR12B':223}
      },
    'p2':{
       1:
       {'TG21A':131,'TG21B':178},
       2:
       {'TG22A':186,'TG22B':234},
       3:
       {'TC21':'S','TC22':15,'TC23':67,'TC24':'N','TREF':'S'},
       4:
       {'TR21A':140,'TR21B':169},
       5:
       {'TR22A':196,'TR22B':224}
      }
     },
'a15':{
     'p1':
      {1:
       {'TG11A':35,'TG11B':84},
       2:
       {'TG12A':91,'TG12B':139},
       3:
       {'TC11':'S','TC12':'S','TC13':0,'TC14':'N','TREF':'S'},
       4:
       {'TR11A':45,'TR11B':73},
       5:
       {'TR12A':101,'TR12B':129}
      },
     'p2':
     {1:
       {'TG21A':[-6,'X'],'TG21B':[32,'X']},
       2:
       {'TG22A':49,'TG22B':97},
       3:
       {'TC21':'S','TC22':'S','TC23':'S','TC24':0,'TREF':'S'},
       4:
       {'TR21A':[3,'X'],'TR21B':[41,'X']},
       5:
       {'TR22A':59,'TR22B':87}
     }
}
}

# concatenates all subsurface sensors per probe (barring TC4) and assigns depth values.

def combine_with_depth(data):
    data_clean={}
    for m in ['a15','a17']:
        data_clean[m]={}
        for p in ['p1','p2']:
            data_clean[m][p]=pd.DataFrame(columns=['Time','T','sensor','depth','flags'])
    for m in data.keys():
        for p in data[m].keys():
            for s in data[m][p].keys():
                for c in data[m][p][s].columns:
                    # is it a temperature value from a sensor we want to include in this set?
                    if c[0]=='T' and c[1]!='i': 
                        if type(depthdict[m][p][s][c])==int and depthdict[m][p][s][c] > 0: 
                            depth_slice=data[m][p][s][['Time',c,'flags']]
                            depth_slice=depth_slice.reindex(['Time','T',c,'sensor','depth','flags'],
                                                            axis=1)
                            depth_slice['depth']=depthdict[m][p][s][c]
                            depth_slice['sensor']=c
                            depth_slice['T']=depth_slice[c]
                            depth_slice=depth_slice.drop(c,axis=1)
                            data_clean[m][p]=data_clean[m][p].append(depth_slice)
            data_clean[m][p]=data_clean[m][p].sort_values(by='Time').reset_index(drop=True)
    return data_clean
