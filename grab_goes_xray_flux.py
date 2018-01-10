import ftplib
import numpy as np
import os
import sys, getopt
from datetime import timedelta,datetime




#Funtion to retrieve events from given day on noao ftp server
def getfiles(day,ftp,sdir):
    files = '{0:%Y%m%d}_Gs_xr_1m.txt'.format(day)
    print files
#create file to write ftp data
    fhandle = open(sdir+'/goes/'+files,'wb')
    try:
        ftp.retrbinary('RETR {0}'.format(files),fhandle.write)
    except:
        print '{0} not in archive'.format(files)

    fhandle.close()

#Funtion to retrieve events from given day on noao ftp server
def getfiles_ace(day,ftp,sdir,swepam=False,mag=False):
    if swepam: files = '{0:%Y%m%d}_ace_swepam_1m.txt'.format(day)
    if mag: files = '{0:%Y%m%d}_ace_mag_1m.txt'.format(day)
#create file to write ftp data
    fhandle = open(sdir+'/ace/'+files,'wb')
    try:
        ftp.retrbinary('RETR {0}'.format(files),fhandle.write)
    except:
        print '{0} not in archive'.format(files)

    fhandle.close()

#Function to retrieve DSCOVR events from swpc json server
def getfiles_dscvor(sdir,swepar=False,mag=False,
    path='http://services.swpc.noaa.gov/products/solar-wind/'):


    import pandas as pd
    #Read json objects from swpc
    ppar = pd.read_json(path+'plasma-7-day.json')
    mpar = pd.read_json(path+'mag-7-day.json')

    #set columns
    ppar.columns = ['time_tag','density','speed','temperature']
    mpar.columns = ["time_tag","bx_gsm","by_gsm","bz_gsm","lon_gsm","lat_gsm","bt"]

    #[Create list of object to use for time formatting
    obj = [ppar,mpar]
    for k in obj:
        k.drop(k.index[0],inplace=True)
        k['time_dt'] = pd.to_datetime(k["time_tag"])
        k.set_index(k.time_dt,inplace=True)
        k['Julian'] = k.index.to_julian_date()
        k['YR'] = k["time_tag"].astype(str).str[0:4]
        k['MO'] = k["time_tag"].astype(str).str[5:7]
        k['DA'] = k["time_tag"].astype(str).str[8:10]
        k['HH'] = k["time_tag"].astype(str).str[11:13]
        k['MM'] = k["time_tag"].astype(str).str[14:16]
        k['HHMM'] = k['HH']+k['MM']
        k['Sec'] = k['HH'].astype(int)*3600+k['MM'].astype(int)*60
        k['S'] = 0
    
    #cut down on columns to only ones required
    ppar = ppar[['YR','MO','DA','HHMM','Julian','Sec','S','density','speed','temperature']]
    mpar = mpar[['YR','MO','DA','HHMM','Julian','Sec','S','bx_gsm','by_gsm','bz_gsm','bt','lon_gsm','lat_gsm']]

    #set name of output plasma and mag file
    outpls = sdir+'/ace/dscovr_swe.csv'
    outmag = sdir+'/ace/dscovr_mag.csv'

    #If plasma already exists add to file
    if os.path.isfile(outpls):
        oppar = pd.read_csv(outpls)
        ppar = pd.concat([oppar,ppar])
        ppar.drop_duplicates(['YR','MO','DA','HHMM'],inplace=True) 
    #If magnetic field already exists add to file
    if os.path.isfile(outmag):
        ompar = pd.read_csv(outmag)
        mpar = pd.concat([ompar,mpar])
        mpar.drop_duplicates(['YR','MO','DA','HHMM'],inplace=True) 


    ppar.to_csv(outpls,index=False)
    mpar.to_csv(outmag,index=False)

def look_xrays(start,end,sdir):
#initialize variables
    ftp = ftplib.FTP('ftp.swpc.noaa.gov','anonymous','jakub.prchlik@cfa.harvard.edu')
#change ftp directory to events directory for a given year
    ftp.cwd('/pub/lists/xray/')
     
    #Get and format DSCOVR files  
    getfiles_dscvor(sdir)

    days =  np.arange(start,end,timedelta(days=1)).astype(datetime)
#loop through days     
    for day in days:
#request file from ftp server to add to local directory
        getfiles(day,ftp,sdir)
#change ftp directory to events directory for a given year
    ftp.cwd('../ace/')
    for day in days:
        getfiles_ace(day,ftp,sdir,swepam=True)
        getfiles_ace(day,ftp,sdir,mag=True)
    
#nicely leave the ftp server
    ftp.quit()

