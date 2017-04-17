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

def look_xrays(start,end,sdir):
#initialize variables
    ftp = ftplib.FTP('ftp.swpc.noaa.gov','anonymous','jakub.prchlik@cfa.harvard.edu')
#change ftp directory to events directory for a given year
    ftp.cwd('/pub/lists/xray/')
    
    days =  np.arange(start,end,timedelta(days=1)).astype(datetime)
#loop through days     
    for day in days:
#request file from ftp server to add to local directory
        getfiles(day,ftp,sdir)
    
#nicely leave the ftp server
    ftp.quit()

