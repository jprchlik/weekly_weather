__version__ = "0.1.0 (2017/09/29"
__authors__ = ["Jakub Prchlik <jakub.prchlik@cfa.harvard.edu>"]
__email__   = "jakub.prchlik@cfa.harvard.edu"

from sunpy.net import hek
import numpy as np
import pandas as pd

def format_string(st):
    """
    formats endtime string into a format of YYYYMMDD_HHMMSS

    Parameters
    ----------
    st: string
        A datetime string in the format YYYY/MM/DD HH:MM

    Returns
    -------
    out: string
        A string in the form of YYYYMMDD_HHMM
    """
    out = st.replace('/','').replace(':','').replace(' ','_')
    return out



def flare_query(tstart,tend,odir=''):
    """
    A program that queries the number of flares in a given time frame, which writes the output into
    a file formatted YYYYMM_HHMM based on tstart into odir.

    Parameters
    ----------
    tstart : string 
        A datetime string in any format accepted by hek.attrs.Time which specifies 
        when to start looking for flares in the HEK.
    tend   : string 
        A datetime string in any format accepted by hek.attrs.Time which specifies 
        when to stop looking for flares in the HEK.
    odir   : string
        Output directory (default = Current directory)

    """

    

    #create base pandas Dataframe
    f_df = pd.DataFrame(columns=['AR','time','time_dt','goes','score'])
    
    client = hek.HEKClient()
    
    #flare event type
    event_type = 'FL'
    
    #query the HEK
    result = client.search(hek.attrs.Time(tstart,tend),hek.attrs.EventType(event_type),(hek.attrs.OBS.Instrument == 'GOES'))
    
    
    #turn dictionary into pandas Data frame
    for j,i in enumerate(result): 
       
        f_df.loc[j] = [i['ar_noaanum'],i['event_peaktime'],pd.to_datetime(i['event_peaktime']),i['fl_goescls'],i['event_score']]
        
    
    #change index to input time
    f_df.set_index(f_df['time_dt'],inplace=True)
    
    #Levels of flares
    f_l = ['X','M','C','B']
    
    #Total number of flares
    f_df['T'] = 0
    
    #make columns for is flare classification 
    for i in f_l: 
        f_df[i] = 0
        #set values to 1 where first character matches value
        f_df[i][f_df.goes.str[0] == i] = 1
        f_df['T'][f_df.goes.str[0] == i] = 1
    
    #day resampling
    d_df = f_df.resample('1D').sum()


    #get index for top percentile rank
    top_ten = int(len(result)*.1)
    if top_ten == 0: top_ten = 1

    #best flares of the week according to HEK event score
    b_df = f_df.sort_values('score',ascending=False)[0:top_ten]
    
    
    #noaa number resampling 
    n_df = f_df.groupby(['AR']).sum().reset_index()
    n_df.set_index(n_df['AR'],inplace=True) # = f_df.groupby(['AR']).sum().reset_index()
    
    
    
    
    #write output to file
    tout = format_string(tend)
    out_f = open(odir+'flares_'+tout+'.txt','w')
    
    #write the output broken up by date
    out_f.write('#######DATE BREAKDOWN##############\n')
    out_f.write(d_df[['X','M','C','B','T']].to_string()+'\n\n')
    #write output broken up by AR
    out_f.write('#########AR BREAKDOWN##############\n')
    out_f.write(n_df[['X','M','C','B','T']].to_string()+'\n\n')
    #top ten percent of flare events
    out_f.write('#########TOP TEN###################\n')
    out_f.write(b_df[['AR','goes','score']].to_string()+'\n\n')
    out_f.close()