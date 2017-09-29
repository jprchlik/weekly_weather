from sunpy.net import hek
import numpy as np
import pandas as pd

def format_string(st):

    out = st.replace('/','').replace(':','').replace(' ','_')
    return out



def flare_query(tstart,tend,odir=''):
     """
     A program that queries the number of flares in a given time frame, which writes the output into
     a file formatted YYYYMM_HHMM based on tstart,




    """

    

    #create base pandas Dataframe
    f_df = pd.DataFrame(columns=['AR','time','time_dt','goes'])
    
    client = hek.HEKClient()
    
    event_type = 'FL'
    
    
    
    #query the HEK
    result = client.search(hek.attrs.Time(tstart,tend),hek.attrs.EventType(event_type),((hek.attrs.AR.NOAANum == '12673') & (hek.attrs.OBS.Instrument == 'GOES')))
    
    
    #turn dictionary into pandas Data frame
    for j,i in enumerate(result): 
        f_df.loc[j] = [i['ar_noaanum'],i['event_peaktime'],pd.to_datetime(i['event_peaktime']),i['fl_goescls']]
        #print i['event_peaktime'],i['fl_goescls']
    
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
    
    
    #noaa number resampling 
    n_df = f_df.groupby(['AR']).sum().reset_index()
    n_df.set_index(n_df['AR'],inplace=True) # = f_df.groupby(['AR']).sum().reset_index()
    
    
    
    
    #write output to file
    tout = format_string(tstart)
    out_f = open('flares_'+tout+'.txt','w')
    
    out_f.write('#######DATE BREAKDOWN##############\n')
    out_f.write(d_df[['X','M','C','B','T']].to_string()+'\n')
    out_f.write('#########AR BREAKDOWN##############\n')
    out_f.write(n_df[['X','M','C','B','T']].to_string()+'\n')
    #set datetime to be the index
    out_f.close()