import urllib
from datetime import datetime,timedelta
import os
import itertools
from multiprocessing import Pool






class download:

    def __init__(self,stime,etime,caden,b_dir,syn_arch='http://jsoc.stanford.edu/data/aia/synoptic/',d_wav=[193],
                 w_fmt='{0:04d}.fits',nproc=8):
    
        """
    
        Parameters:
        -----------
        stime: datetime object
            The time to start downloading AIA files.
      
        etime: datetime object
            The time to end downloading AIA files.
        
        caden: datetime time delta object 
            The time cadence to download AIA files.
        
        b_dir: string
            The base directory for locally storing the AIA archive.
        
        syn_arch: string, optional
            Location of online syntopic archive (default = 'http://jsoc.stanford.edu/data/aia/synoptic/').
        
        d_wav: list, optional
            List of wavelengths to download from the online archive (Default = [193]).
        
        w_fmt: list, optional
            The wavelength format of the wavelength list (defatult = '{0:04d}.fits').
        
        f_dir : string, optional
            Create local directory path format for SDO/AIA files (default = '{0:%Y/%m/%d/H%H00/AIA%Y%m%d_%H%M_}')
        """
    
        self.stime = stime
        self.etime = etime
        self.caden = caden
        self.b_dir = b_dir
        self.syn_arch = syn_arch
        self.d_wav = d_wav
        self.w_fmt = w_fmt
        self.nproc = nproc
    
        
        #desired cadence for the observations
        real_cad = [result for result in self.des_cad(stime,etime,caden)]
        
        #create a list of combination of dates and wavelengths
        inpt_itr = list(itertools.product(real_cad,d_wav))
        
        #Download the files locally in parallel if nproc greater than 1
        if self.nproc < 2:
            for i in intp_itr: self.wrap_download_file(i)
        else:
            pool = Pool(processes=self.nproc)
            outp = pool.map(self.wrap_download_file,inpt_itr)
            pool.close()
            pool.join()

    #retrieve desired cadence from file list
    def des_cad(self,start,end,delta):
        """Create an array from start to end with desired cadence"""
        curr = start
        while curr < end:
            yield curr
            curr += delta
    
    #wrapper for download file for par. processing
    def wrap_download_file(args):
        return self.download_file(*args)
    
    #download files from archive for each wavelength
    def download_file(self,time,wavl):
  
       #Init variables
       w_fmt = self.w_fmt
       f_dir = self.f_dir
       b_dir = self.b_dir

       #format wavelength
       w_fil = w_fmt.format(wavl)
       #format input time
       s_dir = f_dir.format(time)
       #create output file
       o_fil = b_dir+s_dir.split('/')[-1]+w_fil
       #file to download from archive
       d_fil = syn_arch+s_dir+w_fil
    
       #check if output file exists
       if os.path.isfile(o_fil) == False:
           #try to download file if fails continue on
           try:
               urllib.urlretrieve(d_fil,o_fil) 
           except:
               print("Cound not Download {0} from archive".format(d_fil))
    
