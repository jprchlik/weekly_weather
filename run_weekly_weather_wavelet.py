import matplotlib

#fixes multiprocess issue
matplotlib.use('agg')
import sunpy.map
from matplotlib.transforms import Bbox
import sunpy.visualization.colormaps as cm
import matplotlib.dates as mdates
import subprocess
from PIL import Image
import glob
import os
import stat
import numpy as np
from sunpy.net.helioviewer import HelioviewerClient
from datetime import date,datetime
from datetime import timedelta as dt
from multiprocessing import Pool
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import grab_goes_xray_flux as ggxf
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from astropy.io import ascii
from astropy.table import vstack,Table,join
import flare_query as fq
import shutil

from SMEARpy import Scream


################ADDED undecimated isotropic wavelet transform####################
#J. Prchlik 2018/01/16
import pywt
#Interested in the Stationary Wavelet Transformation or a Trous or starlet (swt2)

#get a median filter for wavelet transformation
from scipy.signal import medfilt

#Use skimage to create the background mean filter
#from skimage.morphology import disk #pixel radius
#from skimage.filters import rank #rank filters


#Now need to manually download AIA syntoptic files for weekly weather
#2018/04/23 J. Prchlik
import get_synoptic_files as gsf



def get_file(ind):
#Source ID 11 is AIA 193
#    fils = ind.replace(':','_').replace('/','_').replace(' ','_')+'_AIA_193.jp2'
    meta = hv.get_closest_image(ind,sourceId=11)
    date = meta['date'].strftime('%Y_%m_%d__%H_%M_%S')
    mils = float(meta['date'].strftime('%f'))/1000.
    fstr = date+'_*__SDO_AIA_AIA_193.jp2'.format(mils).replace(' ','0')
#    test = os.path.isfile(sdir+'/raw/'+fstr)
    test = glob.glob(sdir+'/raw/'+fstr)
    if len(test) == 0:
#        filep = hv.download_png(ind,0.3,"[11,1,100]",directory=sdir+'/raw',x0=0,y0=0,width=int(rv*8192),height=8192,watermark=False)#,y1=-1200,y2=1200,x1=-1900,x2=1900)
        filep = hv.download_jp2(ind,sourceId="11",directory=sdir+'/raw',clobber=True)#,y1=-1200,y2=1200,x1=-1900,x2=1900)

#J. Prchlik 2016/10/11
#Added to give physical coordinates
def img_extent(img):
# get the image coordinates in pixels
    px0 = img.meta['crpix1']
    py0 = img.meta['crpix2']
# get the image coordinates in arcsec 
    ax0 = img.meta['crval1']
    ay0 = img.meta['crval2']
# get the image scale in arcsec 
    axd = img.meta['cdelt1']
    ayd = img.meta['cdelt2']
#get the number of pixels
    tx,ty = img.data.shape
#get the max and min x and y values
    minx,maxx = px0-tx,tx-px0
    miny,maxy = py0-ty,ty-py0
#convert to arcsec
    maxx,minx = maxx*axd,minx*axd
    maxy,miny = maxy*ayd,miny*ayd

    return maxx,minx,maxy,miny




#for j,i in enumerate(dayarray):
#reformat file to be in 1900x1200 array and contain timetext
def format_img(i):
        global goes,goesdat,sday,eday
        global aceadat,ace

##    try:
        filep = dayarray[i]
    #output file
    #    outfi = sdir+'/working/seq{0:4d}.png'.format(i).replace(' ','0')
        outfi = filep.replace('raw','working').replace('fits','png')
        test = os.path.isfile(outfi)
    	
    #check image quality
    
        check, img = qual_check(filep)
    
    #test to see if bmpfile exists
        if ((test == False) & (check)):
            print('Modifying file '+filep)
            img = sunpy.map.Map(filep)
            fig,ax = plt.subplots(figsize=(sc*float(w0)/float(dpi),sc*float(h0)/float(dpi)))
            fig.set_dpi(dpi)
            fig.subplots_adjust(left=0,bottom=0,right=1,top=1)
            ax.set_axis_off()
            #ax.imshow(img.data,interpolation='none',cmap=cm.sdoaia193,vmin=0,vmax=255,origin='lower')
    		# J. Prchlik 2016/10/06
            #Modified for fits files 
    #        ax.imshow(np.arcsinh(img.data),interpolation='none',cmap=cm.sdoaia193,vmin=np.arcsinh(70.),vmax=np.arcsinh(7500.),origin='lower')
    #Block add J. Prchlik (2016/10/06) to give physical coordinate values 
    #return extent of image
            maxx,minx,maxy,miny = img_extent(img)

            #Add wavelet transformation J. Prchlik 2018/01/17
            wav_img = img.data
            d_size = 15
            #get median filter
            n_back = medfilt(wav_img,kernel_size=d_size)
            #subtract median filter
            img_sub = wav_img-n_back

            #Use Biorthogonal Wavelet
            wavelet = 'bior2.6'
            
            #use 6 levels
            n_lev = 6
            o_wav = pywt.swt2(img_sub, wavelet, level=n_lev )
            #only use the first 4
            #f_img = pywt.iswt2(o_wav[0:4],wavelet)
            #use last four because PyWavelet switched the ordering
            f_img = pywt.iswt2(o_wav[-4:],wavelet)
            #Add  wavelet back into image
            f_img = f_img+wav_img

            #remove zero values
            f_img[f_img < 0.] = 0. 

    #plot the image in matplotlib
    #        ax.imshow(img.data,interpolation='none',cmap=cm.sdoaia193,vmin=0,vmax=255,origin='lower',extent=[minx,maxx,miny,maxy])
            #ax.imshow(np.arcsinh(f_img),interpolation='none',cmap=cm.sdoaia193,origin='lower',vmin=np.arcsinh(5.),vmax=np.arcsinh(7500.),extent=[minx,maxx,miny,maxy])
            #switch to 0.25 power
            #increase v_min to 33.
            ax.imshow((f_img)**0.25,interpolation='none',cmap=cm.sdoaia193,origin='lower',vmin=(15.)**0.25,vmax=(3500.)**0.25,extent=[minx,maxx,miny,maxy])
    #        ax.set_axis_bgcolor('black')
            ax.text(-2000,-1100,'AIA 193 - '+img.date.strftime('%Y/%m/%d - %H:%M:%S')+'Z',color='white',fontsize=36,zorder=50,fontweight='bold')
            if goes:
#format string for date on xaxis
                myFmt = mdates.DateFormatter('%m/%d')

#only use goes data upto observed time

                use, = np.where((goesdat['time_dt'] < img.date+dt(minutes=150)) & (goesdat['Long'] > 0.0))
                clos,= np.where((goesdat['time_dt'] < img.date) & (goesdat['Long'] > 0.0))
                ingoes = inset_axes(ax,width="27%",height="20%",loc=7,borderpad=-27) #hack so it is outside normal boarders
                ingoes.set_position(Bbox([[0.525,0.51],[1.5,1.48]]))
                ingoes.set_facecolor('black')
#set inset plotting information to be white
                ingoes.tick_params(axis='both',colors='white')
                ingoes.spines['top'].set_color('white')
                ingoes.spines['bottom'].set_color('white')
                ingoes.spines['right'].set_color('white')
                ingoes.spines['left'].set_color('white')
#make grid
                ingoes.grid(color='gray',ls='dashdot')

                ingoes.xaxis.set_major_formatter(myFmt)

                ingoes.set_ylim([1.E-9,1.E-2])
                ingoes.set_xlim([sday,eday])
                ingoes.set_ylabel('X-ray Flux (1-8$\mathrm{\AA}$) [Watts m$^{-2}$]',color='white')
                ingoes.set_xlabel('Universal Time',color='white')
                try:
                    ingoes.plot(goesdat['time_dt'][use],goesdat['Long'][use],color='white')
                    ingoes.scatter(goesdat['time_dt'][clos][-1],goesdat['Long'][clos][-1],color='red',s=10,zorder=1000)
                except:
                    print('No GOES data')
                ingoes.set_yscale('log')
#plot ace information
            if ((ace) & (goes)):
                use, = np.where((aceadat['time_dt'] < img.date+dt(minutes=150)) & (aceadat['S_1'] == 0.0) & (aceadat['S_2'] == 0) & (aceadat['Speed'] > -1000.))
                clos,= np.where((aceadat['time_dt'] < img.date) & (aceadat['S_1'] ==  0) & (aceadat['S_2'] == 0) & (aceadat['Speed'] > -1000))
                
                acetop = inset_axes(ingoes,width='100%',height='100%',loc=9,borderpad=-27)
                acebot = inset_axes(ingoes,width='100%',height='100%',loc=8,borderpad=-27)

#set inset plotting information to be white
                acetop.tick_params(axis='both',colors='white')
                acetop.spines['top'].set_color('white')
                acetop.spines['bottom'].set_color('white')
                acetop.spines['right'].set_color('white')
                acetop.spines['left'].set_color('white')

#set inset plotting information to be white
                acebot.tick_params(axis='both',colors='white')
                acebot.spines['top'].set_color('white')
                acebot.spines['bottom'].set_color('white')
                acebot.spines['right'].set_color('white')
                acebot.spines['left'].set_color('white')
#make grid
                acebot.grid(color='gray',ls='dashdot')
                acetop.grid(color='gray',ls='dashdot')


                acetop.set_facecolor('black')
                acebot.set_facecolor('black')


                acetop.set_ylim([0.,50.])
                acebot.set_ylim([200.,1000.])

                acetop.set_xlim([sday,eday])
                acebot.set_xlim([sday,eday])

                acetop.set_xlabel('Universal Time',color='white')
                acebot.set_xlabel('Universal Time',color='white')
 
                acetop.set_ylabel('B$_\mathrm{T}$ [nT]',color='white')
                acebot.set_ylabel('Wind Speed [km/s]',color='white')

                #skip missing ace solar wind data
                try:
                    acetop.plot(aceadat['time_dt'][use],aceadat['Bt'][use],color='white')
                    acebot.plot(aceadat['time_dt'][use],aceadat['Speed'][use],color='white')
                    
                    acetop.scatter(aceadat['time_dt'][clos][-1],aceadat['Bt'][clos][-1],color='red',s=10,zorder=1000)
                    acebot.scatter(aceadat['time_dt'][clos][-1],aceadat['Speed'][clos][-1],color='red',s=10,zorder=1000)
                except:
                    print('Missing ACE data')
                
                acebot.xaxis.set_major_formatter(myFmt)
                acetop.xaxis.set_major_formatter(myFmt)

                
    ##        ax.set_axis_bgcolor('black')
    #        ax.text(-1000,175,'AIA 193 - '+img.date.strftime('%Y/%m/%d - %H:%M:%S')+'Z',color='white',fontsize=36,zorder=50,fontweight='bold')
            fig.savefig(outfi,edgecolor='black',facecolor='black',dpi=dpi)
            plt.clf()
            plt.close()
##    except:
##        print 'Unable to create {0}'.format(outfi)
        return


#for j,i in enumerate(dayarray):
def qual_check(filep):
#read JPEG2000 file into sunpymap
    img = sunpy.map.Map(filep)
#Level0 quality flag equals 0 (0 means no issues)
    lev0 = img.meta['quallev0'] == 0
#check level1 bitwise keywords (http://jsoc.stanford.edu/doc/keywords/AIA/AIA02840_K_AIA-SDO_FITS_Keyword_Document.pdf)
    lev1_a = ((np.binary_repr(img.meta['quality']) == '1000000000000000000000000000000') & (img.exposure_time.value > 1.85))
    #keep some short exposures as long as they are not too short
    lev1_b = ((img.exposure_time.value > 1.85) & (np.binary_repr(img.meta['quality']) == '1000000000000000000000000000100'))
#check to see if it is a calibration image
#This keyword changed after AIA failure
#    calb = np.binary_repr(img.meta['aiftsid']) == '1010000000000000'
#check that both levels pass and it is not a calibration file
    check = ((lev0) & (lev1_a | lev1_b))# & (calb)) 

    return check,img
    #return lev0,img   


##########################################################
# Phase 1: Set parameters for movie creation             #
##########################################################

#move boolean plotting of ace and goes classifcation to the top
goes = True# overplot goes values
ace = True #overplot ACE wind values

#datasources = hv.get_data_sources()

#Only thing to edit if someone else takes over weekly weather
#running in subdirectory so this step is not need (Prchlik J. 2017/03/03)
#stard = '/Volumes/Pegasus/jprchlik/weekly_weather/'


now = datetime.utcnow()

#day of weekly weather (currently Tuesday)
wday = 1

#dimensions of the image
w0 = 1900
h0 = 1144
#video wall ratio
rv = float(w0)/float(h0)


#hour between aia update
dh = 2

#set dpi
dpi = 300
#scale up the images with increasing dpi
sc = dpi/100

#number of processors for downloading
nproc = 5


#the arcsecond values for the image
x0 = 0
y0 = 0
x1 = -1900
x2 =  1900
y1 = -1200
y2 =  1200

span = 7 #days to run the movie over
cadence = 6. #in minutes
minweek = 10080. #minutes in a week
samples = round(minweek/cadence)#number of dt samples to probe

#if the day is the weekly weather do extra stuff
#difference from Tuesday
dday = (now.weekday()-wday) %7

#put extra padding to be safe but check the same day if the answer is zero then set the span to be 7 days
if dday == 0:  dday = 7

#start day 
sday = now-dt(days=dday)
#set to 12 utc on sday
sday = sday.replace(hour=12,minute=0,second=0)

#end day
eday = sday+dt(days=span)


#remove raw files from previous week if they exist
#Do because I have to download files now 2018/04/23 J. Prchlik
#old_dir = '{0:%Y%m%d}/raw'.format(sday)
#if os.path.exists(old_dir): shutil.rmtree(old_dir)



#create a directory which will contain the raw png files
#sdir = stard+eday.date().strftime('%Y%m%d')
#creating a subdirectory to extra step is not need
#Added wavelet distinctions J. Prchlik 2018/01/17
##########################################################
# Phase 2: Create local directories to for data download #
##########################################################
sdir = eday.date().strftime('%Y%m%d')#+'/wavelet'
try:
    os.mkdir(sdir)
    os.mkdir(sdir+'/raw')
    os.mkdir(sdir+'/working')
    os.mkdir(sdir+'/working/symlinks')
    os.mkdir(sdir+'/final')
    os.mkdir(sdir+'/iris')
    os.mkdir(sdir+'/xrt')
    os.mkdir(sdir+'/goes')
    os.mkdir(sdir+'/ace')
except OSError:
    print('Directories Already Exist. Proceeding to Download')

##########################################################
# Phase 3: Get GOES and Solar wind data from archives    #
##########################################################
#get all days in date time span
if goes: 
    ggxf.look_xrays(sday,now+dt(days=1),sdir)
    goesfil = glob.glob(sdir+'/goes/*txt')
    goesnames = [ 'YR', 'MO', 'DA', 'HHMM', 'JDay', 'Secs', 'Short', 'Long'] 
    goesdat = Table(names=goesnames)

#loop over all day information and add to large array
    for m in goesfil:
        temp = ascii.read(m,guess=True,comment='#',data_start=2,names=goesnames)
        goesdat = vstack([goesdat,temp])

    #create datetime array
    goesdat['time_dt'] = [datetime(int(i['YR']),int(i['MO']),int(i['DA']))+dt(seconds=i['Secs']) for i in goesdat]

#ace = False 
if ace:
    #Updated to DSCOVR parameters using CSV 2018/01/09
    aceb = glob.glob(sdir+'/ace/*mag*csv')
    acep = glob.glob(sdir+'/ace/*swe*csv')

    aceb_names = [ 'YR', 'MO', 'DA', 'HHMM', 'JDay', 'Secs', 'S', 'Bx','By','Bz','Bt','Lat','Long'] 
    acep_names = [ 'YR', 'MO', 'DA', 'HHMM', 'JDay', 'Secs', 'S', 'Den','Speed','Temp'] 
  
    acebdat = Table(names=aceb_names)
    acepdat = Table(names=acep_names)

#put B field in large array
    for m in aceb:
        temp = ascii.read(m,guess=True,comment='#',data_start=2,names=aceb_names)
        acebdat = vstack([acebdat,temp])
#put plasmag in large array
    for m in acep:
        temp = ascii.read(m,guess=True,comment='#',data_start=2,names=acep_names)
        acepdat = vstack([acepdat,temp])


    aceadat = join(acepdat,acebdat,keys=['YR','MO','DA','HHMM'])
    #create datetime array
    aceadat['time_dt'] = [datetime(int(i['YR']),int(i['MO']),int(i['DA']))+dt(seconds=i['Secs_1']) for i in aceadat]

#create a starting time for the weekly movie, which is the previous Tuesday at 12:00:00 utc
nday = sday

#write ffmpeg command
com = open(sdir+'/run_ffmpeg.csh','w')
com.write('/usr/local/bin/ffmpeg -y -f image2 -r 25 -i working/symlinks/seq%4d.png -an -pix_fmt "yuv420p" -vcodec libx264 -level 41 -crf 18.0 -b 8192k -r 25 -bufsize 8192k -maxrate 8192k -g 25 -coder 1 -profile main -preset faster -qdiff 4 -qcomp 0.7 -directpred 3 -flags +loop+mv4 -cmp +chroma -partitions +parti4x4+partp8x8+partb8x8 -subq 7 -me_range 16 -keyint_min 1 -sc_threshold 40 -i_qfactor 0.71 -rc_eq "blurCplx^(1-qComp)" -s "1900x1180" -b_strategy 1 -bidir_refine 1 -refs 6 -deblockalpha 0 -deblockbeta 0 -trellis 1 -x264opts keyint=25:min-keyint=1:bframes=1 -threads 2 final/{0}_weeklyweather.mp4\n'.format(eday.date().strftime('%Y%m%d')))
com.close()





###Added output file recapping the last week of flares
try:
    qfmt = '%Y/%m/%d %H:%M:%S'
    print(sday,eday)
    fq.flare_query(sday.strftime(qfmt),eday.strftime(qfmt),odir=sdir+'/final/')
except:
    print('FLARE SYNPOPSIS DID NOT RUN')

#Added downloading of files
##########################################################
# Phase 4: get AIA file from server                      #
##########################################################
gsf.download(sday,eday,dt(minutes=cadence),sdir+'/raw/',d_wav=[193],nproc=nproc)

#J. Prchlik 2016/10/06
#Updated version calls local files
##########################################################
# Phase 5: get file names                                #
##########################################################
verbose=False
debug = False
#archive = "/data/SDO/AIA/synoptic/"
archive = os.getcwd()+'/'+sdir+'/raw/'
src = Scream(archive=archive,verbose=verbose,debug=debug)
sendspan = "-{0:1.0f}d".format(span) # need to spend current span not total span
paths = src.get_paths(date=eday.strftime("%Y-%m-%d"), time=eday.strftime("%H:%M:%S"),span=sendspan)
fits_files = src.get_filelist(date=eday.strftime("%Y-%m-%d"),time=eday.strftime("%H:%M:%S"),span=sendspan,wavelnth='193')
qfls, qtms = src.run_quality_check(synoptic=True)
fits_files = src.get_sample(files = qfls, sample = '6m', nfiles = 1)
#fits_times = src.get_filetimes(files=fits_files)
#set the file cadence with synoptic 3 minute images being the base
#fits_files = fits_files[::int(cadence)/3]

#create symbolic links for AIA syntopic files in raw directory
for i in fits_files:
    newfile = i.split('/')[-1]
    try:
        os.symlink(i,sdir+'/raw/'+newfile)
    except OSError:
        continue
#        print 'Symlink Already Exists'


##########################################################
# Phase 6: Create png files                              #
##########################################################
#write new file and add text for time
#dayarray = dayarray[0:1]
#dayarray = glob.glob(sdir+'/raw/*jp2')
#J. Prchlik 2016/10/06
#Switched jp2 to fits
dayarray = glob.glob(sdir+'/raw/*fits')
forpool = np.arange(len(dayarray))

#loop is for testing purposes
#for i in forpool: format_img(i)
pool1 = Pool(processes=nproc)
outs = pool1.map(format_img,forpool)
pool1.close()
pool1.join()


startd = sdir+'/' #start from the base directory to create symbolic link
#J. Prchlik 2016/10/06
   # create new symbolic links in order 
fipng = glob.glob(startd+'working/*png')
for i,outfi in enumerate(fipng):
    symli = startd+'/working/symlinks/seq{0:4d}.png'.format(i).replace(' ','0')
    if os.path.islink(symli): os.unlink(symli) # replace existing symbolic link
    os.symlink('../'+outfi.split('/')[-1],symli)

   
##########################################################
# Phase 7: Create movie                                  #
##########################################################
#change to current directory
os.chdir(sdir)

#change file to executable
mod = subprocess.call(['/bin/csh','-c','chmod a+x run_ffmpeg.csh'])

#run ffmpeg
run = subprocess.call(['/bin/csh','-c','./run_ffmpeg.csh'])


