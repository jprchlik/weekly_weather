Weekly Solar Weather
=========================

Summary
-------
One job of a solar support scientist is making a video which spans the last 7 days of solar activity observed by SDO.
This weekly_weather python program fulfills this duty without the need to locate, download, and compile the movie yourself by hand,
which frees up your time to do more exciting things.
All you need to do is schedule a cronjob to run run_weekly_weather.csh and you should not have to worry about making a week long movie every week.

SMEARpy.py
===========
Code originally created by J. Sattleberger to locate and verify AIA files. It is now updated to be python2 and python3 compatible.

Videowall_Notes.docx
====================
A Word Document containing information on how to run the Video Wall. A printed copy exists at the Video Wall station.

flare_query.py
==============
A program that queries the number of flares in a given time frame and writes the output into formatted file.

get_syntoptic_files.py
==============
Contains the class that downloads SDO/AIA synoptic files from jsoc http server.

grab_goes_xray_flux.py
==============
Downloads GOES flux from NOAA ftp directory.

run_weekly_weather.py
==============
Creates unprocessed weekly solar movie in 193 Angstroms with indicators of solar activity on the right.

run_weekly_weather_wavelet.py
==============
Creates a wavelet processed solar movie in 193 Angstroms with indicators of solar activity on the right.
