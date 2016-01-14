#!/usr/bin/python
""" vex2zapints.py [options] <observation.vex> <source name> <telescope name>

Examples: 

  vex2zapints.py -f n15dh01a.fil n15dh01a.vex.difx PSRJ1745-290 Ku
  vex2zapints.py -s 57300.005 -t 0.005 n15dh01a.vex.difx PSRJ1745-290 Ku

Reads the VEX file of a pulsar phase referenced VLBI observation, and
outputs a list of time ranges (PRESTO input data 'intervals') during
which the telescope was off source.

The list is written into 'zapints.cmd'. In principle, it should
be ready for use for PRESTO time range flagging via, e.g.,

  rfifind -blocks 1 -zapints <start1>:<stop1>,...,<startN>:<stopN> [-filterbank file.fb]

Source and telescope names should match the names used in the VEX file."""

#############################################################################################

import datetime, math, optparse, subprocess, sys
import vex #  Mark Kettenis' Python VEX parser (http://www.jive.nl/nexpres/doku.php?id=nexpres:nexpres_wp7)

#############################################################################################

HEADER_PROG = 'header'  # location of SIGPROC program called 'header'

#############################################################################################

parser = optparse.OptionParser(usage=__doc__, version='%prog ' + '1.1  (C) 2015 Jan Wagner')
parser.add_option('--filterbank', '-f',
	type='str', dest='filterbankfile', default=None,
	help='Filter bank file for which to check the start time stamp and integration period length. Requires SIGPROC utility called "header".')
parser.add_option('--startmjd', '-s',
	type='float', dest='startmjd', default=None,
	help='The starting MJD (e.g., 57300.005) of the filterbank data, if no filterbank file is specified.')
parser.add_option('--tint', '-t',
	type='float', dest='tint', default=None,
	help='The integration time in seconds of the filterbank data, if no filterbank file is specified.')
parser.add_option('--blocks', '-b',
	type='int', dest='blocks', default=1,
	help='The value to be later used for the -blocks parameter of rfifind.')

#############################################################################################

def dateVEX2split(vdate):
	"""Split VEX date string (e.g., 2015y262d11h56m15s) into its individual parts, and calculate the MJD"""
	d = {}
	d['vex'] = vdate
	d['yr']  = int(vdate[0:4])
	d['doy'] = int(vdate[5:8])
	d['hh']  = int(vdate[9:11])
	d['mm']  = int(vdate[12:14])
	d['ss']  = int(vdate[15:17])
	# Calculate month, day of month
	dt = datetime.date.fromordinal(datetime.date(d['yr'],1,1).toordinal() + d['doy']-1)
	d['mon']  = dt.month
	d['mday'] = dt.day
	# Calculate MJD
	# http://aa.usno.navy.mil/data/docs/JulianDate.php
	# http://www.csgnetwork.com/julianmodifdateconv.html
	jdfract = (d['hh']-12+d['mm']/60.0+d['ss']/3600.0)/24.0
	d['mjd'] = date_to_jd(dt.year,dt.month,dt.day)+0.5 - 2400000.5 + jdfract
	return d

def dateVEX2sec(vdate):
	"""Return VEX date as seconds-of-year"""
	doy = int(vdate[5:8])
	hh  = int(vdate[9:11])
	mm  = int(vdate[12:14])
	ss  = int(vdate[15:17])
	return ss + 60.0*mm + 3600.0*hh + 86400.0*doy

def dateSec2VEX(sec,year):
	"""Convert seconds-of-year into VEX date"""
	doy = int(sec/86400.0)
	R   = sec - doy*86400.0
	hh  = int(R/3600.0)
	R   = R - hh*3600.0
	mm  = int(R/60.0)
	R   = R - mm*60.0
	ss  = int(R)
	return ('%04dy%03dd%02dh%02dm%02ds' % (int(year),doy,hh,mm,ss)) # 2015y262d11h56m15s

def dateVEXaddsec(vdate,sec):
	yr = vdate[0:4]
	s = dateVEX2sec(vdate)
	s = s + sec
	v = dateSec2VEX(s,yr)
	return v

def doFlag(vex_tstart, vex_tstop, fb_tstart, fb_tint, zapints_list):
	t1 = dateVEX2split(vex_tstart)
	t2 = dateVEX2split(vex_tstop)
	flag_startint = max(0, math.floor((t1['mjd']-fb_Tstart)/fb_Tint))
	flag_stopint  = max(0, math.ceil((t2['mjd']-fb_Tstart)/fb_Tint))
	if flag_startint==flag_stopint:
		pass
	else:
		print ('Flag from %s to %s : ints from %d to %d' % (t1['vex'],t2['vex'],flag_startint,flag_stopint))
		zapints_list.append('%d:%d' % (flag_startint,flag_stopint))

def date_to_jd(year,month,day):
	""" Convert a date to Julian Day.
	Algorithm from 'Practical Astronomy with your Calculator or Spreadsheet', 
	4th ed., Duffet-Smith and Zwart, 2011. """
	# Copied from from https://gist.githubusercontent.com/jiffyclub/1294443/raw/02f4d41b8fa6ba7b39c577490f1e41d841cddd0b/jdutil.py
	if month == 1 or month == 2:
		yearp = year - 1
		monthp = month + 12
	else:
		yearp = year
		monthp = month
	if ((year < 1582) or (year == 1582 and month < 10) or (year == 1582 and month == 10 and day < 15)):
		B = 0 # before start of Gregorian calendar
	else:
		A = math.trunc(yearp / 100.)
		B = 2 - A + math.trunc(A / 4.) #  after start of Gregorian calendar
	if yearp < 0:
		C = math.trunc((365.25 * yearp) - 0.75)
	else:
		C = math.trunc(365.25 * yearp)
	D = math.trunc(30.6001 * (monthp + 1))
	jd = B + C + D + day + 1720994.5
	return jd

#############################################################################################

(options, args) = parser.parse_args()
if len(args)!=3:
	parser.error('Please specify a vex file, source name, and telescope name!')

vex       = vex.Vex(args[0])
source    = args[1]
telescope = args[2].upper()
fb_Tint   = options.tint
fb_Tstart = options.startmjd
zapints_list = []

# Check properties of filterbank file
if options.filterbankfile != None:
	p1 = subprocess.Popen([HEADER_PROG, options.filterbankfile, '-tstart'],   stdout=subprocess.PIPE)
	p2 = subprocess.Popen([HEADER_PROG, options.filterbankfile, '-tsamp'],    stdout=subprocess.PIPE)
	p3 = subprocess.Popen([HEADER_PROG, options.filterbankfile, '-nsamples'], stdout=subprocess.PIPE)
	s1,err = p1.communicate()
	s2,err = p2.communicate()
	s3,err = p3.communicate()
	if err != None:
		print ('Error calling SIGPROC utility "header" (%s) : %s!' % (HEADER_PROG,str(err)))
		os.exit(1)
	fb_Tstart = float(s1)
	fb_Tint   = float(s2)*1e-6
	fb_Nsamp  = int(s3)
	print ('Determined start time (MJD %.12f) and sampling rate (%f usec) of %s' % (fb_Tstart,fb_Tint*1e6,options.filterbankfile))

# Convert integration time from seconds into MJD fraction
fb_Tint = fb_Tint / 86400.0

# The rfifind '-zapints' works on a block level.
# The time duration of an "interval" depends on the block length and rfifind '-blocks' setting.
# For filterbank files the blocksize is hard-coded to 480 (sigproc_fb.c: s->spectra_per_subint = 480),
fb_Tint = fb_Tint * 480 * options.blocks

# Parse the VEX file
scans  = vex['SCHED']
exper  = [x for x in vex['EXPER']][0]
Nscans = len(scans)
if Nscans<=1:
	print ('Number of VEX scans is %d, nothing to do!' % (Nscans))
	sys.exit(0)

# Keep only scans of the selected station
Nscans    = len(scans)
scannames = [x for x in scans]
for ii in range(Nscans):
	stations = scans[scannames[ii]].getall('station')
	if not any([st[0].upper() == telescope for st in stations]):
		del scans[scannames[ii]]
print ('Found %d VEX scans of which %d include %s.' % (Nscans,len(scans),telescope))
if len(scans) <= 1:
	print ('Nothing to do!')
	sys.exit(0)
scannames = [x for x in scans]
Nscans = len(scans)

# Make sure we have a 'fb_Tstart'
if fb_Tstart == None:
	print ('Warning: filterbank file start time unspecified! Assuming it equals start of first VEX scan.')
	fb_Tstart = scans[scannames[0]]['start']

# First and last scan
scan_1 = scans[scannames[0]]
scan_N = scans[scannames[Nscans-1]]

# Throw away all scans that are off-source
for ii in range(Nscans):
	if scans[scannames[ii]]['source'] != source:
		del scans[scannames[ii]]
Nscans    = len(scans)
scannames = [x for x in scans]
if Nscans<1:
	print ('No scans found on source %s. Nothing to do.' % (source))
	sys.exit(0)

# Zap all off-source time ranges now
zap_tstart = scan_1['start']
for ii in range(Nscans):
	stations = scans[scannames[ii]].getall('station')
	for st in stations:
		if not(st[0].upper() == telescope):
			continue
		zap_tstop = scans[scannames[ii]]['start']
		doFlag(zap_tstart, zap_tstop, fb_Tstart, fb_Tint, zapints_list)	

		# Start next zap from end of current on-source scan
		dur_sec = int(st[2].split()[0])
		zap_tstart = dateVEXaddsec(zap_tstop, dur_sec)
if scan_N['start'] != zap_tstart:
	zap_tstop = vex['EXPER'][exper]['exper_nominal_stop'] 
	doFlag(zap_tstart, zap_tstop, fb_Tstart, fb_Tint, zapints_list)

# Write out a 'rfifind -zapints a:b' command file
f = open('zapints.cmd', 'w')
f.write('rfifind -blocks %d -o %s -zapints ' % (options.blocks,exper))
for ii in range(len(zapints_list)):
	f.write(zapints_list[ii])
	if ii<(len(zapints_list)-1):
		f.write(',')

if options.filterbankfile != None:
	f.write(' -filterbank %s' % (options.filterbankfile))

f.close()
print('Wrote new zapints.cmd file')
