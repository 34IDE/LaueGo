#!/usr/bin/env python

import sys
import os
import datetime
import socket
import numpy
import h5py

import epics
from epics import caget

VERSION = 1.06
NEXUS_VERSION = '4.2.1'
IMAGE_GET_TIMEOUT = 40.0			# timeout for getting the image

def getImageDebug(imagePVroot):
	image = numpy.zeros((100, 100), dtype='i2')
	image = image.astype('u2')
#	print 'type of image is', image.dtype
	return image

def getAllPVsDebug(PVroot):
	# for a 100 x 100 detector
	pvValues = {'BeamBad': 1, 'UndGap': 15.290624999999999, 'deviceName': 'Undulator_#34s_3.3cm', 'scanNum': 17247, \
	'wirebaseZ': -707.1234873607259, 'Nx': 100, 'Ny': 100, 'WslitHeight': 0.19987500000000002, 'incident_energy': 17.36679695974342, \
	'UndTaper': -0.0006960578688381958, 'RingCurrent': 102.23893516369534, 'undulatorTaper': -0.0006960578688381958, 'title': '', \
	'MonoMode': 'white slitted', 'I0': 1238, 'ScalerClockFreq': 10000000.0, 'sampleDistance': 0.0, 'endy': 99, 'I0_calc': 1238.0, 'endx': 99, \
	'I_start_calc': 57.0, 'WslitWidth': 0.09993750000000001, 'X1': 509.90000000000003, 'Z1': -3014.5, 'topUp': 1, 'wirebaseY': 707.0855460951844, \
	'ScalerCountTime': 1.0, 'CCDshutter': 1, 'HutchTemperature': 23.144444444444442, 'sampleName': '', 'I_final': 1.0, 'sizeX': 100, 'sizeY': 100, \
	'Y1': -1044.4, 'foil': 0, 'I_final_calc': 1.0, 'gain': 2.0, 'wirebaseX': -51999.950000000004, 'wirescan': 0.0, 'detectorModel': 'XRD0820', \
	'L32_height': 5.0, 'exposure': 3.0, 'userName': 'Xu', 'startx': 0, 'starty': 0, 'ScalerClock_calc': 10000000.0, 'I_start': 57.0, \
	'wireX': -51999.950000000004, 'wireY': 707.0855460951844, 'wireZ': -707.1234873607259, 'top_up': 1, 'LightOn': 0.0, \
	'detectorID': 'PE1621 723-3335', 'biny': 1, 'binx': 1, 'L32_width': 0.20150000000001, 'undulatorGap': 15.290624999999999, \
	'badge': 12345, 'email': 'bob@trump.edu', 'proposal': 'G 12345', 'phone': '630-252-5555', 'affiliation': 'Trump University'}
	return pvValues



def getImage(imagePVroot):
	global IMAGE_GET_TIMEOUT

	try:	port = caget(imagePVroot+'PortName_RBV')
	except:	raise NameError('Could not find \"%r\"' % (imagePVroot+'PortName_RBV',))
	if len(port)<1: raise NameError('Could not find \"%r\"' % (imagePVroot+'PortName_RBV',))

	try:
		NDdimsPV = epics.PV(imagePVroot+'NDimensions_RBV')
		numberOfDims = int(NDdimsPV.get())
	except:
		raise ValueError('Could not get dimensions')

	if not (numberOfDims == 2):
		raise ValueError('number of dimensions for image must be 2')

	try:
		dimsPV = epics.PV(imagePVroot+'Dimensions_RBV')
		dims = dimsPV.get()
		Nx = int(dims[0])
		Ny = int(dims[1])
	except:
		raise ValueError('cannot get the dimensions of the array')

	if Nx<=0 or Ny<=0 or Nx>4096 or Ny>4096:
		raise ValueError('Nx = %r and Ny = %r is not allowed' % (Nx,Ny))
	try:
		imagePV = epics.PV(imagePVroot+'ArrayData')
		image = imagePV.get(timeout=IMAGE_GET_TIMEOUT)
#		print 'type of image is', image.dtype
	except:
		raise ValueError('could not get the image to %s' % imagePVroot+'ArrayData')

	try:
#		image = image.astype('u2')
#		print 'type of image is', image.dtype
		image.resize((Ny,Nx), refcheck=False)
	except:
		raise ValueError('could not resize the image to %d x %d' % (Nx,Ny))

	return image


def getAllPVs(PVroot):
	# epics.caget(pvname[, as_string=False[, count=None[, as_numpy=True[, timeout=None[, use_monitor=False ] ] ] ] ])

	try:	port = caget(PVroot+':cam1:PortName_RBV')
	except:	raise NameError('Could not find \"%r\"' % (PVroot+':cam1:PortName_RBV',))
	if len(port)<1: raise NameError('Could not find \"%r\"' % (PVroot+':cam1:PortName_RBV',))

	pvValues = {}
	pvValues['Nx'] = caget(PVroot+':cam1:MaxSizeX_RBV')
	pvValues['Ny'] = caget(PVroot+':cam1:MaxSizeY_RBV')

	pvValues['startx'] = caget(PVroot+':ROI1:MinX')
	pvValues['sizeX'] = caget(PVroot+':ROI1:SizeX')
	pvValues['binx'] = caget(PVroot+':ROI1:BinX')

	pvValues['starty'] = caget(PVroot+':ROI1:MinY')
	pvValues['biny'] = caget(PVroot+':ROI1:BinY')
	pvValues['sizeY'] = caget(PVroot+':ROI1:SizeY')
	pvValues['detectorModel'] = caget(PVroot+':cam1:Model_RBV')
	pvValues['exposure'] = caget(PVroot+':cam1:AcquireTime_RBV')
	try:	pvValues['gain'] = float(caget(PVroot+':cam1:PEGain'))
	except:	pass

	pvValues['title'] = caget('34ide:title')
	pvValues['userName'] = caget('34ide:userName')
	pvValues['sampleName'] = caget('34ide:sampleName')
	pvValues['scanNum'] = caget('34ide:scanNum')

	# these 5 added January 2017
	pvValues['badge'] = caget('34ide:userBadge')				# facility_user_id
	pvValues['email'] = caget('34ide:userEmail')				# email
	pvValues['proposal'] = caget('34ide:proposal')				# proposal
	pvValues['phone'] = caget('34ide:userPhone')				# telephone_number
	pvValues['affiliation'] = caget('34ide:userAffiliation')	# affiliation

	pvValues['RingCurrent'] = caget('S:SRcurrentAI')
	pvValues['undulatorGap'] = caget('ID34:Gap.VAL')
	pvValues['undulatorTaper'] = caget('ID34:TaperGap.VAL')
	pvValues['topUp'] = caget('Mt:TopUpAutoEnableC')
	pvValues['MonoMode'] = caget('34ide:monoT:transStatus', as_string=True)
	pvValues['HutchTemperature'] = caget('34ide:Geomechanics.J')
	pvValues['BeamBad'] = caget('34ide:beamBad')
	pvValues['CCDshutter'] = caget('34ide:9440:c0:Out0')
	pvValues['LightOn'] = caget('34ide:IPC1:ch8')
	pvValues['foil'] = caget('34ide:foil')
	pvValues['incident_energy'] = caget('34ide:monoE:keV.RBV')
	pvValues['sampleDistance'] = caget('34ide:Keyence:1:ch1')

	pvValues['X1'] = caget('34ide:mxv:c0:m2.RBV')
	pvValues['Y1'] = caget('34ide:mxv:c0:m3.RBV')
	pvValues['Z1'] = caget('34ide:mxv:c0:m1.RBV')

	pvValues['I0'] = caget('34ide:vsc:c0.S3')
	pvValues['I_start'] = caget('34ide:vsc:c0.S4')
	pvValues['I_final'] = caget('34ide:vsc:c0.S2')
	pvValues['I0_calc'] = caget('34ide:vsc:c0_cts1.C')
	pvValues['I_start_calc'] = caget('34ide:vsc:c0_cts1.D')
	pvValues['I_final_calc'] = caget('34ide:vsc:c0_cts1.B')
	pvValues['ScalerCountTime'] = caget('34ide:vsc:c0.T')
	pvValues['ScalerClock_calc'] = caget('34ide:vsc:c0_cts1.A')
	pvValues['ScalerClockFreq'] = caget('34ide:vsc:c0.FREQ')

	pvValues['UndGap'] = caget('ID34:Gap.VAL')
	pvValues['UndTaper'] = caget('ID34:TaperGap.VAL')
	pvValues['top_up'] = caget('Mt:TopUpAutoEnableC')
	pvValues['deviceName'] = caget('ID34:Device')

	pvValues['L32_height'] = caget('34ida:m58:c1:m5.RBV')
	pvValues['L32_width'] = caget('34ida:m58:c1:m4.RBV')
	pvValues['WslitWidth'] = caget('34ide:m58:c0:m5.RBV')
	pvValues['WslitHeight'] = caget('34ide:m58:c0:m7.RBV')

	pvValues['wireX'] = caget('34ide:wireStages.G')
	pvValues['wireY'] = caget('34ide:wireStages.H')
	pvValues['wireZ'] = caget('34ide:wireStages.I')
	pvValues['wirescan'] = caget('34ide:aero:c0:m1.RBV')
	pvValues['wirebaseX'] = caget('34ide:t80:c0:m1.RBV')
	pvValues['wirebaseY'] = caget('34ide:t80:c0:m3.RBV')
	pvValues['wirebaseZ'] = caget('34ide:t80:c0:m2.RBV')

	try:
		endx = pvValues['startx'] + pvValues['sizeX'] - 1
		endy = pvValues['starty'] + pvValues['sizeY'] - 1
		pvValues['endx'] = endx
		pvValues['endy'] = endy
	except:
		pass

	#	key_PVs += "endx=;endy=;xdim=;ydim=;"
	pvValues['detectorID'] = 'PE1621 723-3335'
	return pvValues


################################################################################
################################################################################
#		Start of stuff that actually does some work
#

def writeHDF5(image,destFile,pvValues):
	global VERSION, NEXUS_VERSION

	destFile = os.path.realpath(destFile)
	now = datetime.datetime.now()
	now = now.replace(microsecond=0)

	h5Version = str(h5py.version._h5.get_libversion())
	h5Version = h5Version.replace('(','')	# change (1, 8, 13) to 1.8.13
	h5Version = h5Version.replace(')','')
	h5Version = h5Version.replace(',','.')
	h5Version = h5Version.replace(' ','')

	f = h5py.File(destFile,'w')
	f.attrs.create('HDF5_Version',h5Version)
	f.attrs.create('NeXus_Version',NEXUS_VERSION)
	f.attrs.create('creator',__file__ + ' ' + str(VERSION))
	f.attrs.create('file_name',destFile)
	f.attrs.create('file_time',now.isoformat())

	facility = f.create_group('Facility')
	facility.attrs.create('NX_class','Facility')
	facility.create_dataset('facility_name', data='APS')
	facility.create_dataset('facility_sector', data='XSD/SSM')
	facility.create_dataset('facility_beamline', data='34ID-E')
	facility.create_dataset('facility_station', data='E')
	facility.create_dataset('facility_float', data=3.1415, dtype='f4')
	facility.create_dataset('facility_int', data=12, dtype='i4')

	entry = f.create_group('entry1')
	entry.attrs.create('NX_class','NXentry')
	entry.create_dataset('title', data=pvValues['title'])
	entry.create_dataset('scanNum', data=pvValues['scanNum'], dtype='i4')

	user = entry.create_group('user')
	user.attrs.create('NX_class','NXuser')
	user.create_dataset('name', data=pvValues['userName'])
	user.create_dataset('facility_user_id', data=pvValues['badge'])		# these 5 added January 2017
	user.create_dataset('email', data=pvValues['email'])
	user.create_dataset('proposal', data=pvValues['proposal'])
	user.create_dataset('telephone_number', data=pvValues['phone'])
	user.create_dataset('affiliation', data=pvValues['affiliation'])

	dataGrp = entry.create_group('data')
	dataGrp.attrs.create('NX_class','NXdata')
	dest = dataGrp.create_dataset('data', data=image)
	dest.attrs.create('signal', data=1)

	detector = entry.create_group('detector')
	detector.attrs.create('NX_class','NXdetector')
	detector.create_dataset('Nx', data=pvValues['Nx'], dtype='i4')
	detector.create_dataset('Ny', data=pvValues['Ny'], dtype='i4')
	detector.create_dataset('startx', data=pvValues['startx'], dtype='i4')
	detector.create_dataset('starty', data=pvValues['starty'], dtype='i4')
	detector.create_dataset('endx', data=pvValues['endx'], dtype='i4')
	detector.create_dataset('endy', data=pvValues['endy'], dtype='i4')
	detector.create_dataset('binx', data=pvValues['binx'], dtype='i4')
	detector.create_dataset('biny', data=pvValues['biny'], dtype='i4')
	detector.create_dataset('exposure', data=pvValues['exposure'])
	try:	dest = detector.create_dataset('gain', data=pvValues['gain'])
	except:	pass	;		# not all detectors have a gain
	dest.attrs.create('units', data='pF')
	detector.create_dataset('make', data='Perkin Elmer')
	detector.create_dataset('model', data=pvValues['detectorModel'])
	detector.create_dataset('ID', data=pvValues['detectorID'])

	monitor = entry.create_group('monitor')
	monitor.attrs.create('NX_class','NXmonitor')
	monitor.create_dataset('mode', data='timer')
	monitor.create_dataset('I0', data=pvValues['I0'], dtype='i4')
	monitor.create_dataset('I_start', data=pvValues['I_start'], dtype='i4')
	monitor.create_dataset('I_final', data=pvValues['I_final'], dtype='i4')
	monitor.create_dataset('I0_calc', data=pvValues['I0_calc'])
	monitor.create_dataset('I_start_calc', data=pvValues['I_start_calc'])
	monitor.create_dataset('I_final_calc', data=pvValues['I_final_calc'])
	monitor.create_dataset('ScalerCountTime', data=pvValues['ScalerCountTime'])
	monitor.create_dataset('ScalerClock_calc', data=pvValues['ScalerClock_calc'])
	monitor.create_dataset('ScalerClockFreq', data=pvValues['ScalerClockFreq'])

	instrument = entry.create_group('microDiffraction')
	instrument.attrs.create('NX_class','NXinstrument')
	source = instrument.create_group('source')
	source.attrs.create('NX_class','NXsource')
	dest = source.create_dataset('distance', data=65.0)
	dest.attrs.create('units', data='m')
	dest = source.create_dataset('current', data=pvValues['RingCurrent'])
	dest.attrs.create('units', data='mA')
	dest = source.create_dataset('gap', data=pvValues['UndGap'])
	dest.attrs.create('units', data='mm')
	dest = source.create_dataset('taper', data=pvValues['UndTaper'])
	dest.attrs.create('units', data='mm')
	source.create_dataset('probe', data='x-ray')
	source.create_dataset('type', data='Synchrotron X-ray Source')
	source.create_dataset('top_up', data=pvValues['top_up'], dtype='i2')
	source.create_dataset('name', data=pvValues['deviceName'])
	dest = source.create_dataset('L32_width', data=pvValues['L32_width'])
	dest.attrs.create('units', data='mm')
	dest = source.create_dataset('L32_height', data=pvValues['L32_height'])
	dest.attrs.create('units', data='mm')

	dest = instrument.create_dataset('Wslit_width', data=pvValues['WslitWidth'])
	dest.attrs.create('units', data='mm')
	dest = instrument.create_dataset('Wslit_height', data=pvValues['WslitHeight'])
	dest.attrs.create('units', data='mm')
	filter = instrument.create_group('filter')
	filter.attrs.create('NX_class', data='NXfilter')
	filter.create_dataset('description', data=pvValues['foil'], dtype='i2')
	instrument.create_dataset('CCDshutter', data=pvValues['CCDshutter'], dtype='i2')
	instrument.create_dataset('MonoMode', data=pvValues['MonoMode'])
	dest = instrument.create_dataset('HutchTemperature', data=pvValues['HutchTemperature'])
	dest.attrs.create('units', data='C')
	instrument.create_dataset('BeamBad', data=pvValues['BeamBad'], dtype='i4')
	instrument.create_dataset('LightOn', data=pvValues['LightOn'], dtype='i4')

	sample = entry.create_group('sample')
	sample.attrs.create('NX_class','NXsample')
	sample.create_dataset('name', data=pvValues['sampleName'])
	dest = sample.create_dataset('distance', data=0.0)
	dest.attrs.create('units', data='mm')
	dest = sample.create_dataset('sampleX', data=pvValues['X1'])
	dest.attrs.create('units', data='micron')
	dest = sample.create_dataset('sampleY', data=pvValues['Y1'])
	dest.attrs.create('units', data='micron')
	dest = sample.create_dataset('sampleZ', data=pvValues['Z1'])
	dest.attrs.create('units', data='micron')
	dest = sample.create_dataset('incident_energy', data=pvValues['incident_energy'])
	dest.attrs.create('units', data='keV')
	dest = sample.create_dataset('sampleDistance', data=pvValues['sampleDistance'])
	dest.attrs.create('units', data='mm')

	wire = entry.create_group('wire')
	wire.attrs.create('NX_class','NXsample')
	wire.create_dataset('wireX', data=pvValues['wireX'])
	wire.create_dataset('wireY', data=pvValues['wireY'])
	wire.create_dataset('wireZ', data=pvValues['wireZ'])
	wire.create_dataset('wirescan', data=pvValues['wirescan'])
	wire.create_dataset('wirebaseX', data=pvValues['wirebaseX'])
	wire.create_dataset('wirebaseY', data=pvValues['wirebaseY'])
	wire.create_dataset('wirebaseZ', data=pvValues['wirebaseZ'])

	f.close()

#
#		End of stuff that actually does the work
################################################################################
################################################################################


"""
try testing with:		./getADimage.py 34idePE2:image1 ~/Desktop
 clear ; ./getADimage.py 34idePE2:image1 ~/Developer/AreaDetectorImage/tempImage
"""


if __name__ == '__main__':
	try:
		imagePVroot = sys.argv[1]
		imagePVroot = imagePVroot.rstrip()
		imagePVroot = imagePVroot.rstrip(':')
		imagePVroot += ':'
	except:
		print 'No pv name given for the image, '
		print '  usage: getADimage.py 34idePE2:image1 [detination folder]'
		print '     or: getADimage.py dp_pilatus3:image1 [detination folder]'
		raise NameError('No pv name given')

	try:
		destFile = os.path.abspath(sys.argv[2])	# destination folder
		if len(os.path.splitext(destFile)[1]) == 0: destFile += '.h5'
	except:
		destFile = os.path.join(os.getcwd(),'tempImage.h5')

	if os.path.isdir(destFile):
		raise NameError('destFile = "%s" which is a folder' % destFile)

	destFldr,name = os.path.split(destFile)
	if not os.path.isdir(destFldr):			# check that destination folder exists
		raise NameError('the folder = "%s" does not exist' % destFldr)

	if not os.access(destFldr, os.W_OK):	# check that destination folder is writeable
		raise NameError('the folder = "%s" is not writeable' % destFldr)
	try:
		PVroot = imagePVroot.split(':')[0]
	except:
		PVroot = imagePVroot

	#	print 'imagePVroot = ',imagePVroot
	#	print 'PVroot = ',PVroot
	#	print 'destFile = ',destFile

	if socket.gethostname().find('.xray.aps.anl.gov') < 0:
		print 'Using debug values'
		image = getImageDebug(imagePVroot)
		pvValues = getAllPVsDebug(PVroot)
	else:
		try:
			pvValues =  getAllPVs(PVroot)
			image = getImage(imagePVroot)
		except:
			exit(1)

	writeHDF5(image,destFile,pvValues)



