#!/usr/bin/env python
# -*- coding: utf-8 -*-

#version 1.06

from xdrlib import *
import sys
import os
import string

# fmt0 = '%18.7f'		# this was the old version
fmt0 = '%.12g'

# test with:		clear ; ./mdaAscii2IgorTest.py 2idd_0013.mda outfile.itx

def mdaAscii3d_IGOR(d,outFile=None):
	"""
	mdaAscii3d_IGOR(d) - for input specified MDA scan data structure d, this function
	will extract all 3D arrays and sequentially saved in IGOR ASCII format,
	it returns the output IGOR file name which contains the list of all 3D arrays saved 

	e.g.
	The 2idd_0013.mda contains 3D array:

		from mdaAscii import *
		d = readMDA("2idd_0013.mda",maxdim=3)
		fn = mdaAscii2d_IGOR(d)
	"""
	global fmt0
	fmt = fmt0 + '\t'
	if d[0]['rank']!=3: return ''

#	Variable Nx=d[3].npts, Ny=d[2].npts, Nz=d[1].npts,  Make/N=(Nx,Ny,Nz) data
	Nx = d[3].npts
	Ny = d[2].npts
	Nz = d[1].npts

	zPVs = ''
	for i in range(d[1].np):	# d[i].p is of class scanPositioner
		if i>0: zPVs += ','
		zPVs += str(d[1].p[i].name)
	yPVs = ''
	for i in range(d[2].np):
		if i>0: yPVs += ','
		yPVs += str(d[2].p[i].name)
	xPVs = ''
	for i in range(d[3].np):
		if i>0: xPVs += ','
		xPVs += str(d[3].p[i].name)

	try:									# this try-except is needed in case there are no positioners in d[1]
		if d[1].p[0].step_mode=='LINEAR':	# layers are linearly scaled
			z0 = float(d[1].p[0].data[0])
			dz = float(d[1].p[0].data[1]) - z0
			NzDone = d[1].curr_pt
			if NzDone > 1:
				dz = float(d[1].p[0].data[NzDone-1]) - z0
				dz /= (NzDone-1)

		zUnit = d[1].p[0].unit
		zName = d[1].p[0].desc
	except:
		NzDone = d[1].curr_pt				# if no positioner, then just 0,1,2,...
		z0 = float(0.0)
		dz = float(1.0)
		zUnit = ''
		zName = 'dummy'
		isSaturn = False
		isMCA = False

	try:
		if len(zName)<1: zName = d[1].p[0].name
	except:
		zName = 'Z-axis'

	if d[2].p[0].step_mode=='LINEAR':		# layers are linearly scaled
		y0 = float(d[2].p[0].data[0][0])
		dy = float(d[2].p[0].data[0][1]) - y0
		if NzDone>1:						# the first set of y positions is complete
			dy = float(d[2].p[0].data[0][Ny-1]) - y0
			dy /= (Ny-1)
	yUnit = d[2].p[0].unit
	yName = d[2].p[0].desc
	try:
		if len(yName)<1: yName = d[2].p[0].name
	except:
		yName = 'Y-axis'

	xName = 'X-axis'
	if d[3].np > 0:								# there are real positioners
		if d[3].p[0].step_mode=='LINEAR':		# layers are linearly scaled
			x0 = float(d[3].p[0].data[0][0][0])
			dx = float(d[3].p[0].data[0][0][1]) - x0
			if NzDone>1:						# the first set of x positions is complete
				dx = float(d[3].p[0].data[0][0][Nx-1]) - x0
				dx /= (Nx-1)
		xUnit = d[3].p[0].unit
		xName = d[3].p[0].desc
		try:
			if len(xName)<1: xName = d[3].p[0].name
		except:
			xName = 'X-axis'

	elif d[3].np == 0 and d[3].scan_name.find('scanH') >= 0:	# axis is computed, an MCA spectrum
		try:
			x0 = d[0]['dxpSaturn2idd:mca1.CALO'][2][0]
			dx = d[0]['dxpSaturn2idd:mca1.CALS'][2][0]
			xName = 'Saturn'
			xUnit = 'eV'
			isSaturn = True
			isMCA = False
		except:
			try:
				x0 = d[0]['2idd:mca1.CALO'][2][0]
				dx = d[0]['2idd:mca1.CALS'][2][0]
				xName = 'MCA'
				xUnit = 'eV'
				isSaturn = False
				isMCA = True
			except:
				x0 = 0.0
				dx = 1.0
				xName = 'X-axis'
				xUnit = ''
				isSaturn = False
				isMCA = False
		if dx<0.1:				# change from keV --> eV
			dx *= 1000.0
			x0 *= 1000.0
	else:
		return ''

	if xUnit == 'micron': xUnit = 'µm'
	if yUnit == 'micron': yUnit = 'µm'
	if zUnit == 'micron': zUnit = 'µm'

	if 0:
		print ' '
		print "Z axis, from %g, step=%g,  zName='%s', units='%s'" % (z0,dz,zName,zUnit)
		print "Y axis, from %g, step=%g,  yName='%s', units='%s'" % (y0,dy,yName,yUnit)
		print "X axis, from %g, step=%g,  xName='%s', units='%s'" % (x0,dx,xName,xUnit)
	if 0:
		# d[i] is of class scanDim
		print 'd[1] = ', d[1]
		print 'd[2] = ', d[2]
		print 'd[3] = ', d[3]
		print ' '
		print ' npts[] = ', d[1].npts, d[2].npts, d[3].npts
		print ' curr_pt[] = ', d[1].curr_pt, d[2].curr_pt, d[3].curr_pt
		print ' scan_name[] = ', d[1].scan_name, d[2].scan_name, d[3].scan_name
		print ' time[] = ', d[1].time, d[2].time, d[3].time
		print ' np[] = ', d[1].np, d[2].np, d[3].np
		print ' nt[] = ', d[1].nt, d[2].nt, d[3].nt
		print ' nd[] = ', d[1].nd, d[2].nd, d[3].nd
		print ' d[].dim =',d[1].dim,'   ',d[2].dim,'   ',d[3].dim
		print ' d[].rank =',d[1].rank,'   ',d[2].rank,'   ',d[3].rank
		# d[i].p is of class scanClass, a list of scanPositioner

		print '\nfor d[1]'
		print 'd[1].np =',d[1].np,'   number of positioners'
		for i in range(d[1].np):	# d[i].p is of class scanPositioner
			print 'd[1].p[i] =',str(d[1].p[i])

		print '\nfor d[2]'
		print 'd[2].np =',d[2].np,'   number of positioners'
		for i in range(d[2].np):	# d[i].p is of class scanPositioner
			print 'd[2].p[i] =',str(d[2].p[i])

		print '\nfor d[3]'
		print 'd[3].np =',d[3].np,'   number of positioners'
		for i in range(d[3].np):	# d[i].p is of class scanPositioner
			print 'd[3].p[i] =',str(d[3].p[i])

	if d[0]['rank'] < 3 : return ""
	#	print d[0].keys()

	path,fname = os.path.split(d[0]['filename'])
	froot = string.split(fname,'.mda')
	if outFile is None:
		if os.path.exists('ASCII') == 0 :
			os.mkdir('ASCII')
		dir = os.getcwd()
		outFile = dir+os.sep+'ASCII'+ os.sep + froot[0] + '_IGOR.txt'
	else:
		pp,ff = os.path.split(outFile)
		if len(pp)<1: outFile = os.path.join(path,outFile)


	# make sure there's room for the names, etc.
#	nz = d[1].curr_pt
#	ny = d[2].npts
#	nx = d[3].npts
#	# for i in range(np): print '%d \t%s \t\t%s  (%s)\n' % (i, d[2].p[i].name,d[2].p[i].desc,d[2].p[i].unit)
#	py = d[1].p[0].data
#	px = d[2].p[0].data
#	px = px[0]
#	xunit = d[2].p[0].unit
#	yunit = d[1].p[0].unit
#	xdesc = d[2].p[0].desc
#	ydesc = d[1].p[0].desc
#	xPV = d[2].p[0].name
#	yPV = d[1].p[0].name

	data = d[3].d[0].data
	if 0:
		print type(data),len(data)
		print type(data[0]),len(data[0])
		print type(data[0][0]),len(data[0][0])
		print data[0][0][0]

	prefix = ''				# probably something like '2idd:'
	for item in d[0].keys():
		i = item.find('saveData_fileSystem') 
		if i>=0:
			prefix = item[0:i]
			break

	baseNoteStr = ''						# things common to all wave notes
	try:
		s = d[0]['filename']
		if len(s): baseNoteStr += ('filename=%s;' % s)
	except: pass
	try:
		s = d[0][prefix+'saveData_comment1'][2]
		if len(s)>0: baseNoteStr += ('saveData_comment1=%s;' % s)
	except: pass
	try:
		s = d[0][prefix+'saveData_comment2'][2]
		if len(s)>0: baseNoteStr += ('saveData_comment2=%s;' % s)
	except: pass
	try:
		s = d[0][prefix+'saveData_fileSystem'][2]
		if len(s)>0: baseNoteStr += ('filesystem=%s;' % s)
	except: pass
	try:
		s = d[0][prefix+'scaler1.TP']
		tp = s[2][0]
		if type(tp) is float:
			baseNoteStr += ('TP=%g;' % tp)
	except: pass
	try:
		s = d[0][prefix+'userStringCalc10.AA'][2]
		if len(s)>0: baseNoteStr += ('title=%s;' % s)
	except: pass
	try:
		s = d[0][prefix+'userStringCalc10.BB'][2]
		if len(s)>0: baseNoteStr += ('user=%s;' % s)
	except: pass
	try:
		s = d[0][prefix+'userStringCalc10.CC'][2]
		if len(s)>0: baseNoteStr += ('sample=%s;' % s)
	except: pass

	fo = open(outFile,"w")
	fo.write('IGOR'+'\n')

	# write the full 3D data (probably a lot of MCA spectra)
	nd = d[3].nd			# number of detectors, probably only one MCA, this is not the ROIs
	for m in range(nd):		# loop over each detector (probably one for MCA data)
		wName = '%s_%s_mda' % (d[3].d[m].fieldName, froot[0])
		fo.write('Waves/O/D/N=(%d,%d,%d) %s\n' % (Nx,Ny,NzDone,wName))
		fo.write('BEGIN\n')
		data = d[3].d[m].data
		for iz in range(NzDone):
			for ix in range(Nx):
				for iy in range(Ny):
#					fo.write('%18.7f ' % data[iz][iy][ix])
					fo.write(fmt % data[iz][iy][ix])
				fo.write('\n')			# at end of each row
			fo.write('\n')				# at end of each layer
		fo.write('END\n')

		fo.write('X SetScale/P x %s, %s, "%s", %s\n' % (x0,dx,xUnit,wName))
		fo.write('X SetScale/P y %s, %s, "%s", %s\n' % (y0,dy,yUnit,wName))
		fo.write('X SetScale/P z %s, %s, "%s", %s\n' % (z0,dz,zUnit,wName))
		if len(d[3].d[m].unit)>0:
			fo.write(('X SetScale/I d 0,0,"%s", %s\n' % (d[3].d[m].unit,wName)),)

		waveClass = 'mda'
		if isSaturn or isMCA: waveClass = 'mca'
		noteStr = 'waveClass='+waveClass+';'+baseNoteStr		# things only for the data
		noteStr += ('pv=%s;' % d[3].d[m].name)
		if len(zName)>0: noteStr += ('zAxisName=%s;' % zName)
		if len(yName)>0: noteStr += ('yAxisName=%s;' % yName)
		if len(xName)>0: noteStr += ('xAxisName=%s;' % xName)
		if len(zPVs): noteStr += ('zPVs=%s;' % zPVs)
		if len(yPVs): noteStr += ('yPVs=%s;' % yPVs)
		if len(xPVs): noteStr += ('xPVs=%s;' % xPVs)

#		if len(d[2].d[m].desc)>0: noteStr += ('desc=%s;' % d[2].d[m].desc)
		fo.write(('X Note/K %s, "%s"\n' % (wName,noteStr)),)
		fo.write('\n')

	# done writing the main data, next write the ROI waves
	if isSaturn: mca = 'dxpSaturn'+prefix+'mca1'
	elif isMCA: mca = prefix+'mca1'
	else: 		return outFile 			# no mca, all done
	
	ROIs=[]
	ROIlabels=[]
	for m in range(25):		# loop over each detector (probably one for MCA data)
		label = mca+'.R'+str(m)+'NM'
		try:
			if len(d[0][label][2])>1:
				ROIlabels.append(d[0][label][2])
				ROIs.append(mca+'.R'+str(m))
		except: pass
	if 0:
		print mca
		print ROIs
		print ROIlabels,'   ',len(ROIlabels)


	#	ROIs[] contain the pvs,  want d[2].d[ROIindex].name == ROIs
	#	for each ROIs[i], need to search through d[2].d[j].name to find the j that matches.  then append j to ROIindex
	ROIindex = []
	for i in range(d[2].nd):
		try:
			j = ROIs.index(d[2].d[i].name)
			if j>=0: ROIindex.append(i)
		except: pass
	nd = len(ROIindex)

	wName = '%s_%s_ROIs' % (d[3].d[0].fieldName, froot[0])
	fo.write('Waves/O/D/N=(%d,%d,%d) %s\n' % (Ny,NzDone,nd,wName))	# one layer for each ROI
	fo.write('BEGIN\n')

	if 0:
		print len(d[2].d)				# 65
		print len(d[2].d[0].data)		# 10
		print len(d[2].d[0].data[0])	# 26
		print type(d[2].d[0].data[0][0])
		print 'len(ROIindex) = ',len(ROIindex)
		print 'Nx = ',Nx,'      Ny = ',Ny

	#	d[2].d[nd].data[NzDone][Ny]
	for ilayer in ROIindex:
		for iy in range(Ny):
			for iz in range(NzDone):
#				fo.write('%18.7f ' % d[2].d[ilayer].data[iz][iy])
				fo.write(fmt % d[2].d[ilayer].data[iz][iy])
			fo.write('\n')			# at end of each row
		fo.write('\n')				# at end of each layer
	fo.write('END\n')

	fo.write('X SetScale/P x %s, %s, "%s", %s\n' % (x0,dx,yUnit,wName))
	fo.write('X SetScale/P y %s, %s, "%s", %s\n' % (y0,dy,zUnit,wName))

	noteStr = 'waveClass=ROI;'+baseNoteStr		# things only for the data
	if len(zName)>0: noteStr += ('yAxisName=%s;' % zName)
	if len(yName)>0: noteStr += ('xAxisName=%s;' % yName)
	if len(yPVs): noteStr += ('zPVs=%s;' % zPVs)
	if len(xPVs): noteStr += ('yPVs=%s;' % yPVs)

	fo.write(('X Note/K %s, "%s"\n' % (wName,noteStr)),)
	for i in range(nd):
		fo.write('X SetDimLabel 2, %d, $"%s", %s\n' % (i,ROIlabels[i],wName))

	fo.close()
	return outFile 

def mdaAscii2d_IGOR(d,outFile=None):
	"""
	mdaAscii2d_IGOR(d) - for input specified MDA scan data structure d, this function
	will extract all 2D image arrays and sequentially saved in IGOR ASCII format,
	it returns the output IGOR file name which contains the list of all 2D images saved 

	e.g.
	The 2idd_0004.mda contains 2D array:

		from mdaAscii import *
		d = readMDA("2idd_0004.mda")
		fn = mdaAscii2d_IGOR(d)
	"""
	global fmt0
	fmt = fmt0 + '\t'
	if d[0]['rank'] < 2 : return ""
	#	print d[0].keys()

	path,fname = os.path.split(d[0]['filename'])
	froot = string.split(fname,'.mda')
	if outFile is None:
		if os.path.exists('ASCII') == 0 :
			os.mkdir('ASCII')
		dir = os.getcwd()
		outFile = dir+os.sep+'ASCII'+ os.sep + froot[0] + '_IGOR.txt'
		# print outFile
	else:
		pp,ff = os.path.split(outFile)
		if len(pp)<1: outFile = os.path.join(path,outFile)

	# number of positioners, detectors
	np = d[2].np
	nd = d[2].nd
	min_column_width = 16
	# make sure there's room for the names, etc.
	ny = d[1].curr_pt
	nx = d[2].npts
	# for i in range(np): print '%d \t%s \t\t%s  (%s)\n' % (i, d[2].p[i].name,d[2].p[i].desc,d[2].p[i].unit)

	try:	py = d[1].p[0].data
	except:	py = range(ny)
	try:	yunit = d[1].p[0].unit
	except: yunit = ''
	try:	ydesc = d[1].p[0].desc
	except:	ydesc = ''
	try:	yPV = d[1].p[0].name
	except:	yPV = ''

	px = d[2].p[0].data
	px = px[0]
	xunit = d[2].p[0].unit
	xdesc = d[2].p[0].desc
	xPV = d[2].p[0].name


	if xunit == 'micron': xunit = 'µm'
	if yunit == 'micron': yunit = 'µm'


	fo = open(outFile,"w")
	fo.write('IGOR'+'\n')
	for i in range(nd):
		fo.write(('Waves/O/D/N='),)
		fo.write(('(%s, ' % nx),)
		fo.write(('%s) ' % ny),)
		wName = '%s_%s_mda' % (d[2].d[i].fieldName, froot[0])
		fo.write(('%s' % wName),)
		fo.write('\n')

		fo.write('BEGIN\n')
		data = d[2].d[i].data
		for j in range(nx):
#			for k in range(ny): fo.write(('%18.7f ' % data[k][j]),)
			for k in range(ny): fo.write((fmt % data[k][j]),)
			fo.write('\n')
		fo.write('END\n')

		fo.write(('X SetScale/I x %s, %s, "%s", %s\n' % (px[0],px[nx-1],xunit,wName)),)
		fo.write(('X SetScale/I y %s, %s, "%s", %s\n' % (py[0],py[ny-1],yunit,wName)),)
		if len(d[2].d[i].unit)>0:
			fo.write(('X SetScale/I d 0,0,"%s", %s\n' % (d[2].d[i].unit,wName)),)

		noteStr = ('pv=%s;' % d[2].d[i].name)
		if len(d[2].d[i].desc)>0: noteStr += ('desc=%s;' % d[2].d[i].desc)
		if len(xdesc)>0: noteStr += ('xLabel=%s;' % xdesc)
		if len(ydesc)>0: noteStr += ('yLabel=%s;' % ydesc)
		if len(xPV)>0: noteStr += ('xPV=%s;' % xPV)
		if len(yPV)>0: noteStr += ('yPV=%s;' % yPV)
		fo.write(('X Note/K %s, "%s"\n' % (wName,noteStr)),)
		fo.write('\n')
	fo.close()
	return outFile 

def mdaAscii1d_IGOR(d,outFile=None):
	"""
	mdaAscii1d_IGOR(d) - for input specified MDA scan data structure d, this function
	will extract all 1D arrays and sequentially saved in IGOR ASCII format,
	it returns the output IGOR file name which contains the list of all 1D arrays saved 

	e.g.
	The 2idd_0004.mda contains 2D array:

		from mdaAscii import *
		d = readMDA("2idd_0004.mda")
		fn = mdaAscii1d_IGOR(d)
	"""
	global fmt0
	fmt = fmt0 + '\n'
	if d[0]['rank'] < 1 : return ""
	#	print d[0].keys()

	path,fname = os.path.split(d[0]['filename'])
	froot = string.split(fname,'.mda')
	if outFile is None:
		if os.path.exists('ASCII') == 0 :
			os.mkdir('ASCII')
		dir = os.getcwd()
		outFile = dir+os.sep+'ASCII'+ os.sep + froot[0] + '_IGOR.txt'
		# print outFile
	else:
		pp,ff = os.path.split(outFile)
		if len(pp)<1: outFile = os.path.join(path,outFile)

	# number of positioners, detectors
	np = d[1].np
	nd = d[1].nd
	min_column_width = 16
	# make sure there's room for the names, etc.
	nx = d[1].npts
	# for i in range(np): print '%d \t%s \t\t%s  (%s)\n' % (i, d[2].p[i].name,d[2].p[i].desc,d[2].p[i].unit)
	px = d[1].p[0].data
	xunit = d[1].p[0].unit
	xdesc = d[1].p[0].desc
	xPV = d[1].p[0].name

	if xunit == 'micron': xunit = 'µm'

	fo = open(outFile,"w")
	fo.write('IGOR'+'\n')

	# write the positioner wave first
	fo.write(('Waves/O/D/N='),)
	fo.write(('(%s) ' % nx),)
	wName = '%s_%s_mda' % (d[1].p[0].fieldName, froot[0])
	fo.write(('%s' % wName),)
	fo.write('\n')

	fo.write('BEGIN\n')
	data = d[1].p[0].data
#	for j in range(nx): fo.write(('%18.7f\n' % data[j]),)
	for j in range(nx): fo.write((fmt % data[j]),)
	fo.write('END\n')

	fo.write(('X SetScale/I x %s, %s, "%s", %s\n' % (px[0],px[nx-1],xunit,wName)),)
	if len(xunit)>0:
		fo.write(('X SetScale/I d 0,0,"%s", %s\n' % (xunit,wName)),)

	try:	timeStr = mdaTime2ISOtime(d[1].time)
	except:	timeStr = ''
	noteStr = ('pv=%s;' % d[1].p[0].name)
	if len(xdesc)>0: noteStr += ('desc=%s;' % xunit)
	if len(xPV)>0: noteStr += ('xPV=%s;' % xPV)
	if len(timeStr) > 0: noteStr += ('file_time=%s;' % timeStr)
	fo.write(('X Note/K %s, "%s"\n' % (wName,noteStr)),)
	fo.write('\n')

	# write each of the dectors
	for i in range(nd):
		fo.write(('Waves/O/D/N='),)
		fo.write(('(%s) ' % nx),)
		wName = '%s_%s_mda' % (d[1].d[i].fieldName, froot[0])
		fo.write(('%s' % wName),)
		fo.write('\n')

		fo.write('BEGIN\n')
		data = d[1].d[i].data
#		for j in range(nx): fo.write(('%18.7f\n' % data[j]),)
		for j in range(nx): fo.write((fmt % data[j]),)
		fo.write('END\n')

		fo.write(('X SetScale/I x %s, %s, "%s", %s\n' % (px[0],px[nx-1],xunit,wName)),)
		if len(d[1].d[i].unit)>0:
			fo.write(('X SetScale/I d 0,0,"%s", %s\n' % (d[1].d[i].unit,wName)),)

		noteStr = ('pv=%s;' % d[1].d[i].name)
		if len(d[1].d[i].desc)>0: noteStr += ('desc=%s;' % d[1].d[i].desc)
		if len(xdesc)>0: noteStr += ('xLabel=%s;' % xdesc)
		if len(xPV)>0: noteStr += ('xPV=%s;' % xPV)
		if len(timeStr) > 0: noteStr += ('file_time=%s;' % timeStr)
		fo.write(('X Note/K %s, "%s"\n' % (wName,noteStr)),)
		fo.write('\n')
	fo.close()
	return outFile 


def mdaTime2ISOtime(mdaTime):
	""" convert the time format passed in by mda to an ISO8601 format
	also, round to nearest second """
	if not (type(mdaTime) is str): return ''
	elif len(mdaTime)<1: return ''

	mdaTime = mdaTime.upper()
	mdaTime = mdaTime.replace(',',' ')
	mdaTime = mdaTime.replace('  ',' ')
	mda = mdaTime.split(' ')

	months = ['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
	mon = months.index(mda[0]) + 1
	day = int(mda[1])
	year = int(mda[2])

	hms = mda[3].split(':')
	hr = int(hms[0])
	mn = int(hms[1])
	try:
		sec = int(round(float(hms[2])))
	except:	
		sec = int(0)

	isoStr =  '%04d-%-2d-%02dT%02d:%02d:%02d' % (year,mon,day,hr,mn,sec)
	return isoStr



#!/usr/bin/env python
#
# version 4 Tim Mooney 8/22/02
#
# from xdrlib import *
# import sys
# import os
# import string

class scanDim:
	def __init__(self):
		self.rank = 0
		self.dim = 0
		self.npts = 0
		self.curr_pt = 0
		self.scan_name = ""
		self.time = ""
		self.np = 0				# number of positioners in p[]
		self.p = []				#   list of scanPositioner instances
		self.nd = 0				# number of detectors in d[]
		self.d = []				#   list of scanDetector instances
		self.nt = 0				# number of triggers in t[]
		self.t = []				#   list of scanTrigger instances

	def __str__(self):
		if self.scan_name <> '':
			s = "%dD data from scan_name = \"%s\": %d/%d pts; %d positioners, %d detectors, %d triggers" % (
				self.dim, self.scan_name, self.curr_pt, self.npts, self.np, self.nd,
				self.nt)
		else:
			s = "%dD data (not read in)" % (self.dim)

		return s

class scanClass:
	def __init__(self):
		self.rank = 0
		self.npts = 0
		self.curr_pt = 0
		self.plower_scans = 0
		self.name = ""
		self.time = ""
		self.np = 0
		self.nd = 0
		self.nt = 0
		self.p = []
		self.d = []
		self.t = []

class scanPositioner:
	def __init__(self):
		self.number = 0
		self.fieldName = ""
		self.name = ""
		self.desc = ""
		self.step_mode = ""
		self.unit = ""
		self.readback_name = ""
		self.readback_desc = ""
		self.readback_unit = ""
		self.data = []

	def __str__(self):
		s = "positioner %d (%s), desc:%s, unit:%s\n" % (self.number, self.name,
			self.desc, self.unit)
		s = s + "   step mode: %s, readback:\"%s\"\n" % (self.step_mode,
			self.readback_name)
		s = s + "data:%s" % (str(self.data))
		return s

class scanDetector:
	def __init__(self):
		self.number = 0
		self.fieldName = ""
		self.name = ""
		self.desc = ""
		self.unit = ""
		self.data = []

	def __str__(self):
		s = "detector %d (%s), desc:%s, unit:%s, data:%s\n" % (self.number,
			self.name, self.desc, self.unit, str(self.data))
		return s

class scanTrigger:
	def __init__(self):
		self.number = 0
		self.name = ""
		self.command = 0.0

	def __str__(self):
		s = "trigger %d (%s), command=%f\n" % (self.number,
			self.name, self.command)
		return s


def detName(i,new=0):
	"""
	detName(i,new=0) - this function returns the detector name Di used in sscan record
	where 
	  i - specify the zero based detector sequence number 
	  new - 1 specify the version 5 Di names desired, default 0
	"""
	if new:
		return "D%02d"%(i+1)
		return
	if i < 15:
		return string.upper("D%s"%(hex(i+1)[2]))
	elif i < 85:
		return "D%02d"%(i-14)
	else:
		return "?"

def posName(i):
	"""
	posName(i) - this function returns the positioner name Pi used in sscan record
	where 
	  i - specify the zero based positioner sequence number 
	"""
	if i < 4:
		return "P%d" % (i+1)
	else:
		return "?"

def readScan(file, v, new=0):
	"""
	readScan(file, v, new=0) - internal scan read routine, it unpack a subset of scan data from
	the current position of the file pointer 
	it returns the scan data set extracted
	where
	  file - file pointer of an opened MDA file
	  v - input verbose specified
	  new  - default 0, if 1 specified then version 5 Di name used
	"""
	scan = scanClass()
	buf = file.read(10000) # enough to read scan header
	u = Unpacker(buf)
	scan.rank = u.unpack_int()
	if v: print "scan.rank = ", `scan.rank`
	scan.npts = u.unpack_int()
	if v: print "scan.npts = ", `scan.npts`
	scan.curr_pt = u.unpack_int()
	if v: print "scan.curr_pt = ", `scan.curr_pt`
	if (scan.rank > 1):
		# if curr_pt < npts, plower_scans will have garbage for pointers to
		# scans that were planned for but not written
		scan.plower_scans = u.unpack_farray(scan.npts, u.unpack_int)
		if v: print "scan.plower_scans = ", `scan.plower_scans`
	namelength = u.unpack_int()
	scan.name = u.unpack_string()
	if v: print "scan.name = ", `scan.name`
	timelength = u.unpack_int()
	scan.time = u.unpack_string()
	if v: print "scan.time = ", `scan.time`
	scan.np = u.unpack_int()
	if v: print "scan.np = ", `scan.np`
	scan.nd = u.unpack_int()
	if v: print "scan.nd = ", `scan.nd`
	scan.nt = u.unpack_int()
	if v: print "scan.nt = ", `scan.nt`
	for j in range(scan.np):	# loop over the positioners
		scan.p.append(scanPositioner())
		scan.p[j].number = u.unpack_int()
		scan.p[j].fieldName = posName(scan.p[j].number)
		if v: print "positioner ", j
		length = u.unpack_int() # length of name string
		if length: scan.p[j].name = u.unpack_string()
		if v: print "scan.p[%d].name = %s" % (j, `scan.p[j].name`)
		length = u.unpack_int() # length of desc string
		if length: scan.p[j].desc = u.unpack_string()
		if v: print "scan.p[%d].desc = %s" % (j, `scan.p[j].desc`)
		length = u.unpack_int() # length of step_mode string
		if length: scan.p[j].step_mode = u.unpack_string()
		if v: print "scan.p[%d].step_mode = %s" % (j, `scan.p[j].step_mode`)
		length = u.unpack_int() # length of unit string
		if length: scan.p[j].unit = u.unpack_string()
		if v: print "scan.p[%d].unit = %s" % (j, `scan.p[j].unit`)
		length = u.unpack_int() # length of readback_name string
		if length: scan.p[j].readback_name = u.unpack_string()
		if v: print "scan.p[%d].readback_name = %s" % (j, `scan.p[j].readback_name`)
		length = u.unpack_int() # length of readback_desc string
		if length: scan.p[j].readback_desc = u.unpack_string()
		if v: print "scan.p[%d].readback_desc = %s" % (j, `scan.p[j].readback_desc`)
		length = u.unpack_int() # length of readback_unit string
		if length: scan.p[j].readback_unit = u.unpack_string()
		if v: print "scan.p[%d].readback_unit = %s" % (j, `scan.p[j].readback_unit`)

	for j in range(scan.nd):	# loop over the detectors
		scan.d.append(scanDetector())
		scan.d[j].number = u.unpack_int()
		scan.d[j].fieldName = detName(scan.d[j].number,new=new)
		if v: print "detector ", j
		length = u.unpack_int() # length of name string
		if length: scan.d[j].name = u.unpack_string()
		if v: print "scan.d[%d].name = %s" % (j, `scan.d[j].name`)
		length = u.unpack_int() # length of desc string
		if length: scan.d[j].desc = u.unpack_string()
		if v: print "scan.d[%d].desc = %s" % (j, `scan.d[j].desc`)
		length = u.unpack_int() # length of unit string
		if length: scan.d[j].unit = u.unpack_string()
		if v: print "scan.d[%d].unit = %s" % (j, `scan.d[j].unit`)

	for j in range(scan.nt):	# loop over the triggers
		scan.t.append(scanTrigger())
		scan.t[j].number = u.unpack_int()
		if v: print "trigger ", j
		length = u.unpack_int() # length of name string
		if length: scan.t[j].name = u.unpack_string()
		if v: print "scan.t[%d].name = %s" % (j, `scan.t[j].name`)
		scan.t[j].command = u.unpack_float()
		if v: print "scan.t[%d].command = %s" % (j, `scan.t[j].command`)

	### read data
	# positioners
	file.seek(file.tell() - (len(buf) - u.get_position()))
	buf = file.read(scan.np * scan.npts * 8)
	u = Unpacker(buf)
	for j in range(scan.np):
		if v: print "read %d pts for pos. %d at buf loc %x" % (scan.npts,
			j, u.get_position())
		scan.p[j].data = u.unpack_farray(scan.npts, u.unpack_double)    
		if v: print "scan.p[%d].data = %s" % (j, `scan.p[j].data`)
        
	# detectors
	file.seek(file.tell() - (len(buf) - u.get_position()))
	buf = file.read(scan.nd * scan.npts * 4)
	u = Unpacker(buf)
	for j in range(scan.nd):
		scan.d[j].data = u.unpack_farray(scan.npts, u.unpack_float)
		if v: print "scan.d[%d].data = %s" % (j, `scan.d[j].data`)

	return scan

def readMDA(fname=None, maxdim=2, verbose=1, help=0, new=0):
	"""
	readMDA(fname=None, maxdim=2, verbose=1, help=0) - This fuction reads an MDA
	file and constructs the MDA data structure accordingly 
	it returns the MDA data sturcture constructed where
	  fname - specifies the input mda file name
	  maxdim - specifies the max dimension extract, default 2 
	  verbose - reading info on or off, default 1
	  help - echo help information on or off, default 0 
	  new - 1 specify the version 5 Di names desired, default 0

	e.g.

	from readMDA import *
	d = readMDA('/home/beams/CHA/data/xxx/cha_0001.mda')

	"""
	dim = []

	if (fname == None):
		return dim
	if (not os.path.isfile(fname)): fname = fname + '.mda'
	if (not os.path.isfile(fname)):
		print "ERROR -- ",fname," is not a file"
		return dim

	file = open(fname, 'rb')
	#print "file = ", str(file)
	#file.seek(0,2)
	#filesize = file.tell()
	#file.seek(0)
	buf = file.read(100)		# to read header for scan of up to 5 dimensions
	u = Unpacker(buf)

	# read file header
	version = u.unpack_float()
	scan_number = u.unpack_int()
	rank = u.unpack_int()
	dimensions = u.unpack_farray(rank, u.unpack_int)
	isRegular = u.unpack_int()
	pExtra = u.unpack_int()
	pmain_scan = file.tell() - (len(buf) - u.get_position())

	for i in range(rank):
		dim.append(scanDim())
		dim[i].dim = i+1
		dim[i].rank = rank-i

	file.seek(pmain_scan)
	s0 = readScan(file, max(0,verbose-1),new=new)
	dim[0].npts = s0.npts
	dim[0].curr_pt = s0.curr_pt
	dim[0].scan_name = s0.name
	dim[0].time = s0.time
	dim[0].np = s0.np
	for i in range(s0.np): dim[0].p.append(s0.p[i])
	dim[0].nt = s0.nt
	for j in range(s0.nt): dim[0].t.append(s0.t[j])
	dim[0].nd = s0.nd
	for i in range(s0.nd): dim[0].d.append(s0.d[i])

	if ((rank == 1) and (maxdim > 1) and 0):
		# collect 1D data
		file.seek(s0.plower_scans)
		s = readScan(file, max(0,verbose-1),new=new)
		dim[0].npts = s.npts
		dim[0].curr_pt = s.curr_pt
		dim[0].scan_name = s.name
		dim[0].time = s.time
		# copy positioner, trigger, detector instances
		dim[0].np = s.np
		for j in range(s.np):
			dim[0].p.append(s.p[j])
			tmp = s.p[j].data[:]
			dim[0].p[j].data = []
			dim[0].p[j].data.append(tmp)
		dim[0].nt = s.nt
		for j in range(s.nt): dim[0].t.append(s.t[j])
		dim[0].nd = s.nd
		for j in range(s.nd):
			dim[0].d.append(s.d[j])
			tmp = s.d[j].data[:]
			dim[0].d[j].data = []
			dim[0].d[j].data.append(tmp)

	if ((rank > 1) and (maxdim > 1)):
		# collect 2D data
		for i in range(s0.curr_pt):
			file.seek(s0.plower_scans[i])
			s = readScan(file, max(0,verbose-1),new=new)
			if i == 0:
				dim[1].npts = s.npts
				dim[1].curr_pt = s.curr_pt
				dim[1].scan_name = s.name
				dim[1].time = s.time
				# copy positioner, trigger, detector instances
				dim[1].np = s.np
				for j in range(s.np):
					dim[1].p.append(s.p[j])
					tmp = s.p[j].data[:]
					dim[1].p[j].data = []
					dim[1].p[j].data.append(tmp)
				dim[1].nt = s.nt
				for j in range(s.nt): dim[1].t.append(s.t[j])
				dim[1].nd = s.nd
				for j in range(s.nd):
					dim[1].d.append(s.d[j])
					tmp = s.d[j].data[:]
					dim[1].d[j].data = []
					dim[1].d[j].data.append(tmp)
			else:
				# append data arrays
				for j in range(s.np): dim[1].p[j].data.append(s.p[j].data)
				for j in range(s.nd): dim[1].d[j].data.append(s.d[j].data)

	if ((rank > 2) and (maxdim > 2)):
		# collect 3D data
		for i in range(s0.curr_pt):
			file.seek(s0.plower_scans[i])
			s1 = readScan(file, max(0,verbose-1),new=new)
			for j in range(s1.curr_pt):
				file.seek(s1.plower_scans[j])
				s = readScan(file, max(0,verbose-1),new=new)
				if ((i == 0) and (j == 0)):
					dim[2].npts = s.npts
					dim[2].curr_pt = s.curr_pt
					dim[2].scan_name = s.name
					dim[2].time = s.time
					# copy positioner, trigger, detector instances
					dim[2].np = s.np
					for k in range(s.np):
						dim[2].p.append(s.p[k])
						tmp = s.p[k].data[:]
						dim[2].p[k].data = [[]]
						dim[2].p[k].data[i].append(tmp)
					dim[2].nt = s.nt
					for k in range(s.nt): dim[2].t.append(s.t[k])
					dim[2].nd = s.nd
					for k in range(s.nd):
						dim[2].d.append(s.d[k])
						tmp = s.d[k].data[:]
						dim[2].d[k].data = [[]]
						dim[2].d[k].data[i].append(tmp)
				elif j == 0:
					for k in range(s.np):
						dim[2].p[k].data.append([])
						dim[2].p[k].data[i].append(s.p[k].data)
					for k in range(s.nd):
						dim[2].d[k].data.append([])
						dim[2].d[k].data[i].append(s.d[k].data)
				else:
					# append data arrays
					for k in range(s.np): dim[2].p[k].data[i].append(s.p[k].data)
					for k in range(s.nd): dim[2].d[k].data[i].append(s.d[k].data)

	# Collect scan-environment variables into a dictionary
	dict = {}
	dict['sampleEntry'] = ("description", "unit string", "value")
	dict['filename'] = fname
	dict['rank'] = rank
	dict['dimensions'] = dimensions
	if pExtra:
		file.seek(pExtra)
		buf = file.read()       # Read all scan-environment data
		u = Unpacker(buf)
		numExtra = u.unpack_int()
		for i in range(numExtra):
			name = ''
			n = u.unpack_int()      # length of name string
			if n: name = u.unpack_string()
			desc = ''
			n = u.unpack_int()      # length of desc string
			if n: desc = u.unpack_string()
			unpackType = u.unpack_int()

			unit = ''
			value = ''
			count = 0
			if unpackType != 0:   # not DBR_STRING
				count = u.unpack_int()  # 
				n = u.unpack_int()      # length of unit string
				if n: unit = u.unpack_string()

			if unpackType == 0: # DBR_STRING
				n = u.unpack_int()      # length of value string
				if n: value = u.unpack_string()
			elif unpackType == 32: # DBR_CTRL_CHAR
				#value = u.unpack_fstring(count)
				v = u.unpack_farray(count, u.unpack_int)
				value = ""
				for i in range(len(v)):
					# treat the byte array as a null-terminated string
					if v[i] == 0: break
					value = value + chr(v[i])

			elif unpackType == 29: # DBR_CTRL_SHORT
				value = u.unpack_farray(count, u.unpack_int)
			elif unpackType == 33: # DBR_CTRL_LONG
				value = u.unpack_farray(count, u.unpack_int)
			elif unpackType == 30: # DBR_CTRL_FLOAT
				value = u.unpack_farray(count, u.unpack_float)
			elif unpackType == 34: # DBR_CTRL_DOUBLE
				value = u.unpack_farray(count, u.unpack_double)

			dict[name] = (desc, unit, value)

	dim.reverse()
	dim.append(dict)
	dim.reverse()

	if verbose:
		print "%s is a %d-D file; %d dimensions read in." % (fname, dim[0]['rank'], len(dim)-1)
		print "dim[0] = dictionary of %d scan-environment PVs" % (len(dim[0]))
		print "   usage: dim[0]['sampleEntry'] ->", dim[0]['sampleEntry']
		for i in range(1,len(dim)):
			print "dim[%d] = %s" % (i, str(dim[i]))
		print "   usage: dim[1].p[2].data -> 1D array of positioner 2 data"
		print "   usage: dim[2].d[7].data -> 2D array of detector 7 data"

	if help:
		print " "
		print "   each dimension (e.g., dim[1]) has the following fields: "
		print "      time      - date & time at which scan was started: %s" % (dim[1].time)
		print "      scan_name - name of scan record that acquired this dimension: '%s'" % (dim[1].scan_name)
		print "      curr_pt   - number of data points actually acquired: %d" % (dim[1].curr_pt)
		print "      npts      - number of data points requested: %d" % (dim[1].npts)
		print "      nd        - number of detectors for this scan dimension: %d" % (dim[1].nd)
		print "      d[]       - list of detector-data structures"
		print "      np        - number of positioners for this scan dimension: %d" % (dim[1].np)
		print "      p[]       - list of positioner-data structures"
		print "      nt        - number of detector triggers for this scan dimension: %d" % (dim[1].nt)
		print "      t[]       - list of trigger-info structures"

	if help:
		print " "
		print "   each detector-data structure (e.g., dim[1].d[0]) has the following fields: "
		print "      desc      - description of this detector"
		print "      data      - data list"
		print "      unit      - engineering units associated with this detector"
		print "      fieldName - scan-record field (e.g., 'D01')"


	if help:
		print " "
		print "   each positioner-data structure (e.g., dim[1].p[0]) has the following fields: "
		print "      desc          - description of this positioner"
		print "      data          - data list"
		print "      step_mode     - scan mode (e.g., Linear, Table, On-The-Fly)"
		print "      unit          - engineering units associated with this positioner"
		print "      fieldName     - scan-record field (e.g., 'P1')"
		print "      name          - name of EPICS PV (e.g., 'xxx:m1.VAL')"
		print "      readback_desc - description of this positioner"
		print "      readback_unit - engineering units associated with this positioner"
		print "      readback_name - name of EPICS PV (e.g., 'xxx:m1.VAL')"

	return dim


if __name__ == '__main__':
	if 0:
		print " "
		d = readMDA(sys.argv[1],help=1,verbose=1)
		print " "

	n = len(sys.argv)
	inFile = outFile = None
	if n>1: inFile = sys.argv[1]
	if n>2: outFile = sys.argv[2]

	d = readMDA(inFile,verbose=0,maxdim=3)
#	print d[0]
#	print "1111",d[1]
	try:	rank = d[0]['rank']
	except:	rank = None
	if rank == 1 : fn = mdaAscii1d_IGOR(d,outFile)
	elif rank == 2 : fn = mdaAscii2d_IGOR(d,outFile)
	elif rank == 3 : fn = mdaAscii3d_IGOR(d,outFile)
	else:
		print "ERROR -- "+__file__+" Only understand data of rank 1, 2, or 3, not d[0]['rank'] =",rank
		fn = ""

	print '%s' % fn







