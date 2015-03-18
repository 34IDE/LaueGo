#!/usr/bin/env python

"""
put this into a folder named "Old VersionStatus"  which should be inside "Igor Shared Procedures"
"""

import subprocess
import sys
import os
import time
import shutil
import xml.dom.minidom
import hashlib

isoFmt = '%Y-%m-%dT%H:%M:%S'
VersionXmlName = 'VersionStatus.xml'
archiveExtensions = set(['.ipf', '.txt', '.py' , '.xml', '.app', '.exe', '', '.html', '.ifn', 'ihf', '.png', '.mtl', '.hlp', '.rtf'])	# the '' is only for executables

def MakeVersionStatusDict(FolderPath):
	""" Create a dict containing version status of all ipf, txt, py, xml, etc. files """
	global isoFmt, VersionXmlName, archiveExtensions
	if not os.path.isdir(FolderPath): return ''

	actualFiles = dict()
	for root, dirs, files in os.walk(FolderPath):
		if ignoreName(root): continue						# skip directories with old versions
		for dir in dirs:
			if ignoreName(dir):
				dirs.remove(dir)  							# don't visit old directories

		for dir in dirs:
			if os.path.splitext(dir.lower())[1]=='.app':
				dirs.remove(dir)  		# don't visit files inside of apps
				fullpath = os.path.join(FolderPath,dir)
				mtime = os.stat(fullpath).st_mtime				# modification time (second)
				timeStr = time.strftime(isoFmt,time.localtime(mtime))
				actualFiles[dir] = (float('nan'),timeStr,'','.app')

		relPath = os.path.relpath(root,FolderPath)
		for file in files:
			if len(file)<1: continue
			ext = os.path.splitext(file.lower())[1]
			if os.path.basename(file) == VersionXmlName: continue
			if not(ext in archiveExtensions): continue		# only look at certain files extensions
			fullpath = os.path.join(root,file)				# full path to this file
			if ext=='' and not (is_exe(fullpath)): continue	# empty extensions are only for executables
			if file[0]=='.': continue						# no files starting with '.', no hidden files
			if relPath=='.': rel = file
			else: rel = os.path.join(relPath,file)

			hexHash = ''
			try:
				fullpath = os.path.join(FolderPath,rel)
				statinfo = os.stat(fullpath)				# for file size in bytes & modified time
				flen = statinfo.st_size						# file size in bytes
				mtime = statinfo.st_mtime					# modification time (second)
				f = open(fullpath,'r')
				buf = f.read(flen)							# read the whole file, to get a good hash
				f.close()
				m = hashlib.md5()							# m = hashlib.sha256()
				m.update(buf)
				hexHash = m.hexdigest()
				timeStr = time.strftime(isoFmt,time.localtime(mtime))
			except:
				buf = ''

			version = float('nan')
			if ext=='.ipf':									# ONLY find versions of ipf files
				i = buf.find('#pragma version')+15			# at start of character after '#pragma version'
				if i>=0 :
					j = buf[i:i+30].find('=')				# find the '='
					if j>=0:
						i += j+1
						buf = buf[i:i+100].strip()			# version number must be within 100 of '='
						try:	version = float(buf.split(None,1)[0])
						except:	version = float('nan')

			actualFiles[rel] = (version,timeStr,hexHash,ext)

#			#			if len(file)>0 and version==version:			# valid file name and version number
#			if len(file)>0:									# valid file name and version number
#				if relPath=='.': rel = file
#				else: rel = os.path.join(relPath,file)
	return actualFiles


def WriteVersionStatus(FolderPath,actualFiles,now):
	"""write a list of all files with their versions to an xml file"""
	global isoFmt

	count = len(actualFiles)
	if count<1: return 1

	strDate = time.strftime('%A, %B %d, %Y',now)
	strTime = time.strftime('%I:%M:%S %p %Z',now)
	strIso = time.strftime(isoFmt,now)
	out = '<?xml version="1.0" encoding="UTF-8" ?>\n\n'
	out += '<VersionStatus xmlns="http://sector34.xor.aps.anl.gov/34ide:VersionStatus">\n'
	out += '\t<written date="'+strDate+'" time="'+strTime+'" isoTime="'+strIso+'"></written>\n'
	out += '\t<sourceFolder>'+FolderPath+'</sourceFolder>\n'
	out += '\t<fileCount>'+str(count)+'</fileCount>\n'
	for fName in actualFiles:
		try:
			version,timeStr,hexHash,ext = actualFiles[fName]
			out += '\t<file version="%g" time="%s" md5="%s" ext="%s">%s</file>\n' % (version,timeStr,hexHash,ext,fName)	# append line to out
		except:	pass
	out += "</VersionStatus>\n"

	fileName = os.path.join(FolderPath,VersionXmlName)		# name of file to to write

	# get the date-time from the existing xml file, if I can.  This is use to archive existing file
	try: 	xin = xml.dom.minidom.parse(fileName,parser=None)
	except:	xin = None
	isoTime = 'X'
	if not (xin is None):
		topNode = xin.childNodes[0]
		if topNode.nodeName.find('VersionStatus')!=0:
			print '"'+VersionXmlName+'" is wrong type of file!, Overwriting it.'
		else:
			for node in topNode.childNodes:
				if node.nodeName == 'written':
					isoTime = node.getAttribute('isoTime')
					isoTime = isoTime.replace(':','_')
					isoTime = isoTime.split('.')[0]
					break
			if len(isoTime)<1: isoTime = 'X'

	# location and name of file to archive
	archiveDir = os.path.join(FolderPath,'Old VersionStatus')	# directory to archive the current VersionXmlName
	archiveName = os.path.join(archiveDir,'VersionStatus_'+isoTime+'.xml')	# new name for archived file

	if os.path.isfile(fileName) and os.path.isdir(archiveDir):
		try:
			shutil.copyfile(fileName,archiveName)				# copy existing file to the archive with new name
			print 'archived old Version Status.\n'+str(count)+' files'
		except:
			print 'Failed to backup old '+VersionXmlName+',  copying:'
			print "from: ",fileName
			print "to: ",archiveName
			return 1
	elif not os.path.isfile(fileName):
		print 'No existing '+VersionXmlName+', so not archiving it.'
	elif not os.path.isdir(archiveDir):
		print 'NOT archiving the existing '+VersionXmlName+' since the archive folder does not exist, just overwriting it.'

	try:
		f = open(fileName,'w')							# write new VersionXmlName file
		f.write(out)
		f.close()
	except:
		print "Error writing new Version Status"
		return 1

	return 0


def ignoreName(name):								# ignore folders or files with this name
	low = name.lower()
	ignore = False
	if low.find('old')==0: ignore=True				# starts with old
	if '/old' in low: ignore=True
	if 'old/' in low: ignore=True
	if (len(low)-low.rfind('old'))==3: ignore=True	# ends with old
	if 'old version' in low: ignore=True			# skip definitely old versions
	if '(no copy)' in low: ignore=True				# obviously not this
	if 'obsolete' in low: ignore=True				# ignore obsolete procedures
	if low[0]=='.': ignore=True						# no hidden files or hidden directories

	if low=='new': ignore=True
	return ignore


def is_exe(fullpath):
	return os.path.isfile(fullpath) and os.access(fullpath,os.X_OK)


def getFolderPath():
	""" determine the FolderPath """
	filePath = os.path.realpath(__file__)
	FolderPath = os.path.dirname(filePath)
	if FolderPath.find('/Old VersionStatus')>0: FolderPath = FolderPath + '/..'
	FolderPath = os.path.abspath(FolderPath)			# path to "Igor Shared Procedures"
	return FolderPath


now = time.localtime()

try:	FolderPath = getFolderPath()
except:	FolderPath = ''

try:	actualFiles = MakeVersionStatusDict(FolderPath)
except:	actualFiles = dict()
if len(actualFiles)<1:
	print 'ERROR  --  Could not read info from any files'
	sys.exit(1)
	#	print 'read data on %d files from '+VersionXmlName+',  found data on %d files from disk' % (len(Vfiles), len(actualFiles))

try:		err = WriteVersionStatus(FolderPath,actualFiles,now)
except:	err = 1
if err: print 'ERROR'
sys.exit(err)
