#!/usr/bin/env python

"""
put this into a folder named "Old VersionStatus"  which should be inside "Igor Shared Procedures" or "Local Packages", ...
"""

import subprocess
import sys
import os
import time
import datetime
import xml.etree.ElementTree		# http://docs.python.org/2/library/xml.etree.elementtree.html
import hashlib

isoFmt = '%Y-%m-%dT%H:%M:%S'
VersionXmlName = 'VersionStatus.xml'
archiveExtensions = set(['.ipf', '.txt', '.py' , '.xml', '.app', '.exe', '', '.html', '.ifn', 'ihf', '.png', '.mtl', '.hlp', '.rtf'])	# the '' is only for executables

def compareFiles(Vfiles, actualFiles):		# returns a string detailing mis-matches, '' means a match
#	diffsTimesOK = diffsTimes = diffsVersion = diffsMissingOnDisk = diffsActual = ''
	diffsHash = diffsVersion = diffsMissingOnDisk = diffsActual = ''

	for vf in Vfiles:
		try:	vvers,visoTime,vHash,ext = Vfiles[vf]
		except:
			vvers = float('nan')
			vHash = visoTime = ext = ''
		try:	avers,aisoTime,aHash,ext = actualFiles[vf]
		except:
			avers = float('nan')
			aHash = aisoTime = ext = ''

		if len(vHash)>0 and len(aHash)==0:					# found file in VersionStatus, but not on disk 
			diffsMissingOnDisk += VersionXmlName+' has "%s" version=%g  at %s,  but file missing on disk\n' % (vf,vvers,visoTime)

		elif (vvers != avers) and (vvers==vvers or avers==avers) :	# versions differ
#		elif vvers != avers:								# versions differ
			if not (isNaN(vvers) and isNaN(avers)):
				diffsVersion += 'Mismatch on "%s" \tversion xml=%g, actual=%g,  \txmlTime=%s,  actualTime=%s\n' % (vf,vvers,avers,visoTime,aisoTime)

		elif not(vHash == aHash):							# hashes differ
			diffsHash += 'Hash Mismatch on "%s" \tversion xml=%g  \txmlTime=%s,  actualTime=%s\n' % (vf,vvers,visoTime,aisoTime)

	for af in actualFiles:									# next check that all files on disk are in VersionXmlName
		try:	avers,aisoTime,aHash,ext = actualFiles[af]
		except:
			avers = float('nan')
			aHash = aisoTime = ext = ''
		try:	vvers,visoTime,vHash,ext = Vfiles[af]
		except:
			vvers = float('nan')
			vHash = visoTime = ext = ''

		if len(vHash)==0 and len(aHash)>0:					# file on disk not found in VersionStatus
			diffsActual += 'On Disk, found file "%s" \tversion=%g  at %s,  but no entry in "%s"\n' % (af,avers,aisoTime,VersionXmlName)

	diffs = ''
	if len(diffsVersion)>0:
		diffs += '*** Versions Differ\n'+diffsVersion

	if len(diffsActual)>0:
		if len(diffs)>0: diffs += '\n'
		diffs += '*** Extra file on Disk\n'+diffsActual

	if len(diffsMissingOnDisk)>0:
		if len(diffs)>0: diffs += '\n'
		diffs += '*** Missing file on Disk\n'+diffsMissingOnDisk

	if len(diffsHash)>0:
		if len(diffs)>0: diffs += '\n'
		diffs += '*** Hashes Differ (but versions match!), so content differs\n'+diffsHash

	return diffs


def VersionStatusFile2dict(FolderPath):
	""" read the VersionXmlName file, and use its contents to fill a dict """
	fullFile = os.path.join(FolderPath,VersionXmlName)
	xmlns = '{http://sector34.xor.aps.anl.gov/34ide:VersionStatus}'	# the name space

	try:	tree = xml.etree.ElementTree.parse(fullFile)
	except:	tree = None
	root = tree.getroot()
	if not(root.tag == xmlns+'VersionStatus'): return None

	isoTimeWritten = None
	files = dict()
	fileCount = 0
	for child in root:
		if child.tag == xmlns+'written':
			try:	isoTimeWritten = child.attrib['isoTime']
			except:	isoTimeWritten = None

		elif child.tag == xmlns+'fileCount':
			try:	fileCount = int(child.text)
			except:	fileCount = 0

		elif child.tag == xmlns+'file':
			try:	vers = float(child.attrib['version'])
			except:	vers = None
			try:	isoTime = child.attrib['time']
			except:	isoTime = None
			try:	hexHash = child.attrib['md5']
			except:	hexHash = ''
			try:	fileName = child.text
			except:	fileName = ''
			try:	ext = child.attrib['ext']
			except:	ext = ''
			files[fileName] = (vers,isoTime,hexHash,ext)

	if not(fileCount == len(files)):
		print 'ERROR  --  reading '+VersionXmlName+', fileCount = %d, but I read %d files' % (fileCount,len(files))
		files = dict()

	return files


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

# def MakeVersionStatusDict(IgorSharedPath):
# 	""" Create a dict containing version status of all ipf files """
# 	global  isoFmt
# 
# 	if not os.path.isdir(IgorSharedPath): return ''
# 
# 	cmd = 'cd "'+IgorSharedPath+'" ; grep -lR  "#pragma version" *'
# 	#	print "cmd = ---"+cmd+"---"
# 	proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, )
# 	stdout_value, stderr_value = proc.communicate('through stdin to stdout')
# 	if len(stderr_value)>0:
# 		print >> sys.stderr, 'ERROR --',repr(stderr_value)
# 		return ''										# in case of error, assume 1 core
# 	gresult = stdout_value
# 
# 	files = dict()
# 	for ipfName in gresult.splitlines():
# 		lowName = ipfName.lower()
# 		if not('.ipf' in lowName): continue					# only ipf files
# 		if '/old' in lowName: continue						# skip '/old' directories
# 		if 'old/' in lowName: continue						# skip 'old/' directories
# 		if 'old versions' in lowName: continue					# skip definitely old versions
# 		try:
# 			fpath = os.path.join(IgorSharedPath,ipfName)
# 			statinfo = os.stat(fpath)						# for file size in bytes & modified time
# 			flen = statinfo.st_size							# file size in bytes
# 			mtime = statinfo.st_mtime						# modification time (second)
# 			f = open(fpath,'r')
# 			buf = f.read(flen)								# read the whole file, to get a good hash
# 			f.close()
# #			m = hashlib.sha256()
# 			m = hashlib.md5()
# 			m.update(buf)
# 			hexHash = m.hexdigest()
# 			timeStr = time.strftime(isoFmt,time.localtime(mtime))
# 			# flen = min(4096,flen)
# 			# buf = buf[0:flen]								# trim down buf to seach for #pragma version
# 		except:
# 			buf = ''
# 
# 		i = buf.find('#pragma version')+15					# at start of character after '#pragma version'
# 		if i<0: continue
# 		j = buf[i:i+30].find('=')							# find the '='
# 		if j<0: continue
# 		i += j+1
# 		buf = buf[i:i+100].strip()							# version number must be within 100 of '='
# 		try:	version = float(buf.split(None,1)[0])
# 		except:	version = float('nan')
# 		if len(ipfName)>0 and version==version:				# valid file name and version number
# 			files[ipfName] = (version,timeStr,hexHash)
# 
# 	return files



def ignoreName(name):								# ignore folders or files with this name
	low = name.lower()
	ignore = False
	if low.find('old')==0: ignore=True				# starts with old
	if '/old' in low: ignore=True
	if 'old/' in low: ignore=True
	if (len(low)-low.rfind('old'))==3: ignore=True	# ends with old
	if 'old version' in low: ignore=True			# skip definitely old versions
	if '(no copy)' in low: ignore=True				# obviously not this
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
	FolderPath = os.path.abspath(FolderPath)				# path to "Igor Shared Procedures"
	return FolderPath


def isNaN(num):
	return not(num==num)


try:
	FolderPath = sys.argv[1]				# FolderPath passed as argument
except:
	try:	FolderPath = getFolderPath()	# no argument, use usual search
	except:	FolderPath = ''					# default is where you are, this is very unlikely

try:	Vfiles = VersionStatusFile2dict(FolderPath)
except:	Vfiles = dict()
try:	err = len(Vfiles)<1
except:	err = True
if err:
	print 'ERROR  --  Could not read file '+VersionXmlName
	sys.exit(err)

if not err:
	try:	actualFiles = MakeVersionStatusDict(FolderPath)
	except:	actualFiles = dict()
	err = err or len(actualFiles)<1
if err:
	print 'ERROR  --  Could read file '+VersionXmlName
	sys.exit(err)
#	print 'read data on %d ipf files from '+VersionXmlName+',  found data on %d ipf files from disk' % (len(Vfiles), len(actualFiles))

diffs = compareFiles(Vfiles,actualFiles)
if len(diffs)>0:
	print '\n********************* Mismatch *********************\n'
	print diffs
	print '********************* Mismatch *********************'
else:
	print 'Perfect Match,\n%d files' % len(actualFiles)

if err: print 'ERROR'
sys.exit(err)
