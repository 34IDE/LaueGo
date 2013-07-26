#pragma rtGlobals=1		// Use modern global access method.
#pragma version = 0.04
#include "HDF5images"

// readin Flyscan HDF images, extract single images and export into single files
// requires HDF5images.ipf by Jon Tischler
//Function FlyImageExtract(fileFolder, filePrefix, fileNo_min, fileNo_max, sliceNo_min, sliceNo_max)
Function FlyImageExtract()
	String fileFolder, filePrefix
	String fileRange="1"
	Variable sliceNo_min=0, sliceNo_max=0
	String fileOpenFilter, filePath, fileName, outFolder, outfilePath, imageWaveName
	Variable refNum, fileNumPos, fileNum, fid
	
	//test
	//SetFileFolderInfo /D
	
	// get file path
	fileOpenFilter = "HDF5 Files (*.h5):.h5;All Files (*.*):.*"
	Open /D/R /F=fileOpenFilter /M = "Select one of the fly scan files" refNum
	if(strlen(S_fileName)>0)
		fileFolder = ParseFilePath(1,S_fileName, ":", 1,0)
		fileName = ParseFilePath(3,S_fileName,":",0,0)  // file name without extension
		fileNumPos = strsearch(fileName,"_",Inf,1) // search for "_" from end of name string
		filePrefix = fileName[0,fileNumPos-1]  // prefix of file name
		fileNum = str2num(fileName[fileNumPos+1,Inf])
	else
		return -1
	endif
	
	// get slice range
	STRUCT HDF5DataInfo di
	HDF5images#InitHDF5DataInfo(di)
	HDF5OpenFile/R fid as S_fileName
	HDF5DatasetInfo(fid,"/entry1/data/data",0,di)
	HDF5CloseFile fid
	Variable Nslices, ndims
	ndims = di.ndims
	if(ndims==3)
		Nslices = di.dims[0]
	else
		Nslices = 1
	endif
	
	// ask user for file number range and slice range
	fileRange=num2istr(fileNum)
	Prompt fileRange, "Range of flyscan file numbers (nonexistent files are skipped)"
	DoPrompt "Specify a range of input flyscan file number:", fileRange
	if (V_Flag) // V_Flag is set to 1 by DoPrompt if user canceled
		return -1
	endif
	
	Prompt sliceNo_min, "Min. slice number (first slice = 0):"
	Prompt sliceNo_max, "Max. slice number (last slice = " + num2str(Nslices-1) + "):"
	DoPrompt "Specify range of slices you want to extract from each fly scan", sliceNo_min, sliceNo_max
	if (V_Flag) // V_Flag is set to 1 by DoPrompt if user canceled
		return -1
	endif
	
	// ask user for output folder
	String outFolderOption
	Prompt outFolderOption, "Output options", popup "Use input folder;Select a different output folder"
	DoPrompt "Choose output folder", outFolderOption
	if (V_Flag) // V_Flag is set to 1 by DoPrompt if user canceled
		return -1
	endif
	if(cmpstr(outFolderOption,"Use input folder") != 0)
		do
			GetFileFolderInfo /D/Z=2	// /Z=2 surpresses popup of debug window upon user abort
			if (V_Flag == 0) // folder found correctly
				outFolder = S_path
			endif
		while ((V_Flag != 0) && (V_Flag != -1))				// repeat if folder not found
		if (V_Flag == -1)  // V_Flag is set to -1 by GetFileFolderInfo if user canceled
			return -1
		endif
	else
		outFolder = fileFolder
	endif
	

	// loop over fileRange
	String extras=""
	Variable i,j, nfiles = 0
	for (i=str2num(fileRange);numtype(i)==0;i=NextInRange(fileRange,i))
		filePath = fileFolder + filePrefix + "_" + num2str(i)+".h5"
		// check file existence
		Open /R /Z refNum as filePath
		if (V_flag == 0)
			Close refNum
		else
			continue
		endif
		//loop over slicerange
		for(j=sliceNo_min;j<=sliceNo_max;j+=1)
			// open file, read-in slice
			extras = ReplaceNumberByKey("slice",extras,j)
//			imageWaveName = LoadHDF5imageFile(filePath, slice = j)
//			imageWaveName = LoadGenericImageFile(filePath, extras=extras)
			imageWaveName = LoadHDF5imageFile(filePath, extras=extras)
			// wave = a slice
			// write wave to single file
			outfilePath = outFolder + filePrefix + "_" + num2str(i) + "_S_" + num2str(j)+".h5"
			Write1HDF5imageFile($imageWaveName,outfilePath)
			// kill this wave before memory explosion
			KillWaves $imageWaveName
			nfiles += 1
		endfor
	endfor
	
	print num2str(nfiles)+ " files extracted.\n"
end
