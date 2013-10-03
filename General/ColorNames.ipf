#pragma rtGlobals=1		// Use modern global access method.
#pragma version = 2.04
#pragma IgorVersion = 6.22
#pragma ModuleName=ColorNames
#pragma hide = 1



//	Function test_wavelength2RGB(lam)
//		Variable lam					// wavelength (nm)
//		Make/N=(10,3)/FREE rgbRange
//		rgbRange[0][0]= {0,9510,19275,0,0,65535,65535,65535,32767,0}
//		rgbRange[0][1]= {0,0,0,0,65535,65535,32767,0,0,0}
//		rgbRange[0][2]= {0,21845,33410,65535,0,0,0,0,0,0}
//		String str="", str2=""
//		Variable i,N=10
//		for (i=0;i<10;i+=1)
//			str += str2+"\""+RGBA2name(rgbRange[i][0],rgbRange[i][1],rgbRange[i][2],1,65535)+"\""
//			str2 = ",  "
//		endfor
//		print str
//	
//		Wave ww = wavelength2RGB(lam,65535)
//		printf "%g nm  -->  %s = \"%s\"\r",lam,vec2str(ww), RGBA2name(ww[0],ww[1],ww[2],1,65535)
//	End
//
Static Function/WAVE wavelength2RGB(lam,maxRGB)		// returns a true RGB for the wavelength (nm)
	Variable lam				// wavelength (nm)
	Variable maxRGB			// maximum of an RGB
	maxRGB = maxRGB>0 ? maxRGB : 65535						// default

	Make/N=3/W/U/FREE rgb = {0,0,0}
	if (numtype(lam))
		return rgb
	endif

	Make/N=10/FREE nmRange={0,380,415,472.5,532.5,580,605,685,750,825}	// centers of each color
	Make/N=(10,3)/FREE rgbRange
	rgbRange[0][0]= {0,9510,19275,0,0,65535,65535,65535,32767,0}
	rgbRange[0][1]= {0,0,0,0,65535,65535,32767,0,0,0}
	rgbRange[0][2]= {0,21845,33410,65535,0,0,0,0,0,0}

	Variable i0=BinarySearch(nmRange,lam)
	if (i0<0)					// requested wavelength is out of range, return black
		return rgb
	endif
	Variable mid=BinarySearchInterp(nmRange,lam)
	Variable i1=i0+1			// i0 < mid < i1
	Variable x0 = 1-mid+i0, 	x1 = 1-i1+mid

	Make/N=3/FREE rgb0, rgb1
	rgb0 = rgbRange[i0][p] * maxRGB/65535
	rgb1 = rgbRange[i1][p] * maxRGB/65535
	rgb = limit(round(x0*rgb0 + x1*rgb1),0,maxRGB)
	return rgb
	//	black	0-380 nm
	//	violet	380Ð450Ênm		415
	//	blue	450Ð495Ênm		472.5
	//	green	495Ð570Ênm		532.5
	//	yellow	570Ð590Ênm		580
	//	orange	590Ð620Ênm		605
	//	red		620Ð750Ênm		685
End
//
//	// Example:
//	Function/WAVE MakeRGBwave(wavelength)
//		Wave wavelength
//		Variable i, N=DimSize(wavelength,0)
//		Make/N=(N,3)/W/U/O spectrumRGB=NaN
//		for (i=0;i<N;i+=1)
//			Wave rgb = ColorNames#wavelength2RGB(wavelength[i],65535)
//			spectrumRGB[i][] = rgb[q]
//		endfor
//		return spectrumRGB
//	End



Function/WAVE colorCircleRGB(angle,[maxVal])// returns fully saturated colors
	Variable angle								// angle (rad)
	Variable maxVal
	maxVal = ParamIsDefault(maxVal) ? 1 : maxVal
	maxVal = maxVal>0 && numtype(maxVal)==0 ? maxVal : 1

	Make/N=3/FREE rgb
	Variable a									// a is always in range [0,1)
	angle = mod(angle,2*PI)					// angle must be in range [0,2*PI)
	angle += angle<0 ? 2*PI : 0
	if (numtype(angle))
		rgb = NaN
	elseif (angle< 2/3*PI)						// [0,60¡), red -> green
		a = angle / (2/3*PI)
		rgb = {1-a, a, 0}
	elseif (angle< 4/3*PI)						// [60,120¡)
		a = (angle-2/3*PI) / (2/3*PI)
		rgb = {0, 1-a, a}
	else											// [120,360¡)
		a = (angle-4/3*PI) / (2/3*PI)
		rgb = {a, 0, 1-a}
	endif
	rgb = sqrt(rgb)

	Variable m = WaveMax(rgb)
	rgb = limit(rgb*maxVal/m,0,maxVal)		// saturate the rgb, and clip to allowed range

	if (maxVal==65535)						// for 65535, use unsigned 16bit ints
		rgb = numtype(rgb) ? 0 : rgb
		rgb = round(rgb)
		Redimension/U/W rgb
	endif
	return rgb
End
//	Function test_colorCircleRGB(angle)
//		Variable angle							// (degree)
//		Wave rgb = colorCircleRGB(angle*PI/180)
//		printf "%g¡  ->  %s\r",angle,vec2str(rgb)
//	End
//
//	Function test_colorCircleRGB()
//		Variable N=361
//		Make/N=(N)/O xTest,yTest
//		Make/N=(N,3)/O/W/U rgbTest
//		SetScale/I x 0,2*PI,"", xTest,yTest
//	
//		xTest = cos(x)
//		yTest = sin(x)
//	
//		Variable i, angle
//		for (i=0;i<N;i+=1)
//			angle = i*DimDelta(xTest,0) + DimOffset(xTest,0)
//			Wave rgb = colorCircleRGB(angle,maxVal=65535)
//			rgbTest[i][] = rgb[q]
//		endfor
//	End


Function/WAVE symmetricWhiteRGB(val,[maxVal])	// red -> magenta -> white -> cyan -> blue
	Variable val									// value in range [-1,1]
	Variable maxVal
	maxVal = ParamIsDefault(maxVal) ? 1 : maxVal
	maxVal = maxVal>0 && numtype(maxVal)==0 ? maxVal : 1

	Make/N=3/D/FREE rgb=NaN, red={1,0,0}, magenta={1,0,1}, white={1,1,1}, cyan={0,1,1}, blue={0,0,1}
	val = 2*limit(val,-1,1)					// change range to [-2,2]
	if (val>1)									// go from cyan -> blue, (1,2]
		val += -1								// change range of val to [0,1]
		rgb = (1-val)*cyan + val*blue

	elseif (val>0)								// go from white -> cyan, (0,1]
		//	rgb = (1-val)*white + val*cyan
		rgb = sqrt(abs((1-val)*white + val*cyan))

	elseif (val>-1)								// go from magenta -> white, [-1,0)
		val += 1
		//	rgb = (1-val)*magenta + val*white
		rgb = sqrt(abs((1-val)*magenta + val*white))

	elseif (val>=-2)							// go from red -> magenta, [-2,1)
		val += 2
		rgb = (1-val)*red + val*magenta
	endif

	Variable m = WaveMax(rgb)
	rgb = limit(rgb*maxVal/m,0,maxVal)		// saturate the rgb, and clip to allowed range

	if (maxVal==65535)						// for 65535, use unsigned 16bit ints
		rgb = numtype(rgb) ? 0 : rgb
		rgb = round(rgb)
		Redimension/U/W rgb
	endif
	return rgb
End


Function/WAVE symmetricBlackRGB(val,[maxVal])	// red -> magenta -> black -> cyan -> blue
	Variable val									// value in range [-1,1]
	Variable maxVal
	maxVal = ParamIsDefault(maxVal) ? 1 : maxVal
	maxVal = maxVal>0 && numtype(maxVal)==0 ? maxVal : 1

	Make/N=3/D/FREE rgb=NaN, red={1,0,0}, magenta={1,0,1}, black={0,0,0}, cyan={0,1,1}, blue={0,0,1}
	val = 2*limit(val,-1,1)					// change range to [-2,2]
	Variable valIn = val
	if (val>1)									// go from cyan -> blue, (1,2]
		val += -1								// change range of val to [0,1]
		rgb = (1-val)*cyan + val*blue

	elseif (val>0)								// go from black -> cyan, (0,1]
		//	rgb = (1-val)*black + val*cyan
		rgb = sqrt(abs((1-val)*black + val*cyan))

	elseif (val>-1)								// go from magenta -> black, [-1,0)
		val += 1
		//	rgb = (1-val)*magenta + val*black
		rgb = sqrt(abs((1-val)*magenta + val*black))

	elseif (val>=-2)							// go from red -> magenta, [-2,1)
		val += 2
		rgb = (1-val)*red + val*magenta
	endif

	Variable m = abs(valIn)<1 ? 1 : WaveMax(rgb)	// saturate colors for abs(vals[i]) > 1
	rgb = limit(rgb*maxVal/m,0,maxVal)		// clip to allowed range

	if (maxVal==65535)						// for 65535, use unsigned 16bit ints
		rgb = numtype(rgb) ? 0 : rgb
		rgb = round(rgb)
		Redimension/U/W rgb
	endif
	return rgb
End



Function/T RGBA2name(r,g,b,a,maxVal)// returns the name that is closest to the given rgba
	Variable r,g,b,a						// normalized to a maximum of maxVal (not a)
	Variable maxVal

	Wave/T colorNames=ColorNamesWave()
	Wave colorRGB=ColorRGBsWave()

	Variable N=DimSize(colorRGB,0), i
	Make/N=3/FREE/D rgb0={r,g,b}
	rgb0 = round(rgb0/maxVal*255)	// renormalize to [0,255] range
	a = limit(2*a,0,1)					// the 2 reduces the effect of a
	rgb0 += (1-a)*255
	rgb0 = limit(rgb0,0,255)
	Make/N=(N,3)/FREE rgbs
	rgbs = rgb0[q]
	MatrixOP/FREE diffs = sumRows(magSqr(rgbs - colorRGB))	// this gives magsqr(colorRGB-rgbi)[i]
	WaveStats/Q/M=1 diffs
	return colorNames[V_minloc]
End


Function/WAVE color2RGB(color,maxVal)		// returns a wave with the rgb of a color
	String color									// name of color, e.g. "red", case INsenstive
	Variable maxVal								// desired scaling, e.g. red = {maxVal,0,0}

	Wave/T colorNames=ColorNamesWave()
	Wave colorRGB=ColorRGBsWave()

	color = ReplaceString(" ",color,"")			// remove spaces from color
	strswitch(color)							// deal with synonyms
		case "magenta":
		case "fuchsia":
			color = "fuchsia/magenta"
			break
		case "cyan":
		case "aqua":
			color = "cyan/aqua"
			break
	endswitch

	Variable i=-1,N=DimSize(colorNames,0)
	do
		i += 1
	while(i<N && !stringmatch(ReplaceString(" ",colorNames[i],""),color))
	Make/N=3/D/FREE rgb=NaN				// default color is {NaN,NaN,NaN}
	if (i<N)
		rgb = colorRGB[i][p] * (maxVal/255)	// if color was found, set rgb scaled to maxVal
	endif
	return rgb
End


Static Function/WAVE ColorNamesWave()
	Make/N=(554)/T/FREE colorNames=""

	colorNames[0]= {"indian red","crimson","lightpink","lightpink 1","lightpink 2","lightpink 3","lightpink 4","pink","pink 1","pink 2","pink 3","pink 4","palevioletred","palevioletred 1","palevioletred 2"}
	colorNames[15]= {"palevioletred 3","palevioletred 4","lavenderblush","lavenderblush 2","lavenderblush 3","lavenderblush 4","violetred 1","violetred 2","violetred 3","violetred 4","hotpink","hotpink 1","hotpink 2"}
	colorNames[28]= {"hotpink 3","hotpink 4","raspberry","deeppink","deeppink 2","deeppink 3","deeppink 4","maroon 1","maroon 2","maroon 3","maroon 4","mediumvioletred","violetred","orchid","orchid 1","orchid 2"}
	colorNames[44]= {"orchid 3","orchid 4","thistle","thistle 1","thistle 2","thistle 3","thistle 4","plum 1","plum 2","plum 3","plum 4","plum","violet","fuchsia / magenta","magenta 2","magenta 3","darkmagenta","purple"}
	colorNames[62]= {"mediumorchid","mediumorchid 1","mediumorchid 2","mediumorchid 3","mediumorchid 4","darkviolet","darkorchid","darkorchid 1","darkorchid 2","darkorchid 3","darkorchid 4","indigo","blueviolet"}
	colorNames[75]= {"purple 1","purple 2","purple 3","purple 4","mediumpurple","mediumpurple 1","mediumpurple 2","mediumpurple 3","mediumpurple 4","darkslateblue","lightslateblue","mediumslateblue","slateblue"}
	colorNames[88]= {"slateblue 1","slateblue 2","slateblue 3","slateblue 4","ghostwhite","lavender","blue","blue 2","mediumblue","darkblue","navy","midnightblue","cobalt","royalblue","royalblue 1","royalblue 2"}
	colorNames[104]= {"royalblue 3","royalblue 4","cornflowerblue","lightsteelblue","lightsteelblue 1","lightsteelblue 2","lightsteelblue 3","lightsteelblue 4","lightslategray","slategray","slategray 1","slategray 2"}
	colorNames[116]= {"slategray 3","slategray 4","dodgerblue","dodgerblue 2","dodgerblue 3","dodgerblue 4","aliceblue","steelblue","steelblue 1","steelblue 2","steelblue 3","steelblue 4","lightskyblue","lightskyblue 1"}
	colorNames[130]= {"lightskyblue 2","lightskyblue 3","lightskyblue 4","skyblue 1","skyblue 2","skyblue 3","skyblue 4","skyblue","deepskyblue","deepskyblue 2","deepskyblue 3","deepskyblue 4","peacock","lightblue"}
	colorNames[144]= {"lightblue 1","lightblue 2","lightblue 3","lightblue 4","powderblue","cadetblue 1","cadetblue 2","cadetblue 3","cadetblue 4","turquoise 1","turquoise 2","turquoise 3","turquoise 4","cadetblue"}
	colorNames[158]= {"darkturquoise","azure","azure 2","azure 3","azure 4","lightcyan","lightcyan 2","lightcyan 3","lightcyan 4","paleturquoise 1","paleturquoise","paleturquoise 3","paleturquoise 4","darkslategray"}
	colorNames[172]= {"darkslategray 1","darkslategray 2","darkslategray 3","darkslategray 4","cyan / aqua","cyan 2","cyan 3","darkcyan","teal","mediumturquoise","lightseagreen","manganeseblue","turquoise"}
	colorNames[185]= {"coldgrey","turquoiseblue","aquamarine","aquamarine 2","mediumaquamarine","aquamarine 4","mediumspringgreen","mintcream","springgreen","springgreen 1","springgreen 2","springgreen 3","mediumseagreen"}
	colorNames[198]= {"seagreen 1","seagreen 2","seagreen 3","seagreen","emeraldgreen","mint","cobaltgreen","honeydew","honeydew 2","honeydew 3","honeydew 4","darkseagreen","darkseagreen 1","darkseagreen 2"}
	colorNames[212]= {"darkseagreen 3","darkseagreen 4","palegreen","palegreen 1","lightgreen","palegreen 3","palegreen 4","limegreen","forestgreen","lime","green 2","green 3","green 4","green","darkgreen"}
	colorNames[227]= {"sapgreen","lawngreen","chartreuse","chartreuse 2","chartreuse 3","chartreuse 4","greenyellow","darkolivegreen 1","darkolivegreen 2","darkolivegreen 3","darkolivegreen 4","darkolivegreen"}
	colorNames[239]= {"olivedrab","olivedrab 1","olivedrab 2","yellowgreen","olivedrab 4","ivory","ivory 2","ivory 3","ivory 4","beige","lightyellow","lightyellow 2","lightyellow 3","lightyellow 4","lightgoldenrodyellow"}
	colorNames[254]= {"yellow","yellow 2","yellow 3","yellow 4","warmgrey","olive","darkkhaki","khaki 1","khaki 2","khaki 3","khaki 4","khaki","palegoldenrod","lemonchiffon","lemonchiffon 2","lemonchiffon 3"}
	colorNames[270]= {"lemonchiffon 4","lightgoldenrod 1","lightgoldenrod 2","lightgoldenrod 3","lightgoldenrod 4","banana","gold","gold 2","gold 3","gold 4","cornsilk","cornsilk 2","cornsilk 3","cornsilk 4"}
	colorNames[284]= {"goldenrod","goldenrod 1","goldenrod 2","goldenrod 3","goldenrod 4","darkgoldenrod","darkgoldenrod 1","darkgoldenrod 2","darkgoldenrod 3","darkgoldenrod 4","orange","orange 2","orange 3"}
	colorNames[297]= {"orange 4","floralwhite","oldlace","wheat","wheat 1","wheat 2","wheat 3","wheat 4","moccasin","papayawhip","blanchedalmond","navajowhite","navajowhite 2","navajowhite 3","navajowhite 4"}
	colorNames[312]= {"eggshell","tan","brick","cadmiumyellow","antiquewhite","antiquewhite 1","antiquewhite 2","antiquewhite 3","antiquewhite 4","burlywood","burlywood 1","burlywood 2","burlywood 3","burlywood 4"}
	colorNames[326]= {"bisque","bisque 2","bisque 3","bisque 4","melon","carrot","darkorange","darkorange 1","darkorange 2","darkorange 3","darkorange 4","orange","tan 1","tan 2","peru","tan 4","linen","peachpuff"}
	colorNames[344]= {"peachpuff 2","peachpuff 3","peachpuff 4","seashell","seashell 2","seashell 3","seashell 4","sandybrown","rawsienna","chocolate","chocolate 1","chocolate 2","chocolate 3","saddlebrown"}
	colorNames[358]= {"ivoryblack","flesh","cadmiumorange","burntsienna","sienna","sienna 1","sienna 2","sienna 3","sienna 4","lightsalmon","lightsalmon 2","lightsalmon 3","lightsalmon 4","coral","orangered"}
	colorNames[373]= {"orangered 2","orangered 3","orangered 4","sepia","darksalmon","salmon 1","salmon 2","salmon 3","salmon 4","coral 1","coral 2","coral 3","coral 4","burntumber","tomato","tomato 2","tomato 3"}
	colorNames[390]= {"tomato 4","salmon","mistyrose","mistyrose 2","mistyrose 3","mistyrose 4","snow","snow 2","snow 3","snow 4","rosybrown","rosybrown 1","rosybrown 2","rosybrown 3","rosybrown 4","lightcoral"}
	colorNames[406]= {"indianred","indianred 1","indianred 2","indianred 4","indianred 3","brown","brown 1","brown 2","brown 3","brown 4","firebrick","firebrick 1","firebrick 2","firebrick 3","firebrick 4"}
	colorNames[421]= {"red","red 2","red 3","darkred","maroon","sgi beet","sgi slateblue","sgi lightblue","sgi teal","sgi chartreuse","sgi olivedrab","sgi brightgray","sgi salmon","sgi darkgray","sgi gray 12"}
	colorNames[436]= {"sgi gray 16","sgi gray 32","sgi gray 36","sgi gray 52","sgi gray 56","sgi lightgray","sgi gray 72","sgi gray 76","sgi gray 92","sgi gray 96","white","gray 96","gainsboro","lightgrey"}
	colorNames[450]= {"silver","darkgray","gray","gray 42","black","gray 99","gray 98","gray 97","gray 96","gray 95","gray 94","gray 93","gray 92","gray 91","gray 90","gray 89","gray 88","gray 87","gray 86"}
	colorNames[469]= {"gray 85","gray 84","gray 83","gray 82","gray 81","gray 80","gray 79","gray 78","gray 77","gray 76","gray 75","gray 74","gray 73","gray 72","gray 71","gray 70","gray 69","gray 68","gray 67"}
	colorNames[488]= {"gray 66","gray 65","gray 64","gray 63","gray 62","gray 61","gray 60","gray 59","gray 58","gray 57","gray 56","gray 55","gray 54","gray 53","gray 52","gray 51","gray 50","gray 49","gray 48"}
	colorNames[507]= {"gray 47","gray 46","gray 45","gray 44","gray 43","gray 42","gray 42","gray 40","gray 39","gray 38","gray 37","gray 36","gray 35","gray 34","gray 33","gray 32","gray 31","gray 30","gray 29"}
	colorNames[526]= {"gray 28","gray 27","gray 26","gray 25","gray 24","gray 23","gray 22","gray 21","gray 20","gray 19","gray 18","gray 17","gray 16","gray 15","gray 14","gray 13","gray 12","gray 11","gray 10"}
	colorNames[545]= {"gray 9","gray 8","gray 7","gray 6","gray 5","gray 4","gray 3","gray 2","gray 1"}
	return colorNames
End

Static Function/WAVE ColorRGBsWave()
	Make/N=(554,3)/B/U/FREE colorRGB=0

	colorRGB[0][0]= {176,220,255,255,238,205,139,255,255,238,205,139,219,255,238,205,139,255,238,205,139,255,238,205,139,255,255,238,205,139,135,255,238,205,139,255,238,205,139,199,208,218,255,238,205,139,216}
	colorRGB[47][0]= {255,238,205,139,255,238,205,139,221,238,255,238,205,139,128,186,224,209,180,122,148,153,191,178,154,104,75,138,155,145,125,85,147,171,159,137,93,72,132,123,106,131,122,105,71,248,230,0}
	colorRGB[95][0]= {0,0,0,0,25,61,65,72,67,58,39,100,176,202,188,162,110,119,112,198,185,159,108,30,28,24,16,240,70,99,92,79,54,135,176,164,141,96,135,126,108,74,135,0,0,0,0,51,173,191,178,154,104,176,152}
	colorRGB[150][0]= {142,122,83,0,0,0,0,95,0,240,224,193,131,224,209,180,122,187,174,150,102,47,151,141,121,82,0,0,0,0,0,72,32,3,64,128,0,127,118,102,69,0,245,0,0,0,0,60,84,78,67,46,0,189,61,240,224,193,131}
	colorRGB[209][0]= {143,193,180,155,105,152,154,144,124,84,50,34,0,0,0,0,0,0,48,124,127,118,102,69,173,202,188,162,110,85,107,192,179,154,105,255,238,205,139,245,255,238,205,139,250,255,238,205,139,128,128}
	colorRGB[260][0]= {189,255,238,205,139,240,238,255,238,205,139,255,238,205,139,227,255,238,205,139,255,238,205,139,218,255,238,205,139,184,255,238,205,139,255,238,205,139,255,253,245,255,238,205,139,255}
	colorRGB[306][0]= {255,255,255,238,205,139,252,210,156,255,250,255,238,205,139,222,255,238,205,139,255,238,205,139,227,237,255,255,238,205,139,255,255,238,205,139,250,255,238,205,139,255,238,205,139,244}
	colorRGB[352][0]= {199,210,255,238,205,139,41,255,255,138,160,255,238,205,139,255,238,205,139,255,255,238,205,139,94,233,255,238,205,139,255,238,205,139,138,255,238,205,139,250,255,238,205,139,255,238,205}
	colorRGB[399][0]= {139,188,255,238,205,139,240,205,255,238,139,205,165,255,238,205,139,178,255,238,205,139,255,238,205,139,128,142,113,125,56,113,142,197,198,85,30,40,81,91,132,142,170,183,193,234,244,255}
	colorRGB[447][0]= {245,220,211,192,169,128,105,0,252,250,247,245,242,240,237,235,232,229,227,224,222,219,217,214,212,209,207,204,201,199,196,194,191,189,186,184,181,179,176,173,171,168,166,163,161,158,156}
	colorRGB[494][0]= {153,150,148,145,143,140,138,135,133,130,127,125,122,120,117,115,112,110,107,105,102,99,97,94,92,89,87,84,82,79,77,74,71,69,66,64,61,59,56,54,51,48,46,43,41,38,36,33,31,28,26,23,20,18}
	colorRGB[548][0]= {15,13,10,8,5,3}
	colorRGB[0][1]= {23,20,182,174,162,140,95,192,181,169,145,99,112,130,121,104,71,240,224,193,131,62,58,50,34,105,110,106,96,58,38,20,18,16,10,52,48,41,28,21,32,112,131,122,105,71,191,225,210,181,123,187}
	colorRGB[52][1]= {174,150,102,160,130,0,0,0,0,0,85,102,95,82,55,0,50,62,58,50,34,0,43,48,44,38,26,112,130,121,104,71,61,112,104,90,111,103,89,60,248,230,0,0,0,0,0,25,89,105,118,110,95,64,149,196,225,210}
	colorRGB[110][1]= {181,123,136,128,226,211,182,123,144,134,116,78,248,130,184,172,148,100,206,226,211,182,123,206,192,166,112,206,191,178,154,104,161,216,239,223,192,131,224,245,229,197,134,245,229,197}
	colorRGB[156][1]= {134,158,206,255,238,205,139,255,238,205,139,255,238,205,139,79,255,238,205,139,255,238,205,139,128,209,178,168,224,138,199,255,238,205,139,250,255,255,238,205,139,179,255,238,205,139}
	colorRGB[202][1]= {201,252,145,255,238,205,139,188,255,238,205,139,251,255,238,205,139,205,139,255,238,205,139,128,100,128,252,255,238,205,139,255,255,238,205,139,107,142,255,238,205,139,255,238,205,139}
	colorRGB[248][1]= {245,255,238,205,139,250,255,238,205,139,128,128,183,246,230,198,134,230,232,250,233,201,137,236,220,190,129,207,215,201,173,117,248,232,200,136,165,193,180,155,105,134,185,173,149,101}
	colorRGB[294][1]= {165,154,133,90,250,245,222,231,216,186,126,228,239,235,222,207,179,121,230,180,102,153,235,239,223,192,131,184,211,197,170,115,228,213,183,125,168,145,140,127,118,102,69,128,165,154,133}
	colorRGB[341][1]= {90,240,218,203,175,119,245,229,197,134,164,97,105,127,118,102,69,36,125,97,54,82,130,121,104,71,160,149,129,87,127,69,64,55,37,38,150,140,130,112,76,114,106,91,62,51,99,92,79,54,128,228}
	colorRGB[393][1]= {213,183,125,250,233,201,137,143,193,180,155,105,128,92,106,99,58,85,42,64,59,51,35,34,48,44,38,26,0,0,0,0,0,56,113,158,142,198,142,193,113,85,30,40,81,91,132,142,170,183,193,234,244,255}
	colorRGB[447][1]= {245,220,211,192,169,128,105,0,252,250,247,245,242,240,237,235,232,229,227,224,222,219,217,214,212,209,207,204,201,199,196,194,191,189,186,184,181,179,176,173,171,168,166,163,161,158,156}
	colorRGB[494][1]= {153,150,148,145,143,140,138,135,133,130,127,125,122,120,117,115,112,110,107,105,102,99,97,94,92,89,87,84,82,79,77,74,71,69,66,64,61,59,56,54,51,48,46,43,41,38,36,33,31,28,26,23,20,18}
	colorRGB[548][1]= {15,13,10,8,5,3}
	colorRGB[0][2]= {31,60,193,185,173,149,101,203,197,184,158,108,147,171,159,137,93,245,229,197,134,150,140,120,82,180,180,167,144,98,87,147,137,118,80,179,167,144,98,133,144,214,250,233,201,137,216,255,238}
	colorRGB[49][2]= {205,139,255,238,205,139,221,238,255,238,205,139,128,211,255,238,205,139,211,204,255,238,205,139,130,226,255,238,205,139,219,255,238,205,139,139,255,238,205,255,238,205,139,255,250,255}
	colorRGB[95][2]= {238,205,139,128,112,171,225,255,238,205,139,237,222,255,238,205,139,153,144,255,238,205,139,255,238,205,139,255,180,255,238,205,139,250,255,238,205,139,255,238,205,139,235,255,238,205}
	colorRGB[141][2]= {139,201,230,255,238,205,139,230,255,238,205,139,255,238,205,139,160,209,255,238,205,139,255,238,205,139,255,238,205,139,79,255,238,205,139,255,238,205,139,128,204,170,158,208,135,140}
	colorRGB[187][2]= {212,198,170,116,154,250,127,118,102,69,113,159,148,128,87,87,201,64,240,224,193,131,143,193,180,155,105,152,154,144,124,84,50,34,0,0,0,0,0,0,20,0,0,0,0,0,47,112,104,90,61,47,35,62,58}
	colorRGB[242][2]= {50,34,240,224,193,131,220,224,209,180,122,210,0,0,0,0,105,0,107,143,133,115,78,140,170,205,191,165,112,139,130,112,76,87,0,0,0,0,220,205,177,120,32,37,34,29,20,11,15,14,12,8,0,0,0,0,240}
	colorRGB[299][2]= {230,179,186,174,150,102,181,213,205,173,161,139,94,201,140,31,18,215,219,204,176,120,135,155,145,125,85,196,183,158,107,105,33,0,0,0,0,0,0,79,73,63,43,230,185,173,149,101,238,222,191}
	colorRGB[350][2]= {130,96,20,30,36,33,29,19,33,64,3,15,45,71,66,57,38,122,114,98,66,80,0,0,0,0,18,122,105,98,84,57,86,80,69,47,36,71,66,57,38,114,225,210,181,123,250,233,201,137,143,193,180,155,105,128}
	colorRGB[406][2]= {92,106,99,58,85,42,64,59,51,35,34,48,44,38,26,0,0,0,0,0,142,198,192,142,113,56,170,113,85,30,40,81,91,132,142,170,183,193,234,244,255,245,220,211,192,169,128,105,0,252,250,247,245,242}
	colorRGB[460][2]= {240,237,235,232,229,227,224,222,219,217,214,212,209,207,204,201,199,196,194,191,189,186,184,181,179,176,173,171,168,166,163,161,158,156,153,150,148,145,143,140,138,135,133,130,127,125}
	colorRGB[506][2]= {122,120,117,115,112,110,107,105,102,99,97,94,92,89,87,84,82,79,77,74,71,69,66,64,61,59,56,54,51,48,46,43,41,38,36,33,31,28,26,23,20,18,15,13,10,8,5,3}
	return colorRGB
End
