#pragma rtGlobals= 2
Constant JZTalwaysFirst_Version=2.4
//#pragma version = JZTalwaysFirst_Version
//#pragma IndependentModule=JZTalwaysFirst
#include "GeneralFirst", version>=3.0
#include "LaueGoFirst", version>=1.0
#pragma ModuleName=JZTalwaysFirst
#pragma hide = 1


Static Function AfterFileOpenHook(refNum,file,pathName,type,creator,kind)
	Variable refNum, kind
	String file,pathName,type,creator
	if ((kind==1) || (kind==2))		// an experiment (packed or unpacked)
		PathInfo $pathName				// expand the path name, "/Users/name/data/Copper" is better than "home'
		pathName = SelectString(V_flag,pathName+":",S_path)
		printf "\r%s  %s  restarting this file on '%s' from '%s%s'\r\r",date(),time(),getHostName(1),pathName,file
	endif
End
Static Function IgorStartOrNewHook(IgorApplicationNameStr)
	String IgorApplicationNameStr
	printf "%s  %s  starting new Untitled Igor Experiment on '%s'\r\r",date(),time(),getHostName(1)
	//	GetWindow kwCmdHist wsize
	//	V_right = (V_right>900) ? V_right-100 : V_right
	//	MoveWindow /C V_left+20,V_top-40,V_right,V_bottom-2
End

Static Function/T getHostName(local)	// return unix hostname
	Variable local							// 0-> network name (bob.world.com), 1-> local name (Bob's MacBook)
	if (!stringmatch(StringByKey("OS",IgorInfo(3)),"Macintosh OS X"))
		print "Only know how to get user name from a Mac"
		return ""							// cannot get answer
	endif
	if (local)								// get the name in double quotes
		ExecuteScriptText "do shell script \"scutil --get ComputerName\""	// gets local computer name
	else									// get the network name (hostname) in double quotes
		ExecuteScriptText "do shell script \"hostname\""
	endif
	Variable i=strlen(S_value)-2
	if (i<=0)								// no name,
		DoAlert 0, "Unable to get user name, message is '"+S_value+"'"
	endif
	return S_value[1,i]					// first & last char are double-quotes, remove
End
