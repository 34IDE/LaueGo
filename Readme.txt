LaueGo - X-ray Microdiffraction (XMD) data analysis package based on Igor Pro.

by Beamline 34ID-E, Advanced Photon Source.

=======================
  Installation Guide
(Last updated 12/2018)
=======================
The following applies to both MS Windows and Mac computers. "[MyDocuments]" refers to the special user documents folder designated by the operating system, namely "%USERPROFILE%\My documents" on Windows, or "$HOME:Documents" on Mac.

1.	If you've just freshly installed Igor Pro, run it at least once. This will create a "Wavemetrics" folder inside "[MyDocuments]".

2.	Download and unzip the LaueGo package from 34ID-E website.
	https://sector33.xray.aps.anl.gov/~tischler/

3.	unzip (if necessary), and you should now have a foler named "LaueGo" somewhere convienent. It should be named "LaueGo" (change if necessary)



for WINDOWS:
"[MyDocuments]" refers to the special user documents folder designated by the operating system, namely "%USERPROFILE%\My documents" on Windows

4.	Copy the unzipped "LaueGo" folder into "[MyDocuments]\WaveMetrics\Igor Pro 8 User Files\User Procedures"; if an old version of the "LaueGo" folder exists in "User Procedures", (optionally backup) and replace it. Inside the "User Procedures" folder you only have ONE copy of the "LaueGo" procedures.

5.	(First-time installation only) Inside the "User Procedures\LaueGo" folder,  find the file "always\always first.ipf". Create a shortcut to this file and put it into "[MyDocuments]\WaveMetrics\Igor Pro 8 User Files\Igor Procedures".

6.	(First-time installation only) Find the files:
"C:\Program Files (x86)\WaveMetrics\Igor Pro Folder\More Extensions\File Loaders\HDF5.xop"
and
"C:\Program Files (x86)\WaveMetrics\Igor Pro Folder\More Extensions\File Loaders\HDF5 Help.ihf"
Create shortcuts to both these files and put the aliases inside "[MyDocuments]\WaveMetrics\Igor Pro 8 User Files\Igor Extensions (64-bit)".
(for Igor 7, you can use the 32 bit, but 64 is best, for Igor6, use "Igor Extensions")
7.	Done. You only need to repeat step 2 and 3 the next time you update your microdiffraction package.



for MAC:
the folder "Documents:..." will be in your home folder. You can open "Documents" from the Finder "Go" menu.

4.	Copy the unzipped "LaueGo" folder into:
"Documents:WaveMetrics:Igor Pro 8 User Files:User Procedures"
If an old version of the "LaueGo" folder exists in "User Procedures", (optionally backup) and replace it. Inside the "User Procedures" folder you only have ONE copy of the "LaueGo" procedures.

5.	(First-time installation only) Inside the "User Procedures:LaueGo" folder,  
find the file "always:always first.ipf". Create a alias to this file and put the alias into "Documents:WaveMetrics:Igor Pro 8 User Files:Igor Procedures".

6.	(First-time installation only) Find the files:
Applications:Igor Pro 8 Folder:More Extensions (64-bit):File Loaders/HDF5-64.xop
and
Applications:Igor Pro 8 Folder:More Extensions:File Loaders:HDF5 Help.ihf 
Create shortcuts to both these files and put the aliases inside "Documents:WaveMetrics"Igor Pro 8 User Files:Igor Extensions (64-bit)".

7.	Done. You only need to repeat step 2 and 3 the next time you update your microdiffraction package.
