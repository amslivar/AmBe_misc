
This is a repo with several pieces of code used to calibrate and work with data collected from the Ge detectors at UA, and to complement the other AmBe directories. 

# Convert files 
- Convert the *.Chn files to *.dat files using the get_asci_alpha executable attached. The script should work on neutrino3 if you just copy paste it to your folder.
- *.Chn files are saved under neutrino3.ph.ua.edu:/mnt/raid5/raid4b/users/labge/GeIII/<year>/<folderName>
- make a list of files named as “in_read” using command:
```
   ls -1 <path/to/ChnFiles>/*.Chn > in_read
```
then run get_asci_alpha script placing it in the same directory as “in_read”: 
```
./get_asci_alpha
```
- Once you have *.dat files you can use them as an input to extract activity from their data.

- (A known issue is that there is a cap on the number of characters allowed on the file name string. A workaround is to copy the files in your folder momentarily for a shorter path)

# Root script for GeI, GeII, or GeIII calibration.

Steps:
- Convert the *.Chn file from MAESTRO into *.dat file using the “get_asci_alpha” executable described above.
- The calibrate_Ge.C script takes that *. dat file (for example Th228_5min.dat)  as an input, It also asks for the calibration source name and detector name as an input. The output pdf contains the calibration result will be named after the detector name and source name. Open root in the directory where the code is, then run the code using the format:

```
calibrate_Ge("Th228_Cal.dat","GeIII","Th228")
```

- The script use the source name to search for corresponding gamma peaks and fit them to extract their gaussian mean and width.
    - You may need to adjust the fit range or modify the SetParamters(….) function to get a good fit.
    - To get the approximate locations for the corresponding energy peaks the fit uses ADC-Energy maps from earlier calibration. Those values does’t usually need to be adjusted unless there is a major change in the detector such as change in the detector gain due to the usage of attenuator. 
- This script can also be used for Ba133 calibration as well. Just need to provide the source name.
