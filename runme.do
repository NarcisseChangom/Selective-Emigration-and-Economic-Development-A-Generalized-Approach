

clear all
cls
set more off
set rmsg on

* Replace the following line with the appropriate directory to designate the folder for the replication package
	gl mydir "add directory here\Replication-ier"
	cd "$mydir"
	
	gl dataSTATA "$mydir\stata\data\"
	gl dataMATLAB "$mydir\matlab\output\"
	gl FIG "$mydir\main\figures\"	
	gl TAB "$mydir\main\tables\"	



****************** Install Required Packages ***********************************

local pkg "estout ppmlhdfe hdfe ranktest ivreghdfe etime coefplot wbopendata kdens etime moremata labutil"
foreach pk of local pkg{
	cap which `pk'
	if _rc !=0 ssc install `pk', replace 
}



* esttab not in ssc
cap ado uninstall ftools
net install ftools, from("https://raw.githubusercontent.com/sergiocorreia/ftools/master/src/") replace

cap ado uninstall reghdfe
net install reghdfe, from("https://raw.githubusercontent.com/sergiocorreia/reghdfe/master/src/") replace


cap ado uninstall ivreg2
ssc install ivreg2, replace

cap ado uninstall ivreghdfe
net install ivreghdfe, from("https://raw.githubusercontent.com/sergiocorreia/ivreghdfe/master/src/") replace

etime, start

****************** Tables, Figures and Results in the main text  ***************

do "$mydir/main/main_results.do"

****************** Tables, Figures and Results in the Appendix  ****************

do "$mydir/appendix/appendix_results.do"

set rmsg off

etime // display the time duration of the program 
/* END
