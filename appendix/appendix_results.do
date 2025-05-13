

********************************************************************************
*		Tables and figures appendix 
********************************************************************************
	clear all
	gl mydir "U:\CHANGOM\CDDM2021_Completefile\CHAP1\Replication-ier"
	cd "$mydir"
	
	gl dataSTATA "$mydir\stata\data\"
	gl dataMATLAB "$mydir\matlab\output\"
	gl FIG "$mydir\appendix\figures\"	
	gl TAB "$mydir\appendix\tables\"	

*%%%------------------------------------------------------------------------%%%*
*%%%							A - TABLES 									%%%*
*%%%------------------------------------------------------------------------%%%*
	
	
/*------------------------------------------------------------------------------
	Table A.1 : Pseudo gravity regression for predicting bilateral migration stocks
*------------------------------------------------------------------------------*/

use "${dataSTATA}data_gravity", clear 

	g logpop =ln(population)
	g logdist =ln(distw)
	g logNetw_10 =ln(1+MIG_10)
	g logNetw_20 =ln(1+MIG_20)
	g logwj =ln(wj)
	g loggdpc_o = ln(GDPc)
	g loggdpc_d = ln(GDPc_d)		
	
	local mig NijL NijH
	local dyadic_cntrs logNetw_20 comlang_off colony contig 
	local orig_cntrs logpop 
	local dest_cntrs logwj
	local fes id_d year 
	local cls id_d 

	lab var logpop "Population size origin (log)"
	lab var logdist "Distance (log)"	
	lab var logwj "Wage at destination (log)"
	lab var logNetw_20 "Network 20 years ago (log)"	
	lab var comlang_off "Common language"		
	lab var colony "Colony"	
	lab var contig "Contiguity"
	

	
	eststo clear 
	eststo: quiet  ppmlhdfe NijH `orig_cntrs' logdist `dyadic_cntrs' `dest_cntrs', abs(`fes') cl(`cls') d 
	eststo m1
	predict NijH_hatT20 if NijH >0, mu
	estadd local dFE "\checkmark", replace
	estadd local tFE "\checkmark", replace
	estadd local R2 = `e(r2_p)' 
	
	eststo: quiet  ppmlhdfe NijL `orig_cntrs' logdist `dyadic_cntrs' `dest_cntrs', abs(`fes') cl(`cls') d 
	eststo m2 
	predict NijL_hatT20 if NijL >0, mu
	estadd local dFE "\checkmark", replace
	estadd local tFE "\checkmark", replace
	estadd local R2 = `e(r2_p)' 
	
	eststo: quiet  ppmlhdfe NijH `orig_cntrs' logdist i.year#c.logdist `dyadic_cntrs' `dest_cntrs', abs(`fes') cl(`cls') d 
	eststo m3
	predict NijH_hatTD20 if NijH >0, mu
	estadd local dFE "\checkmark", replace
	estadd local tFE "\checkmark", replace
	estadd local R2 = `e(r2_p)' 
	
	eststo: quiet  ppmlhdfe NijL `orig_cntrs' logdist i.year#c.logdist `dyadic_cntrs' `dest_cntrs', abs(`fes') cl(`cls') d 
	eststo m4 
	predict NijL_hatTD20 if NijL >0, mu
	estadd local dFE "\checkmark", replace
	estadd local tFE "\checkmark", replace
	estadd local R2 = `e(r2_p)'
	
	esttab m1 m2 m3 m4, label nomtitles se star(* 0.10 ** 0.05 *** 0.01) ///
		   drop(1990.year#c.logdist _cons) ///
		   order(logpop logdist 2000.year#c.logdist 2010.year#c.logdist logwj logNetw_20) ///
		   stats(N R2 dFE tFE, fmt(%9.0fc %9.1fc  %9.1fc %9.1fc) ///
			  label("Observations" "Pseudo R-sq" "Destination FE" "Time FE"))

esttab m1 m2 m3 m4 using  "$TAB\Table_A1.tex", ///
			  prehead(\begin{tabular}{l*{@M}{c}} \toprule) ///
			  b(3) se(3) drop(1990.year#c.logdist _cons)  ///
			  star(* 0.10 ** 0.05 *** 0.01) ///
			  label nonotes nomtitle collabels(none) compress ///
			  stats(N R2 dFE tFE, fmt(%9.0fc %9.3fc %9.0fc) ///
			  label("Observations" "Pseudo R-sq" "Destination FE" "Time FE")) ///
			  title("Pseudo-Gravity model for dyadic migration stocks (M_{ijs})") booktabs ///
			  addnotes("Pseudo Maximum Likelihood (PPML) coefficients with standard errors clustered at oigin-destination pair in parentheses. \sym{***}, \sym{**}, and \sym{*} denote significance at the 1\%, 5\% and 10\% level respectively.") ///
			  mgroups("M_{ijht}" "M_{ijlt}" "M_{ijht}" "M_{ijlt}", ///
			  pattern(1 1 1 1) ///
			  prefix(\multicolumn{@span}{c}{) suffix(}) span ///
			  erepeat(\cmidrule(lr){@span})) alignment(D{.}{.}{-1}) replace	




/*------------------------------------------------------------------------------
	Table A.2 : First-stage regression (instrumenting emigration differentials in 1990 and 2000)
*-------------------------------------------------------------------------------*/

use "${dataSTATA}data_empirics", clear 

	xtset id year, delta(10)		
		g dlH =ln(h_n)-ln(l.h_n)
		g lH_a =ln(l.h_n)
		g p =m_h-m_l
		g p_hat =mhhat20 - mlhat20
		g income_gp=1 if inrange(GDPc,0,1000) // Income groups based on GDPpc
		replace income_gp=2 if inrange(GDPc,1001,3999) 
		replace income_gp=3 if inrange(GDPc,4000,12499) 
		replace income_gp=4 if GDPc>=12500 
		g dvping = (inlist(income_gp,1,2,3))
		g relgap = p/p_hat

		gen admit=0 // Identifying outliers 
		replace admit=1 if inlist(id,106,110,116,68,134,100,101,90) 
		replace admit=1 if inlist(id,37,121,2)


** Variable definition and macros 

		g l_p = l.p
		g l_p_hat  = l.p_hat 
		g l_DEN = l.DEN 
		gl controls lH_a l_DEN
		gl fes income_gp year
		gl endog l_p c.l_p#income_gp 
		gl instr l_p_hat c.l_p_hat#income_gp 
		gl yvar dlH
		gl fsvar l_p
			
		g sample1 = 1 // full sample 
		g sample2 = 1 if dvping==1 // developing countries 
		
		lab var lH_a "$ ln(H_{i,t}) $"
		lab var l_p "$ m_{i,h,t}-m_{i,l,t} \equiv \delta_{i,t} $"
		lab var l_p_hat "$ \widehat{\delta_{i,t}} $"
		lab var l_DEN "Population density"
		lab define income_gp 1 "Low-Income" 2 "Lower-Middle" 3 "Upper-Middle" 4 "High-Income"
		lab values income_gp income_gp 
		lab var admit "Outliers"

		

***  First stage regression per se 

eststo clear
local tex_file_name "Table_A2.tex"
local latex_header ///
	" & \multicolumn{7}{c}{Emigration differential} \\ \cmidrule(lr){2-7} & \multicolumn{3}{c}{Developing countries} & \multicolumn{3}{c}{Full sample} \\ \cmidrule(lr){2-4} \cmidrule(lr){5-7} & (1) & (2) & (3) & (4) & (5) & (6) \\" 
eststo fs1: quiet reghdfe $fsvar  $instr $controls if sample2==1, abs($fes) vce(cl id)
	estadd local dFE "\checkmark", replace
	estadd local tFE "\checkmark", replace
	estadd local oFE "-", replace
	estadd local out "Included", replace
	
eststo fs2: quiet reghdfe $fsvar  $instr $controls admit if sample2==1, abs($fes) vce(cl id)
	estadd local dFE "\checkmark", replace
	estadd local tFE "\checkmark", replace
	estadd local oFE "\checkmark", replace
	estadd local out "Included", replace
	
eststo fs3: quiet reghdfe $fsvar  $instr $controls if sample2==1 & admit!=1, abs($fes) vce(cl id)
	estadd local dFE "\checkmark", replace
	estadd local tFE "\checkmark", replace
	estadd local oFE "-", replace
	estadd local out "Excluded", replace
	
eststo fs4: quiet reghdfe $fsvar  $instr $controls if sample1==1, abs($fes) vce(cl id)
	estadd local dFE "\checkmark", replace
	estadd local tFE "\checkmark", replace
	estadd local oFE "-", replace
	estadd local out "Included", replace
	
eststo fs5: quiet reghdfe $fsvar  $instr $controls admit if sample1==1, abs($fes) vce(cl id)
	estadd local dFE "\checkmark", replace
	estadd local tFE "\checkmark", replace
	estadd local oFE "\checkmark", replace
	estadd local out "Included", replace
	
eststo fs6: quiet reghdfe $fsvar  $instr $controls if sample1==1 & admit!=1, abs($fes) vce(cl id)
	estadd local dFE "\checkmark", replace
	estadd local tFE "\checkmark", replace
	estadd local oFE "-", replace
	estadd local out "Excluded", replace
	

estout fs* ///
	using "${TAB}/`tex_file_name'", replace ///
	drop(_cons 1.income_gp#c.l_p_hat admit) ///
	order(lH_a l_p_hat 2.income_gp#c.l_p_hat 3.income_gp#c.l_p_hat 4.income_gp#c.l_p_hat l_DEN) ///
	style(tex) mlabel(none) collabel(none) label ///
	prehead("\begin{tabular}{lcccccc}" "\toprule") ///
	posthead("`latex_header'" "\midrule") ///
	cells("b(fmt(%9.3f) star)" "se(fmt(%9.3f) par)")  ///
	starlevels(* 0.10 ** 0.05 *** 0.01) ///
	prefoot("\midrule") ///
	stats(N r2 dFE tFE oFE out, labels("N" "\$R^2\$" "Income group FE" "Decade FE" "Outlier FE" "Outlier") ///
	fmt(%9.0fc %10.3f %9.0g %9.0g  %9.0g  %9.0g)) ///
	postfoot("\bottomrule" "\end{tabular}")	




/*------------------------------------------------------------------------------
	Table B.1 : Validating the calibrated migration costs
*------------------------------------------------------------------------------*/
	use "${dataSTATA}migcost_validation", clear 

	
** Index of country skill selectivity (leave on out approach)	

	foreach mig in NijH_00 NijL_00 {
		bys isod: egen I`mig' = total(`mig')
	}
	
	g N_ijH = INijH_00 - NijH_00
	g N_ijL = INijL_00 - NijL_00 
	egen N_ij = rowtotal(N_ijH N_ijL)
	g selec_ind = N_ijH/N_ij
	g selec_ind2= 100*selec_ind
	
	g logdist =ln(distw)
	g logNetw =ln(1+MIG90)
	g logselec =ln(selec_ind)

	lab var logdist "Distance (log)"	
	lab var logNetw "Network 20 years ago (log)"	
	lab var comlang_off "Common language"		
	lab var colony "Colonial link"
	lab var logselec " Skill selection index 10 years ago (log)"	
	lab var selec_ind "selection index"	
	lab var selec_ind2 "selection index"	
	lab var avgVisaall "Visa requirements"	
	lab var GUEST_90 "Guestwork prog. 90s"	
	
	
	local contr1 logdist comlang_off colony logNetw
	local contr2 avgVisaall GUEST_90 
	local contr3 logselec 

*** Validation exercise per se 
	eststo clear
	local tex_file_name "Table_B1.tex"
	local latex_header ///
	" & \multicolumn{3}{c}{High skilled migration costs ($c_{ijh}$)} & \multicolumn{3}{c}{Low skilled migration costs ($c_{ijl}$)} \\ \cmidrule(lr){2-4} \cmidrule(lr){5-7} & (1) & (2) & (3) & (4) & (5) & (6) \\ \cmidrule(lr){2-4} \cmidrule(lr){5-7}" 
	eststo m1: quiet reghdfe cijha `contr1', abs(isoo isod) vce(rob)
			estadd local oFE "\checkmark", replace
			estadd local dFE "\checkmark", replace
	eststo m2: quiet reghdfe cijha `contr1' `contr2', abs(isoo isod) vce(rob)
			estadd local oFE "\checkmark", replace
			estadd local dFE "\checkmark", replace	
	eststo m3: quiet reghdfe cijha `contr1' `contr2' `contr3', abs(isoo isod) vce(rob)
			estadd local oFE "\checkmark", replace
			estadd local dFE "\checkmark", replace
	eststo m4: quiet reghdfe cijla `contr1', abs(isoo isod) vce(rob)
			estadd local oFE "\checkmark", replace
			estadd local dFE "\checkmark", replace
	eststo m5: quiet reghdfe cijla `contr1' `contr2', abs(isoo isod) vce(rob)
			estadd local oFE "\checkmark", replace
			estadd local dFE "\checkmark", replace
	eststo m6: quiet reghdfe cijla `contr1' `contr2' `contr3', abs(isoo isod) vce(rob)
			estadd local oFE "\checkmark", replace
			estadd local dFE "\checkmark", replace
		
estout m1 m2 m3 m4 m5 m6 ///
	using "${TAB}/`tex_file_name'", replace ///
	style(tex) mlabel(none) collabel(none) label ///
	prehead("\begin{tabular}{lcccccc}" "\toprule\toprule") ///
	posthead("`latex_header'" "\midrule") ///
	cells("b(fmt(%9.3f) star)" "se(fmt(%9.3f) par)")  ///
	starlevels(* 0.10 ** 0.05 *** 0.01) ///
	prefoot("\midrule") ///
	stats(N r2 oFE dFE, labels("N" "\$R^2\$" "Origin country FE" "Destination country FE") ///
	fmt(%9.0fc %10.3f %9.0g %9.0g  %9.0g  %9.0g)) ///
	postfoot("\bottomrule" "\end{tabular}")	



/*------------------------------------------------------------------------------
	Table B.2 : Validating the scale factor (Gi) of Access to education 
*------------------------------------------------------------------------------*/


***	Data preparation 

	
* Controls access in January 04, 2025 (data subject to updates)
wbopendata, indicator(NY.GDP.PCAP.PP.KD; SP.POP.TOTL;SP.POP.1564.TO;BX.TRF.PWKR.CD.DT; ///
					  NY.GDP.PCAP.KD;NY.GNP.PCAP.PP.KD; NY.GNP.PCAP.KD; ///
					  SP.URB.TOTL.IN.ZS;SE.XPD.TOTL.GD.ZS) ///
					  year(2000:2020) clear long  
drop if incomelevel == "NA" // Get rid of regional aggregates 
ta countrycode if missing(ny_gdp_pcap_pp_kd)		  
g iso3 = countrycode				  

g included =0
replace included = 1 if inlist(iso3,"AFG","AGO","ALB","ARE","ARG","ARM","AUS","AUT")
replace included = 1 if inlist(iso3,"AZE","BDI","BEL","BEN","BFA","BGR","BHR","BHS")					  
replace included = 1 if inlist(iso3,"BIH","BLR","BLZ","BOL","BRA","BRB","BRN","BTN")
replace included = 1 if inlist(iso3,"BWA","CAF","CAN","CHL","CHN","CIV","CMR","COG")
replace included = 1 if inlist(iso3,"COL","COM","CPV","CRI","CUB","CYP","CZE","DEU")
replace included = 1 if inlist(iso3,"DJI","DNK","DOM","DZA","ECU","EGY","ERI","ESP")
replace included = 1 if inlist(iso3,"EST","ETH","FIN","FJI","FRA","FSM","GAB","GBR")
replace included = 1 if inlist(iso3,"GEO","GHA","GIN","GMB","GNB","GNQ","GRC","GRD")
replace included = 1 if inlist(iso3,"GTM","GUY","HND","HRV","HTI","HUN","IDN","IND")
replace included = 1 if inlist(iso3,"IRL","IRN","IRQ","ISL","ISR","ITA","JAM","JOR")
replace included = 1 if inlist(iso3,"JPN","KAZ","KEN","KGZ","KHM","KWT","LAO","LBN")
replace included = 1 if inlist(iso3,"LBR","LBY","LCA","LKA","LSO","LTU","LUX","LVA")

replace included = 1 if inlist(iso3,"MAR","MDA","MDG","MDV","MEX","MKD","MLI","MLT")
replace included = 1 if inlist(iso3,"MMR","MNG","MOZ","MRT","MUS","MWI","MYS","NAM")
replace included = 1 if inlist(iso3,"NER","NGA","NIC","NLD","NOR","NPL","NZL","OMN")
replace included = 1 if inlist(iso3,"PAK","PAN","PER","PHL","PNG","POL","PRT","PRY")
replace included = 1 if inlist(iso3,"QAT","ROU","RWA","SAU","SDN","SEN","SGP","SLB")
replace included = 1 if inlist(iso3,"SLE","SLV","SOM","SRB","STP","SUR","SVK","SVN")
replace included = 1 if inlist(iso3,"SWE","SWZ","SYR","TCD","TGO","THA","TJK","TKM")
replace included = 1 if inlist(iso3,"TON","TTO","TUN","TUR","TZA","UGA","UKR","URY")
replace included = 1 if inlist(iso3,"USA","UZB","VCT","VEN","VNM","VUT","WSM","YEM")
replace included = 1 if inlist(iso3,"ZAF","COD","ZMB","ZWE","BGD","CHE")
replace iso3 = "ZAR" if iso3 == "COD"
replace iso3 = "ROM" if iso3 == "ROU"
codebook iso3 if included==1
keep if included==1  

g isoo = iso3 
keep if year==2010 
drop year included 

	merge 1:1 isoo using "${dataSTATA}Gstar"
	drop _m
	
	g urban = ln(sp_urb_totl_in_zs/100)
	g gdppc =  ln(gnic)
	g exp  = ln(se_xpd_totl_gd_zs/100)
	g Gstar = ln(G)
	
	g devcat =.
	replace devcat = 1 if lic==1 
	replace devcat = 2 if lmic==1 
	replace devcat = 3 if umic==1 
	replace devcat = 4 if hic ==1 
	
	lab define devcat 1 "Low" 2 "Lower-middle" 3 "Upper-middle" 4 "High-income"
	lab values devcat devcat 
	
	lab var urban "Log urban population (as \% of total population)"
	lab var exp "Log Public expenditure (as \% of GDP)"
	lab var gdppc "Log GDP per capita"
	
*** Validation exercise per se 
local controls exp urban gdppc 
local controlsa exp urban c.exp#devcat
	eststo clear 
	eststo : quiet regress Gstar `controls', vce(rob)
	eststo m1
	estadd local oFE "-", replace

	eststo : quiet reghdfe Gstar `controlsa', abs(devcat) vce(rob)
	eststo m2
	estadd local oFE "\checkmark", replace
	
	esttab m1 m2, label nomtitles se ar2 star(* 0.10 ** 0.05 *** 0.01) drop(1.devcat#c.exp)  
	
file open mytexfile using "$TAB\Table_B2.tex", write replace
file write mytexfile "\begin{tabular}{lcc}" _n
file write mytexfile "\hline \hline" _n
file write mytexfile "&\multicolumn{2}{c}{\$ln G_{i}\$} \\" _n
file write mytexfile "\cmidrule(lr){2-3}" _n
file write mytexfile "& (1) & (2) \\" _n
file close mytexfile 

estout m1 m2  using "$TAB\Table_B2.tex", append style(tex) ///
	   cells(b(star fmt(%9.3f)) se(par fmt(%9.3f))) ///
	   starlevels(\$^*\$ 0.10 \$^{**}\$ 0.05 \$^{***}\$ 0.01) ///
	   drop(1.devcat#c.exp) /// 
	   stats(N r2 oFE, fmt(%9.0fc %9.3fc  %9.0fc) ///
	   labels("\hline Observations" "Adj. R-squared" "Income-group FE")) ///
	   label coll(none) mlab(none)  extracols(4) prehead(\hline) ///
	   postfoot(\hline \hline \end{tabular})  
	   
	   
	
/*------------------------------------------------------------------------------
	Table B.3 : Human capital response to skill biased emigration 
*------------------------------------------------------------------------------*/

***	Data preparation 

	import delimited "$dataMATLAB\WAGES_25.csv", clear case(preserve) 
	g id = _n 
	tempfile Wages 
	save `Wages'

	import delimited "$dataMATLAB\Benchmark_IER_25.csv", clear 
	keep iso nir6_nm nir6_pesc nir6_sc nir7_nm nir7_pesc nir7_sc
	ren (nir6_nm nir6_pesc nir6_sc) (nir6_nmbenc nir6_pescbenc nir6_scbenc)
	ren (nir7_nm nir7_pesc nir7_sc) (nir7_nmbenc nir7_pescbenc nir7_scbenc)
	g id = _n 
	tempfile Benchmark 
	save `Benchmark'
	
	import delimited "$dataMATLAB\cc_channelsNM_b_25.csv", clear case(preserve)
	g id = _n 
	ren (HumCap TecExt DiaExt FisExt MktExt RemEff NatEff) (HumCap_b TecExt_b DiaExt_b FisExt_b MktExt_b RemEff_b NatEff_b)
	tempfile Channels 
	save `Channels'
	
	import delimited "$dataMATLAB\pessimistic_IER_25.csv", clear 
	keep iso nir6_nm nir6_pesc nir6_sc nir7_nm nir7_pesc nir7_sc
	ren (nir6_nm nir6_pesc nir6_sc) (nir6_nmp nir6_pescp nir6_scp)
	ren (nir7_nm nir7_pesc nir7_sc) (nir7_nmp nir7_pescp nir7_scp)	
	g id = _n 
	tempfile pessimistic 
	save `pessimistic'
	
	import delimited "$dataMATLAB\HumanCap_SR_25.csv", clear  case(preserve)
	keep iso HrNMp HrPESCp HrSCp Lambdap LambdaNMp LambdaPESCp LambdaSCp
	tempfile shortrun 
	save `shortrun'

	import delimited "$dataMATLAB\HumanCap_Bench_25.csv", clear  case(preserve)
	merge 1:1 iso using `shortrun'
	drop _m
	merge 1:1 iso using `Benchmark'
	drop _m	
	merge 1:1 iso using `pessimistic'
	drop _m	
	merge 1:1 id using `Channels'
	drop _m
	merge 1:1 id using `Wages'
	drop _m
	
	g llmic = LIC+LMIC
	g umhic = UMIC+HIC
	g hti = Niih/(Niih+Niil)
	g p = mH-mL

	
	local gammaLLMIC_SR  = 1.311
	local gammaLLMIC_LR  = 3.238
	g gammaSR = `gammaLLMIC_SR'*llmic 
	g gammaLR = `gammaLLMIC_LR'*llmic 
	g Hti = (Niih+Mh)/(Niih+Mh+Niil+Ml)
	g HtiNM_SRe = Hn*exp(-gammaSR*p)
	g HtiNM_LRe = Hn*exp(-gammaLR*p)
	
	g dhSRe = hti-HtiNM_SRe
	g dhLRe = hti-HtiNM_LRe
	g dhNMLR = hti-HrNM 
	g dhNMSR = hti-HrNMp

*** Table export 
local vlist Lambda Hn hti p HtiNM_LRe dhLRe HtiNM_SRe dhSRe LambdaNM HrNM dhNMLR HrNMp dhNMSR	
table (iso) var (result),  stat(mean `vlist')  notot nformat(%9.3fc) 
collect export "$TAB\Table_B3.xlsx", sheet("HumanCap") modify 	


/*------------------------------------------------------------------------------
	Table B.4 : Welfare implications for those left behind 
*------------------------------------------------------------------------------*/

***	Data preparation 

	import delimited "$dataMATLAB\WAGES_25.csv", clear case(preserve) 
	g id = _n 
	tempfile Wages 
	save `Wages'

	import delimited "$dataMATLAB\Benchmark_IER_25.csv", clear 
	keep iso nir6_nm nir6_pesc nir6_sc nir7_nm nir7_pesc nir7_sc
	ren (nir6_nm nir6_pesc nir6_sc) (nir6_nmbenc nir6_pescbenc nir6_scbenc)
	ren (nir7_nm nir7_pesc nir7_sc) (nir7_nmbenc nir7_pescbenc nir7_scbenc)
	g id = _n 
	tempfile Benchmark 
	save `Benchmark'
	
	import delimited "$dataMATLAB\cc_channelsNM_b_25.csv", clear case(preserve)
	g id = _n 
	ren (HumCap TecExt DiaExt FisExt MktExt RemEff NatEff) (HumCap_b TecExt_b DiaExt_b FisExt_b MktExt_b RemEff_b NatEff_b)
	tempfile Channels 
	save `Channels'
	
	import delimited "$dataMATLAB\pessimistic_IER_25.csv", clear 
	keep iso nir6_nm nir6_pesc nir6_sc nir7_nm nir7_pesc nir7_sc
	ren (nir6_nm nir6_pesc nir6_sc) (nir6_nmp nir6_pescp nir6_scp)
	ren (nir7_nm nir7_pesc nir7_sc) (nir7_nmp nir7_pescp nir7_scp)	
	g id = _n 
	tempfile pessimistic 
	save `pessimistic'
	
	import delimited "$dataMATLAB\HumanCap_SR_25.csv", clear  case(preserve)
	keep iso HrNMp HrPESCp HrSCp Lambdap LambdaNMp LambdaPESCp LambdaSCp
	tempfile shortrun 
	save `shortrun'

	import delimited "$dataMATLAB\HumanCap_Bench_25.csv", clear  case(preserve)
	merge 1:1 iso using `shortrun'
	drop _m
	merge 1:1 iso using `Benchmark'
	drop _m	
	merge 1:1 iso using `pessimistic'
	drop _m	
	merge 1:1 id using `Channels'
	drop _m
	merge 1:1 id using `Wages'
	drop _m
	
	g llmic = LIC+LMIC
	g umhic = UMIC+HIC
	g hti = Niih/(Niih+Niil)
	g p = mH-mL

	
	local gammaLLMIC_SR  = 1.311
	local gammaLLMIC_LR  = 3.238
	g gammaSR = `gammaLLMIC_SR'*llmic 
	g gammaLR = `gammaLLMIC_LR'*llmic 
	g Hti = (Niih+Mh)/(Niih+Mh+Niil+Ml)
	g HtiNM_SRe = Hn*exp(-gammaSR*p)
	g HtiNM_LRe = Hn*exp(-gammaLR*p)
	
	g dhSRe = hti-HtiNM_SRe
	g dhLRe = hti-HtiNM_LRe
	g dhNMLR = hti-HrNM 
	g dhNMSR = hti-HrNMp

*** Table export 	
local vlist p nir6_nmbenc nir6_nmp HumCap_b TecExt_b DiaExt_b FisExt_b MktExt_b RemEff_b  
table (iso) var (result),  stat(mean `vlist')  notot nformat(%9.3fc) 
collect export "$TAB\Table_B4.xlsx", sheet("Welfare") modify 




/*------------------------------------------------------------------------------
	Table B.4 : Theil index of income inequality 
*------------------------------------------------------------------------------*/

***	Data preparation

	import delimited "$dataMATLAB\THEIL_IER_25.csv", clear  case(preserve)
	g type = _n
	lab define type 1 "Total" 2 "Across" 3 "Within"
	lab values type type 
	order type 
	keep type THEIL_Obs THEIL_NMc THEIL_NM THEIL_NMc_p THEIL_NM_p
	
	rename (THEIL_Obs THEIL_NMc THEIL_NM ///
			THEIL_NMc_p THEIL_NM_p) ///
		   (T_Obs T_NMc T_NM T_NMc_p T_NM_p)
			
foreach var in T_Obs T_NMc T_NM T_NMc_p T_NM_p {
	g dl`var' = (T_Obs/`var'-1)
}

foreach var in T_Obs T_NMc T_NM T_NMc_p T_NM_p {
	g S`var' = `var' if type==1
	egen m`var' = mean(S`var')
	g sh`var' = (`var'/m`var')
}




*** Inequality table per se 
local vlist T_Obs T_NM_p T_NMc_p T_NM T_NMc 	
local dlvlist dlT_Obs dlT_NM_p dlT_NMc_p dlT_NM dlT_NMc	
local shvlist shT_Obs shT_NM_p  shT_NMc_p shT_NM	shT_NMc

* Panel A - Inequality in level 
table (type) var (result),  stat(mean `vlist')  notot nformat(%9.3fc) 
collect export "$TAB\Table_D1.xlsx", sheet("Levels") modify 	

*Panel B - Relative change in inequality in response to selective migration 
table (type) var (result),  stat(mean `dlvlist')  notot nformat(%9.3fc) 
collect export "$TAB\Table_D1.xlsx", sheet("RelDev") modify 

*Panel C - Relative size of within and across components of the Theil 
table (type) var (result),  stat(mean `shvlist')  notot nformat(%9.3fc) 
collect export "$TAB\Table_D1.xlsx", sheet("Shares") modify 







*%%%------------------------------------------------------------------------%%%*
*%%%							B - FIGURES 								%%%*
*%%%------------------------------------------------------------------------%%%*


/*-----------------------------------------------------------------------------*
 Figure B1: Calibration of z by country income group
				Directly plotted from Matlab 
*-----------------------------------------------------------------------------*/


/*-----------------------------------------------------------------------------*
 Figure B2: Effect of selective emigration on human capital accumulation (hi)
				Insights from the empirical model 
*-----------------------------------------------------------------------------*/


***	Data preparation 

	import delimited "$dataMATLAB\HumanCap_SR_25.csv", clear  case(preserve)
	keep iso HrNMp HrPESCp HrSCp Lambdap LambdaNMp LambdaPESCp LambdaSCp
	tempfile shortrun 
	save `shortrun'

	import delimited "$dataMATLAB\HumanCap_Bench_25.csv", clear  case(preserve)
	merge 1:1 iso using `shortrun'
	drop _m 
	g llmic = LIC+LMIC
	g umhic = UMIC+HIC
	g hti = Niih/(Niih+Niil)
	g p = mH-mL
	
	local gammaLLMIC_SR  = 1.311
	local gammaLLMIC_LR  = 3.238
	g gammaSR = `gammaLLMIC_SR'*llmic 
	g gammaLR = `gammaLLMIC_LR'*llmic 
	g Hti = (Niih+Mh)/(Niih+Mh+Niil+Ml)
	g HtiNM_SRe = Hn*exp(-gammaSR*p)
	g HtiNM_LRe = Hn*exp(-gammaLR*p)
	
	g dhSRe = hti-HtiNM_SRe
	g dhLRe = hti-HtiNM_LRe
	g dhNMLR = hti-HrNM 
	g dhNMSR = hti-HrNMp
	

*** Plots 

***	Short-term (SR) effect (∆hi)
	sc dhSRe HtiNM_SRe [aw=N], ms(Oh) mc(navy) mlw(medthick) ///
	   xla(, labs(medium)) yla(-.1(.05).05, labs(large)) ///
	   yli(0, lpattern(dash) lcolor(red)) ///
	   yti("Absolute deviation from closed" "economy", size(large)) ///
	   xti("No-emigration resident's human capital (SR)", size(large))
	graph export "$FIG\Figure_B2a.png", replace
	

***	Long-term (LR) effect (∆hi)	
	sc dhLRe HtiNM_LRe [aw=N], ms(Oh) mc(navy) mlw(medthick) ///
	   yla(-.1(.05).1, labs(medium)) xla(, labs(large))	///
	   yli(0, lpattern(dash) lcolor(red)) ///
	   yti("Absolute deviation from closed" "economy", size(large)) ///
	   xti("No-emigration resident's human capital (LR)", size(large))
	graph export "$FIG\Figure_B2b.png", replace	   

	   	   

***	Generalized vs. Empirical (∆hi in SR)   
	twoway sc dhNMSR dhSRe [aw=N], ms(Oh) mc(navy) mlw(medthick) ///
			  yla(-.1(.05).1, labs(medium)) ||  ///
			  function y=x, range(-.1 .1) xlab(, labs(large)) /// 
			  xti("Econometric approach (SR)", size(large)) ///
			  yti("Micro-founded approach (SR)", size(large)) leg(off)
	graph export "$FIG\Figure_B2c.png", replace
	

***	Generalized vs. Empirical (∆hi in LR)
	twoway sc dhNMLR dhLRe [aw=N], ms(Oh) mc(navy) mlw(medthick) ///
			  xla(, labs(large)) yla(-.1(.05).1, labs(large)) || ///
			  function y=x, range(-.1 .1) ///
			  xti("Econometric approach (LR)", size(large)) ///
			  yti("Micro-founded approach (LR)", size(large)) leg(off)
	graph export "$FIG\Figure_B2d.png", replace			 


***	Density of ∆hi (SR)			 
	twoway kdensity dhSRe, bw(.002)  lpattern(solid) lcolor(black) lwidth(thick) || ///
		   kdensity dhNMSR, bw(.0025) lpattern(dash) lcolor(black) lwidth(thick)  ///
			 xti("Absolute deviation from closed economy (SR)", size(large)) ///
			 yti("Density", size(large)) ///
			 xlab(-.1(.05).1, labs(large)) ///
			 xline(0, lcolor(red) lpattern(dash) lwidth(.3pt)) ///
			 note("") title("") ylab(0(20)100, labs(large)) legend(pos(6) cols(3) ///
			 label(1 "Econometric") ///
			 label(2 "Micro-founded") size(large))
	graph export "$FIG\Figure_B2e.png", replace
	

***	Density of ∆hi (LR)	
	twoway kdensity dhLRe, lpattern(solid) lcolor(black) lwidth(thick) || ///
		   kdensity dhNMLR, lpattern(dash) lcolor(black) lwidth(thick)  ///
			 xti("Absolute deviation from closed economy (LR)", size(large)) ///
			 yti("Density", size(large)) ///
			 xlab(-.1(.05).1, labs(large)) yla(,labs(large)) ///
			 xline(0, lcolor(red) lpattern(dash) lwidth(.3pt)) ///
			 note("") title("") legend(pos(6) cols(3) ///
			 label(1 "Econometric") ///
			 label(2 "Micro-founded") size(large))
	graph export "$FIG\Figure_B2f.png", replace	
	
	
	
/*-----------------------------------------------------------------------------*
 Figure B3: Effect of selective emigration on human capital accumulation (hi)
 Generalized approach with conservative parameter set 
					Insights from the empirical model 
*-----------------------------------------------------------------------------*/
/*
	Data preparation 
*/
	import delimited "$dataMATLAB\HumanCap_SR_25.csv", clear  case(preserve)
	keep iso HrNMp HrPESCp HrSCp Lambdap LambdaNMp LambdaPESCp LambdaSCp
	tempfile shortrun 
	save `shortrun'
	
	import delimited "$dataMATLAB\HumanCap_Bench_25.csv", clear  case(preserve)
	merge 1:1 iso using `shortrun'
	drop _m 
	g llmic = LIC+LMIC
	g umhic = UMIC+HIC
	g hti = Niih/(Niih+Niil)
	g p = mH-mL
	
	local gammaLLMIC_SR  = 1.311
	local gammaLLMIC_LR  = 3.238
	g gammaSR = `gammaLLMIC_SR'*llmic 
	g gammaLR = `gammaLLMIC_LR'*llmic 
	 
	g HtiNM_SRe = Hn*exp(-gammaSR*p)
	g HtiNM_LRe = Hn*exp(-gammaLR*p)
	
	g dhSRe = hti-HtiNM_SRe
	g dhLRe = hti-HtiNM_LRe
	g dhNMLR = hti-HrNM 
	g dhNMSR = hti-HrNMp
	
/*
	Short-term (SR) effect (∆hi)
*/
	scatter dhNMSR HtiNM_SR [aw=N], ms(Oh) mc(navy) mlw(medthick) ///
	        yla(-.04(.02).08,labs(large)) xla(,labs(large)) ///
			yli(0, lpattern(dash) lcolor(red)) ///
			yti("Absolute deviation from closed economy", size(large)) ///
			xti("No-emigration resident's human capital (SR)", size(large))   
	graph export "$FIG\Figure_B3a.png", replace
	
/*
	Observed vs. no-migration (SR)
*/
	twoway scatter hti HrNMp [aw=N], ms(Oh) mc(navy) yla(0(.1).5, labs(large)) mlw(medthick) || ///
		   function y=x, range(0 .5) xlab(, labs(large)) ///
		   xti("Econometric approach (SR)", size(large)) ///
		   yti("Micro-founded approach (SR)", size(large)) leg(off)
	graph export "$FIG\Figure_B3b.png", replace	
	   	   

		   
/*-----------------------------------------------------------------------------*
 Figure C1: Effect of selective emigration on disposable income (yi)
		Results obtained under the conservative scenario 
(short-term human capital responses and conservative elasticity values) 
*-----------------------------------------------------------------------------*/	

/*
	Data preparation 
*/
	import delimited "$dataMATLAB\WAGES_25.csv", clear case(preserve) 
	g id = _n 
	tempfile Wages 
	save `Wages'

	import delimited "$dataMATLAB\pessimistic_IER_25.csv", clear 
	g id = _n 
	tempfile pessimistic 
	save `pessimistic'
	
	import delimited "$dataMATLAB\cc_channelsNM_p_25.csv", clear case(preserve)
	g id = _n 
	tempfile Channels 
	save `Channels'
	
	
	use `Wages', clear 
	merge 1:1 id using `pessimistic'
	drop _m	
	merge 1:1 id using `Channels'
	drop _m 
	
/*
	Figure preparation  
*/

/*
	Figure C.1a
*/
	kdens NIR6_NMp, gen(d1 x1) ci(ci1_1 ci1_2) ///
		  xla(-.2(.1).4, labs(large)) yla(, labs(large)) lcolor(black) ///
		  lwidth(thick) xli(0, lpattern(dash) lcolor(red)) ///
		  xti("Net per worker disposable income response", size(large)) ///
		  yti("Density", size(large)) ///
		  legend(position(6) cols(2) label(1 "95% CI") ///
		  label(2 "Kernel density") size(large))
	graph export "$FIG\Figure_C1a.png", replace 	


/*
	Figure C1.b
*/	
twoway kdensity HumCap_p, bwidth(.004) lpattern(solid) lcolor(black) lwidth(medthick) || ///
	   kdensity TecExt_p, bwidth(.002) lpattern(dash) lcolor(black) lwidth(medthick)  || ///
	   kdensity DiaExt_p, bwidth(.005) lpattern(solid) lcolor(blue) lwidth(medthick)  || ///
	   kdensity FisExt_p, bwidth(.010) lpattern(vshortdash) lcolor(black) lwidth(medthick)  || ///
	   kdensity MktExt_p, bwidth(.005) lpattern(solid) lcolor(black)  lwidth(small)  || ///
	   kdensity RemEff_p, bwidth(.010) lpattern(dash) lcolor(blue) lwidth(medthick)  ///
			 xti("Percentage change in per worker disposable income", size(large)) ///
			 yti("Density", size(large)) ///
			 xlab(-.1(.1).4, labs(large)) ///
			 xline(0, lcolor(red) lpattern(dash) lwidth(.3pt)) ///
			 note("") title("") ylab(0(30)90, labs(large)) legend( ring(0) pos(2) cols(1) ///
			 label(1 "Human capital") ///
			 label(2 "Technological ext.") ///
			 label(3 "Diaspora ext.") ///
			 label(4 "Fiscal ext.") ///
			 label(5 "Market size ext.") ///
			 label(6 "Remittances") size(large))
	graph export "$FIG\Figure_C1b.png", replace
 
 /*
	Figure C1.c
*/	


	g ywd_nmp = 	ywd/(1+NIR6_NMp)
	g lywd_nmp = ln(ywd_nmp)

	twoway lfitci NIR6_NMp lywd_nmp, legend(off) lcolor(red) alwidth(none) || ///
	       scatter NIR6_NMp lywd_nmp [aw=LNM], ms(Oh) mcolor(navy) mlw(medthick) ///
			   mlabpos(0) mlabcolor(black) mlabsize(vsmall) ///
			   xti("Disposable income per worker (no-emigration)", size(large)) ///
			   yti("Net per worker disposable" "income response", size(large)) ///
			   xla(,labs(medium)) yla(-.2(.1).3,labs(large)) yli(0, lpattern(shortdash))
	graph export "$FIG\Figure_C1c.png", replace 
	
	
 /*
	Figure C1.d
*/		

	twoway (kdensity NIR6_NMp, lcolor(black) lpattern(solid) lwidth(thick)) ///
		   (kdensity NIR7_NMp, lcolor(black) lpattern(dash) lwidth(thick) ///
		    legend(ring(0) position(2) cols(1) label(1 "Per worker") label(2 "Per natural") size(large)) ///
			xla(-.2(.2)1, labs(large)) yla(0(5)15, labs(large)) xti("Net per worker disposable income response", size(large)) ///
			yti("Density", size(large)) xli(0, lpattern(shortdash) lcolor(red)))
	graph export "$FIG\Figure_C1d.png", replace 
	

	
/*-----------------------------------------------------------------------------*
 Figure C2: Average disposable income responses to selective emigration
				(conservative variant)
*-----------------------------------------------------------------------------*/	

				   
	import delimited "$dataMATLAB\AVG_NIR_25.csv", clear case(preserve)  
	g devcat = _n 
	
	lab define devcat 1 "World" 2 "LIC" 3 "LMIC" 4 "UMIC" 5 "HIC"
	lab values devcat devcat 	
	
	foreach var in AVG_NIR6_NMp AVG_NIR7_NMp {
		replace `var' = 100*`var'
	}
	
	graph bar AVG_NIR6_NMp AVG_NIR7_NMp, over(devcat, label(labs(large))) ///
		  yla(0(5)30,labs(large)) ///
		  bar(1, fcolor(navy) fintensity(50) color(none)) ///
		  bar(2, fcolor(maroon) fintensity(50) color(none)) ///
		  legend(position(6) cols(2) label(1 "Per worker") ///
		  label(2 "Per natural") size(large)) ///
		  yti("Net disposable income response", size(large))  
	graph export "$FIG\Figure_C2.png", replace 	
	

/*-----------------------------------------------------------------------------*
 Figure C3: Impact of selection per se (NS)
*-----------------------------------------------------------------------------*/

/*
	Data preparation 
*/
	preserve 
	import delimited "$dataMATLAB\WAGES_25.csv", clear case(preserve) 
	keep iso mH mL mi L Ll Lh LSC LlSC LhSC LPESC LlPESC LhPESC LNS LlNS LhNS LNM LlNM LhNM wH wL w wLd wHdd wLdd wHGEd wLGEd
	g id = _n  
	tempfile metrics_A 
	save `metrics_A'
	
	import delimited "$dataMATLAB\HumanCap_Bench_25.csv", clear case(preserve)
	keep iso Lambda LambdaNM LambdaPESC LambdaSC LIC LMIC UMIC HIC
	g id = _n  
	tempfile metrics_B 
	save `metrics_B'
	
	use `metrics_A', clear 
	merge 1:1 iso using `metrics_B'
	drop _m 
	
	tempfile  metrics 
	save `metrics'
	restore 
	
	import delimited "$dataMATLAB\Benchmark_IER_25.csv", clear  case(preserve)
	g id = _n  
	tempfile Benchmark 
	save `Benchmark'

	import delimited "$dataMATLAB\Pessimistic_IER_25.csv", clear  case(preserve)
	g id = _n  
	tempfile pessimistic 
	save `pessimistic'

 
	import delimited "$dataMATLAB\WAGES_25.csv", clear case(preserve) 
	g id = _n 
	tempfile Wages 
	save `Wages'

	use `metrics', clear 
	merge 1:1 id using `Benchmark'
	drop _m 

	merge 1:1 id using `pessimistic'
	drop _m 
	merge 1:1 id using `Wages'
	drop _m 

/*
	Figure preparation  
*/
		  
	kdens NIR6_PESC4,  gen(d1 x1) ci(ci1_1 ci1_2) ///
		  xla(-.1(.05).2, labs(large)) yla(0(10)40,labs(large)) lcolor(black) lwidth(thick) ///
		  xli(0, lpattern(dash) lcolor(red)) ///
		  xti("Net per worker disposable income response - Benchmark", size(large)) ///
		  yti("Density", size(large)) legend(position(6) cols(2) ///
		  label(1 "95% CI") label(2 "Kernel density") size(large))
	graph export "$FIG\Figure_C3a.png", replace 
	

	g ywd_NS = 	ywd/(1+NIR6_PESC4)
	g lywd_NS = ln(ywd_NS)

 			   
	twoway lfitci NIR6_PESC4 lywd_NS, legend(off) lcolor(red) alwidth(none) || ///
	       scatter NIR6_PESC4 lywd_NS [aw=LNS], ms(Oh) mcolor(navy) mlw(medthick) ///
			   mlabpos(0) mlabcolor(black) mlabsize(large) ///
			   xti("Disposable income per worker (no-selection)", size(large)) ///
			   yti("Net per worker disposable" "income response", size(large)) ///
			   xla(,labs(large)) yla(-.1(.1).2,labs(large)) yli(0, lpattern(shortdash))
	graph export "$FIG\Figure_C3b.png", replace 	   
	   

			   
	import delimited "$dataMATLAB\AVG_NIR_25.csv", clear case(preserve)  
	g devcat = _n 
	
	lab define devcat 1 "World" 2 "LIC" 3 "LMIC" 4 "UMIC" 5 "HIC"
	lab values devcat devcat 	
	
	ren (AVG_NIR6_NM AVG_NIR7_NM AVG_NIR6_NS AVG_NIR7_NS AVG_NIR6_NMp AVG_NIR7_NMp AVG_NIR6_NSp AVG_NIR7_NSp) ///
		(AVG_NIR_NM6 AVG_NIR_NM7 AVG_NIR_NS6 AVG_NIR_NS7 AVG_NIR_NMp6 AVG_NIR_NMp7 AVG_NIR_NSp6 AVG_NIR_NSp7)
		
	reshape long AVG_NIR_NM AVG_NIR_NS AVG_NIR_NMp AVG_NIR_NSp, i(devcat) j(nature)
	lab define nature 6 "Per worker" 7 "Per natural"
	lab values nature nature 
	
	
	graph bar AVG_NIR_NM AVG_NIR_NS AVG_NIR_NMp AVG_NIR_NSp if devcat==1 & nature==6, yla(0(.01).04,labs(large)) ///
		  bar(1, fcolor(black) fintensity(80) color(none)) ///
		  bar(2, fcolor(black) fintensity(40) color(none)) ///
		  bar(3, fcolor(blue) fintensity(80) color(none)) ///
		  bar(4, fcolor(blue) fintensity(40) color(none)) ///
		  legend(position(6) cols(4) label(1 "NM (Bench.)") ///
									 label(2 "NS (Bench.)") ///
									 label(3 "NM (Cons.)") ///
									 label(4 "NS (Cons.)") size(large)) ///
									 yti("Net disposable income response", size(large))
	graph export "$FIG\Figure_C3c.png", replace 		
	
	graph bar AVG_NIR_NM AVG_NIR_NS AVG_NIR_NMp AVG_NIR_NSp, over(devcat, label(labs(large) angle(60))) by(nature, note("")) yla(0(.1).3, labs(large)) ///
		  bar(1, fcolor(black) fintensity(80) color(none)) ///
		  bar(2, fcolor(black) fintensity(40) color(none)) ///
		  bar(3, fcolor(blue) fintensity(80) color(none)) ///
		  bar(4, fcolor(blue) fintensity(40) color(none)) ///
		  legend(position(6) cols(1) label(1 "NM (Bench.)") ///
									 label(2 "NS (Bench.)") ///
									 label(3 "NM (Cons.)") ///
									 label(4 "NS (Cons.)") size(large)) ///
									 yti("Net disposable income response", size(large))
	graph export "$FIG\Figure_C3d.png", replace 	
	
	
		
/*-----------------------------------------------------------------------------*
 Figure C4: Selection and world distribution of income
*-----------------------------------------------------------------------------*/
/*
	Data preparation 
*/
	
	import delimited "$dataMATLAB\WAGES_25.csv", clear case(preserve) 
	keep iso Lh LhSC LhPESC LhNS wHSC4d wHGE4d wHSCd wHGEd wHdd wH mH
		ren (Lh LhSC LhPESC LhNS wHSC4d wHGE4d wHSCd wHGEd wHdd wH mH) ///
			(L  LSC  LPESC  LSC4 wSC4d wGE4d wSCd wGEd wdd w m)
	
	g skill = 1 
	tempfile HighSkill 
	save `HighSkill'
	
	import delimited "$dataMATLAB\WAGES_25.csv", clear case(preserve) 
	keep iso Ll LlSC LlPESC LlNS wLSC4d wLGE4d wLSCd wLGEd wLdd wL mL
		ren (Ll LlSC LlPESC LlNS wLSC4d wLGE4d wLSCd wLGEd wLdd wL mL) ///
			(L  LSC  LPESC  LSC4 wSC4d wGE4d wSCd wGEd wdd w m)
	
	g skill = 2
	append using `HighSkill'	
	lab define skill 1 "High Skilled" 2 "Low Skilled"
	lab values skill skill 
	
	foreach var in w wdd wGEd wSCd wGE4d wSC4d {
		g ln`var' = ln(`var')
	}

	g LSCa = int(LSC)
	g LSC4a = int(LSC4)
	g La = int(L)


/*
	Figure preparation  
*/			
	
	twoway (kdensity lnwdd [fw=La], bw(.17) lcolor(black) lpattern(solid) lwidth(large)) ///
		   (kdensity lnwSC4d [fw=La], bw(.17) lcolor(black) lpattern(dash) lwidth(large)) ///
		   (kdensity lnwSC4d [fw=LSC4a], bw(.17) lcolor(blue) lpattern(solid) lwidth(large) /// 
		    legend(position(6) cols(3) label(1 "Observed") ///
									   label(2 "NS(Initial labor force)") ///
									   label(3 "NS(New labor force)") size(large)) ///
			xla(7(1)12, labs(large)) xla(, labs(large)) ///
			xti("Disposable income per worker", size(large)) ///
			yti("Density", size(large)) xli(7.603 9.426, lpattern(shortdash) lcolor(red)) ///
			text(0.3 7.47 "International poverty line ($2,004)", place(e) orientation(vertical) size(.3cm)) ///
			text(0.25 9.3 "Median disposable income ($12,404)", place(e) orientation(vertical) size(.3cm)))
	graph export "$FIG\Figure_C4.png", replace 	

	
	
/*-----------------------------------------------------------------------------*
 Figure C5: Sensitivity to mu and sigma
*-----------------------------------------------------------------------------*/
/*
	Data preparation 
*/	

	import delimited "$dataMATLAB\Benchmark_IER_25.csv", clear 
	keep iso nir6_nm nir6_pesc nir6_sc nir7_nm nir7_pesc nir7_sc
	ren (nir6_nm nir6_pesc nir6_sc) (nir6_nmbenc nir6_pescbenc nir6_scbenc)
	ren (nir7_nm nir7_pesc nir7_sc) (nir7_nmbenc nir7_pescbenc nir7_scbenc)
	g id = _n 
	tempfile Benchmark 
	save `Benchmark'
	
	import delimited "$dataMATLAB\Benchmarksigma15_IER_25.csv", clear 
	keep iso nir6_nm nir6_pesc nir6_sc nir7_nm nir7_pesc nir7_sc
	ren (nir6_nm nir6_pesc nir6_sc) (nir6_nmsigma15 nir6_pescsigma15 nir6_scsigma15)
	ren (nir7_nm nir7_pesc nir7_sc) (nir7_nmsigma15 nir7_pescsigma15 nir7_scsigma15)
	g id = _n 
	tempfile sigma15 
	save `sigma15'	
	
	
	import delimited "$dataMATLAB\Benchmarkmu06_IER_25.csv", clear 
	keep iso nir6_nm nir6_pesc nir6_sc nir7_nm nir7_pesc nir7_sc
	ren (nir6_nm nir6_pesc nir6_sc) (nir6_nmmu06 nir6_pescmu06 nir6_scmu06)
	ren (nir7_nm nir7_pesc nir7_sc) (nir7_nmmu06 nir7_pescmu06 nir7_scmu06)
	g id = _n 
	tempfile mu06 
	save `mu06'		
	
	import delimited "$dataMATLAB\Benchmarkmu142_IER_25.csv", clear 
	keep iso nir6_nm nir6_pesc nir6_sc nir7_nm nir7_pesc nir7_sc
	ren (nir6_nm nir6_pesc nir6_sc) (nir6_nmmu142 nir6_pescmu142 nir6_scmu142)
	ren (nir7_nm nir7_pesc nir7_sc) (nir7_nmmu142 nir7_pescmu142 nir7_scmu142)
	g id = _n 
	tempfile mu142 
	save `mu142'		
	

	use `Benchmark', clear 
	merge 1:1 id using `sigma15'
	drop _m 
	merge 1:1 id using `mu06'
	drop _m 	
	merge 1:1 id using `mu142'
	drop _m 
	

/*
	Figure preparation  
*/			


twoway kdensity nir6_pescbenc, lpattern(solid) lcolor(black) lwidth(medthick) || ///
	   kdensity nir6_pescmu142, lpattern(dash) lcolor(blue) lwidth(medthick)  ||  ///
	   kdensity nir6_pescmu06, lpattern(dash) lcolor(maroon) lwidth(medthick) || ///
	   kdensity nir6_pescsigma15, lpattern(dash) lcolor(navy) lwidth(medthick) ///
			 xti("Net per worker disposable income response", size(large)) ///
			 yti("Density", size(large)) ///
			 xlab(-.3(.1).5, labs(large)) ///
			 xline(0, lcolor(red) lpattern(dash) lwidth(.3pt)) ///
			 note("") title("") ylab(0(5)15, labs(large)) legend(pos(6) cols(2) ///
			 label(1 "Benchmark (1/{&mu} = 0.7, {&sigma} = 2)") ///
			 label(2 "1/{&mu} = 1.42") ///
			 label(3 "1/{&mu} = 0.6") ///
			 label(4 "{&sigma} = 1.5") size(large))
	graph export "$FIG\Figure_C5.png", replace		
	
	
	
/*-----------------------------------------------------------------------------*
 Figure D1: Emigration and the world income distribution: Conservative variant
*-----------------------------------------------------------------------------*/
/*
	Data preparation 
*/	
	import delimited "$dataMATLAB\WAGESp_25.csv", clear case(preserve) 
	keep iso wHdd wHSCd LhSC Lh
	rename (wHdd wHSCd LhSC Lh) (wdd wSCd LSC L)
	g skill = 1 
	tempfile HighSkill 
	save `HighSkill'
	
	import delimited "$dataMATLAB\WAGES_25.csv", clear case(preserve) 
	keep iso wLdd wLSCd LlSC Ll
	rename (wLdd wLSCd LlSC Ll) (wdd wSCd LSC L)
	
	g skill = 2
	append using `HighSkill'	
	lab define skill 1 "High Skilled" 2 "Low Skilled"
	lab values skill skill 
	
	foreach var in wdd wSCd LSC L {
		g ln`var' = ln(`var')
	}

	g LSCa = int(LSC)
	g La = int(L)
	

/*
	Figure preparation  
*/			


	twoway (kdensity lnwdd [fw=La], bw(.17) lcolor(black) lpattern(solid) lwidth(large)) ///
		   (kdensity lnwSCd [fw=La], bw(.17) lcolor(black) lpattern(dash) lwidth(large)) ///
		   (kdensity lnwSCd [fw=LSCa], bw(.17) lcolor(blue) lpattern(solid) lwidth(large) /// 
		    legend(position(6) cols(3) label(1 "Observed") ///
									   label(2 "NM(Initial labor force)") ///
									   label(3 "NM(New labor force)") size(large)) ///
			xla(7(1)12, labs(large)) yla(0(.2).8, labs(large)) ///
			xti("Disposable income per worker", size(large)) ///
			yti("Density", size(large)) xli(7.603 9.426, lpattern(shortdash) lcolor(red)) ///
			text(0.3 7.47 "International poverty line ($2,004)", place(e) orientation(vertical) size(.3cm)) ///
			text(0.25 9.3 "Median disposable income ($12,404)", place(e) orientation(vertical) size(.3cm)))
	graph export "$FIG\Figure_D1.png", replace 		
	

	
	
	
/*-----------------------------------------------------------------------------*
 Figure D2: Income per worker relative to the US
*-----------------------------------------------------------------------------*/	

/*
	Data preparation 
*/
	import delimited "$dataMATLAB\WAGES_25.csv", clear case(preserve) 
	g id = _n 
	tempfile Wages 
	save `Wages'

	import delimited "$dataMATLAB\Benchmark_IER_25.csv", clear 
	keep iso nir6_nm nir6_pesc nir6_sc nir7_nm nir7_pesc nir7_sc
	ren (nir6_nm nir6_pesc nir6_sc) (nir6_nmbenc nir6_pescbenc nir6_scbenc)
	ren (nir7_nm nir7_pesc nir7_sc) (nir7_nmbenc nir7_pescbenc nir7_scbenc)
	g id = _n 
	tempfile Benchmark 
	save `Benchmark'

	import delimited "$dataMATLAB\Pessimistic_IER_25.csv", clear 
	keep iso nir6_nm nir6_pesc nir6_sc nir7_nm nir7_pesc nir7_sc
	ren (nir6_nm nir6_pesc nir6_sc) (nir6_nmp nir6_pescp nir6_scp)
	ren (nir7_nm nir7_pesc nir7_sc) (nir7_nmp nir7_pescp nir7_scp)
	g id = _n 
	tempfile Pessimistic 
	save `Pessimistic'
	
	
	import delimited "$dataMATLAB\cc_channelsNM_b_25.csv", clear case(preserve)
	g id = _n 
	ren (HumCap TecExt DiaExt FisExt MktExt RemEff NatEff) (HumCap_b TecExt_b DiaExt_b FisExt_b MktExt_b RemEff_b NatEff_b)
	tempfile Channels 
	save `Channels'
	
	import delimited "$dataMATLAB\cc_channelsNM_p_25.csv", clear case(preserve)
	g id = _n 
	tempfile Channelsp 
	save `Channelsp'
	
	
	
	use `Wages', clear 
	merge 1:1 id using `Benchmark'
	drop _m	
	merge 1:1 id using `Pessimistic'
	drop _m
	merge 1:1 id using `Channels'
	drop _m 
	merge 1:1 id using `Channelsp'
	drop _m 	
	
	
	
	g Hn_a = Nh/N 
	g Hn_aUS = Hn_a if iso=="USA" 
	egen Hn_aUSx = mean(Hn_aUS)
	drop Hn_aUS
	ren Hn_aUSx Hn_aUS
	local sigma = 2 
	local rho = (`sigma'-1)/`sigma'
	
	g QiT_US = (thetH*Hn_aUS^`rho' + (1-thetH)*(1-Hn_aUS)^`rho')^(1/`rho')
	g tau = (TFP*LiT - ywd)/(TFP*LiT)
	g ydUS = TFP*QiT_US*(1-tau)
	
	
/* Counterfactual disposable income per worker */	
	g ypcd = ywd+rpc
	g ypcd_CFb = ywd/(1+NIR6_NMb) // from benchmark
	g ypcd_CFp = ywd/(1+NIR6_NMp) // from conservative variant  	
	g ypcd_CFUS = ydUS
	
	g ypcd_BGb = ywd/(1+HumCap_b) // Benchmark Brain gain mechanism only 
	g ypcd_BGp = ywd/(1+HumCap_p) // Conservative Brain gain mechanism only 	
	
	foreach var in ypcd ypcd_CFb ypcd_CFp ypcd_BGb ypcd_BGp ypcd_CFUS {
		g `var'US = `var' if iso=="USA"
		replace `var'US = 0 if missing(`var'US)
		egen `var'USx = max(`var'US)
		replace `var'US = `var'USx 
		drop `var'USx 
	}

	
	g gap_to_US = (ypcd/ypcdUS)  
	g gap_to_US_CFb = (ypcd_CFb/ypcd_CFbUS)  
	g gap_to_US_CFp = (ypcd_CFp/ypcd_CFpUS)  
	g gap_to_US_BGb = (ypcd_BGb/ypcd_BGbUS)  
	g gap_to_US_BGp = (ypcd_BGp/ypcd_BGpUS) 
	g gap_to_US_Q =  (ypcd_CFb/ypcd_CFUS) 
	gsort ypcd  
	g rank = _n 
	gsort ypcd_CFb 
	g rank2 = _n 

	
/*
	Figure preparation  
*/

*sc rank ypcd_CFb
g lypcd_CFb = ln(ypcd_CFb)
g lypcd = ln(ypcd)
g lypcd_CFb_inv = -lypcd_CFb
g lypcd_inv = -lypcd
g lypcd_BGb = ln(ypcd_BGb)
g lypcd_BGb_inv = -lypcd_BGb

g llmic =lic+lmic 

drop rank2 
sum ypcd, de 
local med = r(p50)
g line1 = (ypcd<=`med')



preserve
sum ypcd, de 
local med = r(p50)
keep if ypcd<=`med'
	gsort ypcd_CFb 
	g rank2 = _n 
	labmask rank2, val(iso)	
twoway line gap_to_US rank2, lpattern(longdash) lcolor(blue) lwidth(thick) || ///
	   line gap_to_US_BGb rank2, lpattern(dash) lcolor(maroon) lwidth(thick) || ///
	   line gap_to_US_CFb rank2, lpattern(solid) lcolor(black) lwidth(thick) ///
	   yli(1, lpattern(dash) lcolor(red))  name(ConvSpeed, replace) ///
	   yti("Disposable income per worker" "Relative to the US (US = 1)", size(large)) ///
	   xti("No migration disposable income per worker (USD)", size(large)) ///
	   xlab(1(1)87, valuelabel labsize(vsmall) angle(90) nogrid) ///
	   yla(0(.05).25, labsize(large)) legend(ring(0) position(11) cols(1) ///
	   label(1 "Observed (data)") ///
	   label(2 "No migration with brain gain") ///
	   label(3 "No migration with all channels" )size(medium))
	graph export "$FIG2\Figure_D2a.png", replace	   
restore 



	   
preserve 
sum ypcd, de 
local med = r(p50)
keep if ypcd>`med'
	gsort ypcd_CFb 
	g rank2 = _n 
	labmask rank2, val(iso)	
twoway line gap_to_US rank2, lpattern(longdash) lcolor(blue) lwidth(thick) || ///
	   line gap_to_US_BGb rank2, lpattern(dash) lcolor(maroon) lwidth(thick) || ///
	   line gap_to_US_CFb rank2, lpattern(solid) lcolor(black) lwidth(thick) ///
	   yli(1, lpattern(dash) lcolor(red))  name(ConvSpeed, replace) ///
	   yti("Disposable income per worker" "Relative to the US (US = 1)", size(large)) ///
	   xti("No migration disposable income per worker (USD)", size(large)) ///
	   xlab(1(1)87, labsize(vsmall) angle(90) valuelabel nogrid) ///
	   yla(0(.2)1.2, labsize(large)) legend(ring(0) position(9) cols(1) ///
	   label(1 "Observed (data)") ///
	   label(2 "No migration with brain gain") ///
	   label(3 "No migration with all channels" )size(medium))
	graph export "$FIG2\Figure_D2b.png", replace	  
restore 	  

lab define line1 1 "< Median income" 0 "> Median income" 
lab values line1 line1 


g devcat = lic 
replace devcat = 2 if lmic==1
replace devcat = 3 if umic==1
replace devcat = 4 if hic==1 

lab define devcat 1 "LIC" 2 "LMIC" 3 "UMIC" 4 "HIC" 
lab values devcat devcat 


graph bar gap_to_US gap_to_US_BGb gap_to_US_CFb, over(devcat, label(labsize(large))) ///
		  bar(1, color(navy*0.7)) bar(2, color(maroon*0.7)) bar(3, color(black*0.4)) ///
		  legend(ring(0) position(11) cols(1) ///
		  label(1 "Observed (data)") ///
	      label(2 "No migration with brain gain") ///
	      label(3 "No migration with all channels" )size(medium)) ///
		  bargap(10) ylabel(0(.1).5, /*gstyle(dot)*/ labs(large)) ///
		  blabel(bar,format(%9.3f) color(navy) pos(outside) size(medsmall)) ///
		  yti("Income per worker relative to the US (US=1)", size(large))
	graph export "$FIG\Figure_D2c.png", replace	
	
/* relative variations */ // relative between the observed income gap and the counterfactual one 
g dgapBG = .
replace dgapBG = -5.7 if lic==1 
replace dgapBG = -9.9 if lmic==1 
replace dgapBG = -2.4 if umic==1 
replace dgapBG = -1.7 if hic==1 

g dgapCF = .
replace dgapCF = -11.4 if lic==1 
replace dgapCF = -15.4 if lmic==1 
replace dgapCF = -5.3 if umic==1 
replace dgapCF = -3.6 if hic==1 

graph bar dgapBG dgapCF, over(devcat, label(labsize(large))) ///
		  bar(1, color(navy*0.7)) bar(2, color(maroon*0.7)) ///
		  legend(ring(0) position(11) cols(1) ///
	      label(1 "No migration with brain gain") ///
	      label(2 "No migration with all channels" )size(medium)) ///
		  bargap(10) ylabel(-20(5)5, /*gstyle(dot)*/ labs(large)) ///
		  blabel(bar,format(%9.1f) color(navy) pos(outside) size(medsmall)) ///
		  yti("Relative change in per worker income gap" "with the US (%)", size(large))
	graph export "$FIG\Figure_D2d.png", replace	
