

/*
	Table 2: Selective emigration and education incentives: short- and long-run effects
*/

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

		
	
eststo clear
local tex_file_name "Table_2.tex"
local latex_header ///
	" & \multicolumn{4}{c}{Ordinary Least Squares} & \multicolumn{4}{c}{Instrumental Variables} \\ \cmidrule(lr){2-5} \cmidrule(lr){6-9}  & \multicolumn{2}{c}{Developing countries} & \multicolumn{2}{c}{Full sample}  & \multicolumn{2}{c}{Developing countries} & \multicolumn{2}{c}{Full sample} \\ & (1) & (2) & (3) & (4) & (5) & (6) & (7) & (8) \\" 


	/* Panel A. Short run estimates */
*** OLS 
	local panel_text "Panel A- Short-term estimates"
eststo m1: quiet reghdfe $yvar  $endog $controls i.income_gp if sample2==1, abs(year) cl(id)
			estadd local dFE "\checkmark", replace
			estadd local tFE "\checkmark", replace
			estadd local oFE "-", replace
eststo m2: quiet reghdfe $yvar  $endog $controls  i.income_gp  admit if sample2==1, abs(year) cl(id)
			estadd local dFE "\checkmark", replace
			estadd local tFE "\checkmark", replace
			estadd local oFE "\checkmark", replace	
eststo m3: quiet reghdfe $yvar  $endog $controls  i.income_gp  if sample1==1, abs(year) cl(id)
			estadd local dFE "\checkmark", replace
			estadd local tFE "\checkmark", replace
			estadd local oFE "-", replace	
eststo m4: quiet reghdfe $yvar  $endog $controls  i.income_gp  admit if sample1==1, abs(year) cl(id)
			estadd local dFE "\checkmark", replace
			estadd local tFE "\checkmark", replace
			estadd local oFE "\checkmark", replace		
*** IV 
eststo m5: quiet ivreghdfe $yvar  ($endog = $instr) $controls  i.income_gp  if sample2==1, abs(year) cl(id) first
			estadd local dFE "\checkmark", replace
			estadd local tFE "\checkmark", replace
			estadd local oFE "-", replace
eststo m6: quiet ivreghdfe $yvar  ($endog = $instr) $controls  i.income_gp  admit if sample2==1, abs(year) cl(id) first
			estadd local dFE "\checkmark", replace
			estadd local tFE "\checkmark", replace
			estadd local oFE "\checkmark", replace
eststo m7: quiet ivreghdfe $yvar  ($endog = $instr) $controls  i.income_gp  if sample1==1, abs(year) cl(id) first
			estadd local dFE "\checkmark", replace
			estadd local tFE "\checkmark", replace
			estadd local oFE "-", replace
eststo m8: quiet ivreghdfe $yvar  ($endog = $instr) $controls  i.income_gp  admit if sample1==1, abs(year) cl(id) first
			estadd local dFE "\checkmark", replace
			estadd local tFE "\checkmark", replace
			estadd local oFE "\checkmark", replace	

	estout m1 m2 m3 m4 m5 m6 m7 m8 ///
	using "${TAB}\`tex_file_name'", replace ///
	drop(_cons 1.income_gp 1.income_gp#c.l_p admit) ///
	order(lH_a l_p 2.income_gp#c.l_p 3.income_gp#c.l_p 4.income_gp#c.l_p l_DEN) ///
	style(tex) mlabel(none) collabel(none) label ///
	prehead("\begin{tabular}{lcccccccc}" "\toprule\toprule") ///
	posthead("`latex_header'" "\midrule" "\multicolumn{4}{l}{\textbf{`panel_text'}} \\") ///
	cells("b(fmt(%9.3f) star)" "se(fmt(%9.3f) par)")  ///
	starlevels(* 0.10 ** 0.05 *** 0.01) ///
	prefoot("\midrule")	
	
	
	/* Panel B. Long run estimates */	
	local panel_text "Panel B- Long-term estimates"
*** Long term OLS 
		quiet reghdfe $yvar  $endog $controls if sample2==1, abs($fes) cl(id)
		eststo m1_lmic : nlcom (_b[l_p]) / (-_b[lH_a]), post
		quiet reghdfe $yvar  $endog $controls if sample2==1, abs($fes) cl(id)
		eststo m1_umic : nlcom (_b[l_p] + _b[3.income_gp#c.l_p]) / (-_b[lH_a]), post
		quiet reghdfe $yvar  $endog $controls if sample2==1, abs($fes) cl(id)
			estadd local dFE "\checkmark", replace
			estadd local tFE "\checkmark", replace
			estadd local oFE "-", replace
		eststo m1_hic : nlcom (0)/ (-_b[lH_a]), post
		
		
		quiet reghdfe $yvar  $endog $controls admit if sample2==1, abs($fes) cl(id)
		eststo m2_lmic : nlcom (_b[l_p]) / (-_b[lH_a]), post
		quiet reghdfe $yvar  $endog $controls admit if sample2==1, abs($fes) cl(id)
		eststo m2_umic : nlcom (_b[l_p] + _b[3.income_gp#c.l_p]) / (-_b[lH_a]), post
		quiet reghdfe $yvar  $endog $controls admit if sample2==1, abs($fes) cl(id)
		eststo m2_hic : nlcom (0) / (-_b[lH_a]), post
			estadd local dFE "\checkmark", replace
			estadd local tFE "\checkmark", replace
			estadd local oFE "\checkmark", replace
		
		quiet reghdfe $yvar  $endog $controls if sample1==1, abs($fes) cl(id)
		eststo m3_lmic : nlcom (_b[l_p]) / (-_b[lH_a]), post
		quiet reghdfe $yvar  $endog $controls if sample1==1, abs($fes) cl(id)
		eststo m3_umic : nlcom (_b[l_p] + _b[3.income_gp#c.l_p]) / (-_b[lH_a]), post
		quiet reghdfe $yvar  $endog $controls if sample1==1, abs($fes) cl(id)
			estadd local dFE "\checkmark", replace
			estadd local tFE "\checkmark", replace
			estadd local oFE "-", replace
		eststo m3_hic : nlcom (_b[l_p] + _b[4.income_gp#c.l_p])/ (-_b[lH_a]), post
		
		quiet reghdfe $yvar  $endog $controls admit if sample1==1, abs($fes) cl(id)
		eststo m4_lmic : nlcom (_b[l_p]) / (-_b[lH_a]), post
		quiet reghdfe $yvar  $endog $controls admit if sample2==1, abs($fes) cl(id)
		eststo m4_umic : nlcom (_b[l_p] + _b[3.income_gp#c.l_p]) / (-_b[lH_a]), post
		quiet reghdfe $yvar  $endog $controls admit if sample2==1, abs($fes) cl(id)
			estadd local dFE "\checkmark", replace
			estadd local tFE "\checkmark", replace
			estadd local oFE "\checkmark", replace
		eststo m4_hic : nlcom (_b[l_p] + _b[3.income_gp#c.l_p]) / (-_b[lH_a]), post
		
*** Long term IV 
		quiet ivreghdfe $yvar  ($endog = $instr) $controls if sample2==1, abs($fes) cl(id) first
		eststo m5_lmic : nlcom( _b[l_p]) / (-_b[lH_a]), post
		quiet ivreghdfe $yvar  ($endog = $instr) $controls if sample2==1, abs($fes) cl(id) first
		eststo m5_umic : nlcom (_b[l_p] + _b[3.income_gp#c.l_p]) / (-_b[lH_a]), post
		quiet ivreghdfe $yvar  ($endog = $instr) $controls if sample2==1, abs($fes) cl(id) first
			estadd local dFE "\checkmark", replace
			estadd local tFE "\checkmark", replace
			estadd local oFE "-", replace
		eststo m5_hic : nlcom (0) / (-_b[lH_a]), post
	
		quiet ivreghdfe $yvar  ($endog = $instr) $controls admit if sample2==1, abs($fes) cl(id) first
		eststo m6_lmic : nlcom (_b[l_p]) / (-_b[lH_a]), post
		quiet ivreghdfe $yvar  $endog $controls admit if sample2==1, abs($fes) cl(id) first
		eststo m6_umic : nlcom (_b[l_p] + _b[3.income_gp#c.l_p]) / (-_b[lH_a]), post
		quiet ivreghdfe $yvar  $endog $controls admit if sample2==1, abs($fes) cl(id) first
			estadd local dFE "\checkmark", replace
			estadd local tFE "\checkmark", replace
			estadd local oFE "\checkmark", replace
		eststo m6_hic : nlcom (0) / (-_b[lH_a]), post
		
		quiet ivreghdfe $yvar  ($endog = $instr) $controls if sample1==1, abs($fes) cl(id) first
		eststo m7_lmic : nlcom (_b[l_p]) / (-_b[lH_a]), post
		quiet ivreghdfe $yvar  ($endog = $instr) $controls if sample1==1, abs($fes) cl(id) first
		eststo m7_umic : nlcom (_b[l_p] + _b[3.income_gp#c.l_p]) / (-_b[lH_a]), post
		quiet ivreghdfe $yvar  ($endog = $instr) $controls if sample1==1, abs($fes) cl(id) first
			estadd local dFE "\checkmark", replace
			estadd local tFE "\checkmark", replace
			estadd local oFE "-", replace
		eststo m7_hic : nlcom (_b[l_p] + _b[4.income_gp#c.l_p]) / (-_b[lH_a]), post
		
		quiet ivreghdfe $yvar  ($endog = $instr) $controls admit if sample1==1, abs($fes) cl(id) first
		eststo m8_lmic : nlcom (_b[l_p]) / (-_b[lH_a]), post
		quiet ivreghdfe $yvar  ($endog = $instr) $controls admit if sample1==1, abs($fes) cl(id) first
		eststo m8_umic : nlcom (_b[l_p] + _b[3.income_gp#c.l_p]) / (-_b[lH_a]), post
		quiet ivreghdfe $yvar  ($endog = $instr) $controls admit if sample1==1, abs($fes) cl(id) first
			estadd local dFE "\checkmark", replace
			estadd local tFE "\checkmark", replace
			estadd local oFE "\checkmark", replace	
		eststo m8_hic : nlcom (_b[l_p] + _b[4.income_gp#c.l_p]) / (-_b[lH_a]), post		
		
	estout m1_lmic m2_lmic m3_lmic m4_lmic m5_lmic m6_lmic m7_lmic m8_lmic ///
	using "${TAB}\`tex_file_name'", append ///
	style(tex) mlabel(none) collabel(none) label ///
	prehead(" ") ///
	posthead("\multicolumn{9}{l}{\textbf{`panel_text'}} \\") ///
	cells("b(fmt(%9.3f) star)" "se(fmt(%9.3f) par)")  ///
	varlabels(_nl_1  "Low and Lower-Middle", elist(_nl_1)) ///
	starlevels(* 0.10 ** 0.05 *** 0.01) ///
	prefoot("\midrule") ///	

	
	estout m1_umic m2_umic m3_umic m4_umic m5_umic m6_umic m7_umic m8_umic ///
	using "${TAB}\`tex_file_name'", append ///
	style(tex) mlabel(none) collabel(none) label ///
	prehead(" ") ///
	cells("b(fmt(%9.3f) star)" "se(fmt(%9.3f) par)")  ///
	varlabels(_nl_1  "Upper-Middle", elist(_nl_1)) ///
	starlevels(* 0.10 ** 0.05 *** 0.01) ///	
	
	estout m1_hic m2_hic m3_hic m4_hic m5_hic m6_hic m7_hic m8_hic ///
	using "${TAB}\`tex_file_name'", append ///
	style(tex) mlabel(none) collabel(none) label ///
	prehead(" ") ///
	cells("b(fmt(%9.3f) star)" "se(fmt(%9.3f) par)")  ///
	varlabels(_nl_1  "High-Income", elist(_nl_1)) ///
	starlevels(* 0.10 ** 0.05 *** 0.01) ///
	
	
	estout m1 m2 m3 m4 m5 m6 m7 m8 /// 
	using "${TAB}/`tex_file_name'", append ///
	drop(_cons 1.income_gp#c.l_p 2.income_gp#c.l_p 3.income_gp#c.l_p 4.income_gp#c.l_p admit lH_a l_DEN l_p 1.income_gp 2.income_gp 3.income_gp 4.income_gp) ///
	style(tex) mlabel(none) collabel(none) label ///
	prehead(" ") ///
	cells("b(fmt(%9.3f) star)" "se(fmt(%9.3f) par)")  ///
	starlevels(* 0.10 ** 0.05 *** 0.01) ///
	stats(N r2 tFE oFE rkf,  labels("N" "\$R^2\$" "Decade FE" "Outlier FE" "K. Paap F-Stat" ) ///
	fmt(%9.0fc %9.3fc %9.0g %10.0f %10.0g %10.2fc)) ///
	prefoot("\midrule") ///
	postfoot("\bottomrule\bottomrule" "\end{tabular}")  
	
	
	
/*-----------------------------------------------------------------------------*
 Figure 1: Effect of selective emigration on human capital accumulation (hi).
				Insights from the generalized approach
*-----------------------------------------------------------------------------*/


***	Data preparation 

	preserve 
	import delimited "$dataMATLAB\Benchmark_IER_25.csv", clear 
	keep iso nir6_nm nir6_pesc nir6_sc nir7_nm nir7_pesc nir7_sc
	ren (nir6_nm nir6_pesc nir6_sc) (nir6_nmbenc nir6_pescbenc nir6_scbenc)
	ren (nir7_nm nir7_pesc nir7_sc) (nir7_nmbenc nir7_pescbenc nir7_scbenc)
	g id = _n  
	tempfile Benchmark 
	save `Benchmark'
	restore 
	
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
	
	preserve 
	import delimited "$dataMATLAB\HumanCap_Bench_25.csv", clear case(preserve)
	g id = _n  
	tempfile Bench 
	save `Bench'
	restore 
	
	preserve 
	import delimited "$dataMATLAB\HumanCap_SR_25.csv", clear case(preserve)
	g id = _n  
	tempfile ShortRun 
	save `ShortRun'
	restore 
	

	use `metrics', clear 
	merge 1:1 id using `Bench'
	drop _m 
	merge 1:1 id using `Benchmark'
	drop _m
	merge 1:1 id using `ShortRun'
	drop _m 
	
	
	g dLambNM = Lambda/LambdaNM-1
	g hti = Niih/(Niih+Niil)
	g dh = hti-HrNM
	

****	Figure 1 preparation:   
* 1-a   
* 1-b 
* 1-c 
* 1-d 

	kdens dLambNM, gen(d1 x1) ci(ci1_1 ci1_2) ///
		  xla(-.5(.5)2.0, labs(medium)) yla(0(1)5, labs(medium)) lcolor(black) ///
		  lwidth(thick) xli(0, lpattern(dash) lcolor(red)) ///
		  xti("Relative deviation in expected education premium", size(large)) ///
		  yti("Density", size(large)) ///
		  legend(position(6) cols(2) ///
		  label(1 "95% CI") ///
		  label(2 "Kernel density") size(large))
	graph export "$FIG\Figure_1a.png", replace 
	

	twoway lfitci dLambNM LambdaNM, legend(off) lcolor(red) alwidth(none) || ///
	       scatter dLambNM LambdaNM, ms(none) ml(iso) ///
			   mlabpos(0) mlabcolor(navy) mlabs(vsmall) ///
			   xti("No-emigration expected education premium", size(large)) ///
			   yti("Relative deviation in education premium", size(large)) ///
			   xla(1(1)4, labs(medium)) yla(-.5(.5)2.0, labs(large))
	graph export "$FIG\Figure_1b.png", replace 	   
	   

	scatter dh HrNM [aw=LNM], ms(Oh) mcolor(navy) mlw(medthick) ///
			   xti("No-emigration human capital for resident (LR)", size(large)) ///
			   yti("Absolute deviation from closed economy", size(large)) ///
			   yli(0, lpattern(dash) lcolor(red)) ///
			   yla(-.04(.02).08, labs(medium)) xla(,labs(large))
	graph export "$FIG\Figure_1c.png", replace 	
	
	
	scatter hti HrNM [aw=LNM], ms(Oh) mcolor(navy) mlw(medthick) || ///
				function y=x, range(0 .5) ///
			   xti("No-emigration human capital (LR)", size(large)) ///
			   yti("Observed human capital"	, size(large)) leg(off) ///
			   xla(,labs(medium)) yla(,labs(large))
	graph export "$FIG\Figure_1d.png", replace 	
	
	
	
/*-----------------------------------------------------------------------------*
 Figure 2: Effect of selective emigration on disposable income (yi)
*-----------------------------------------------------------------------------*/



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
	
	
	use `Wages', clear 
	merge 1:1 id using `Benchmark'
	drop _m	
	merge 1:1 id using `Channels'
	drop _m 
	

***	Figure preparation  


/*
	Figure 2.a
*/
	kdens NIR6_NMb, ///
		  gen(d1 x1) ci(ci1_1 ci1_2) ///
		  xla(-.2(.1).5, labs(medium)) yla(0(5)15, labs(medium)) lcolor(black) ///
		  lwidth(thick) xli(0, lpattern(dash) lcolor(red)) ///
		  xti("Net per worker disposable income response", size(medium)) ///
		  yti("Density", size(medium)) ///
		  legend(position(6) cols(2) ///
		  label(1 "95% CI") ///
		  label(2 "Kernel density") size(medium))
	graph export "$FIG\Figure_2a.png", replace 	
		  
		  
/*
	Figure 2.b
*/	

twoway kdensity HumCap_b, bwidth(.004) lpattern(solid) lcolor(black) lwidth(medthick) || ///
	   kdensity TecExt_b, bwidth(.003) lpattern(dash) lcolor(black) lwidth(medthick)  || ///
	   kdensity DiaExt_b, bwidth(.004) lpattern(solid) lcolor(blue) lwidth(medthick)  || ///
	   kdensity FisExt_b, bwidth(.010) lpattern(vshortdash) lcolor(black) lwidth(medthick)  || ///
	   kdensity MktExt_b, bwidth(.004) lpattern(solid) lcolor(black)  lwidth(small)  || ///
	   kdensity RemEff_b, bwidth(.012) lpattern(dash) lcolor(blue) lwidth(medthick)  ///
			 xti("Percentage change in per worker disposable income", size(large)) ///
			 yti("Density", size(large)) ///
			 xlab(-.1(.1).4, labs(large)) ///
			 xline(0, lcolor(red) lpattern(dash) lwidth(.3pt)) ///
			 note("") title("") ylab(, labs(large)) legend(ring(0) pos(2) cols(1) ///
			 label(1 "Human capital") ///
			 label(2 "Technological ext.") ///
			 label(3 "Diaspora ext.") ///
			 label(4 "Fiscal ext.") ///
			 label(5 "Market size ext.") ///
			 label(6 "Remittances") size(large))
	graph export "$FIG\Figure_2b.png", replace
	
	
/*
	Figure 2.c
*/		

g ywd_nm = 	ywd/(1+NIR6_NMb)
g lywd_nm = ln(ywd_nm)

	twoway lfitci NIR6_NMb lywd_nm, legend(off) lcolor(red) alwidth(none) || ///
	       scatter NIR6_NMb lywd_nm [aw=LPE], ms(Oh) mcolor(navy) mlw(medthick) ///
			   mlabpos(0) mlabcolor(navy) mlabs(large) ///
			   xti("Log disposable income per worker (no-emigration)", size(large)) ///
			   yti("Net per worker disposable" "income response", size(large)) ///
			   xla(,labs(medium)) yla(-.2(.1).5, labs(large)) yli(0, lpattern(shortdash))
	graph export "$FIG\Figure_2c.png", replace 	


/*
	Figure 2.d
*/	

	twoway (kdensity NIR6_NMb, bwidth(.020) lcolor(black) lpattern(solid) lwidth(thick)) ///
		   (kdensity NIR7_NMb, bwidth(.060) lcolor(black) lpattern(dash) lwidth(thick) ///
		    legend( ring(0) position(1) cols(1) ///
			label(1 "Per worker") ///
			label(2 "Per natural") size(large)) ///
			xla(-.2(.2)1, labs(medium)) yla(0(2)10, labs(large)) ///
			xti("Net per worker disposable income response", size(large)) ///
			yti("Density", size(large)) ///
			xli(0, lpattern(shortdash) lcolor(red)))
	graph export "$FIG\Figure_2d.png", replace 			
	
	
	

/*-----------------------------------------------------------------------------*
 Figure 3: Average disposable income responses to selective emigration
*-----------------------------------------------------------------------------*/	

*** Data preparation 
	
	import delimited "$dataMATLAB\AVG_NIR_Bench_25.csv", clear case(preserve)  
	g devcat = _n 
	
	lab define devcat 1 "World" 2 "LIC" 3 "LMIC" 4 "UMIC" 5 "HIC"
	lab values devcat devcat 	
	foreach var in AVG_NIR6_NM AVG_NIR7_NM {
		replace `var' = 100*`var'
	}
	
*** Figure preparation 
	graph bar AVG_NIR6_NM AVG_NIR7_NM, over(devcat, label(labs(large))) yla(0(5)30, labs(large)) ///
		  bar(1, fcolor(navy) fintensity(50) color(none)) ///
		  bar(2, fcolor(maroon) fintensity(50) color(none)) ///
		  legend(position(6) cols(2) label(1 "Per worker") label(2 "Per natural") size(large)) ///
		  yti("Net disposable income response", size(large))  
	graph export "$FIG\Figure_3.png", replace 			
	

	
	

/*------------------------------------------------------------------------------*
	Figure 4: Emigration and world income distribution				   *					 
*------------------------------------------------------------------------------*/

*** Data preparation 
	
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

		
*** Figure preparation 			
	twoway (kdensity lnwdd [fw=La], bw(.17) lcolor(black) lpattern(solid) lwidth(medthick)) ///
		   (kdensity lnwSCd [fw=La], bw(.17) lcolor(black) lpattern(dash) lwidth(medthick)) ///
		   (kdensity lnwSCd [fw=LSCa], bw(.17) lcolor(blue) lpattern(solid) lwidth(medthick) /// 
		    legend(position(6) cols(3) label(1 "Observed") ///
									   label(2 "NM(Initial labor force)") ///
									   label(3 "NM(New labor force)") size(large)) ///
			xla(7(1)13, labs(medium)) xti("Disposable income per worker", size(large)) yla(0(.2).8, labs(large)) /// 
			yti("Density", size(large)) xli(7.603 9.426, lpattern(shortdash) lcolor(red)) ///
			text(0.3 7.47 "International poverty line ($2,004)", place(e) orientation(vertical) size(.3cm)) ///
			text(0.25 9.3 "Median disposable income ($12,404)", place(e) orientation(vertical) size(.3cm)))
	graph export "$FIG\Figure_4.png", replace 			
		
