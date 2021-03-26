This folder contains all code in the NA monsoon project, from processing data, running back trajectory model, analyzing model output, to making figures in the paper. 


Process data 
data_process/ 
	mls_monthly_reading.py
		Read MLS data and save as .nc file over North American region 
		Line 145: outfile, change save directory 
	MLS_GRID.py
		Gridding the MLS data 
	ERAi_tropopause.py
		Calculate tropopause height, and save as .nc file 
	ERAi_100.py 
		Interpolate wind and temperature field to 100 hPa 
	NEXRAD_reading.py
		Calculate the convection occurrence in NEXRAD data 
			Remember to change output file name in line 333


Files for back trajectory model 
traj_back/
	call_model.pro 
		The program that initiate 27 parcels around MLS observations and then call the back traj model
		Change: Line 2: infile 
			     line 3: savFolder (add _l when -5 in traj_erai_back_multiple.pro line 215, _h when +5, see the introduction of traj_erai_back_multiple.pro)
		 
	traj_erai_back_multiple.pro:
		The back trajectory model 
		Line 122: outDir: output file directory 
		Line 48: Recday: record x days 
		Line 215: oth = theta1x 
				 or oth = theta1x+5, or oth = theta1x-5
		Line 69: b_ERAi_run = 1  choose reanalysis data that drives the back trajectory model 
	
	call_model2: 
		The program that initiate parcels evenly over NA and call the back traj model 
			

Analyze the back trajectory model result 
back_traj_analysis/
	backtraj_calculation.py
		Read back trajectory model result, and compare it with deep convection location 
		Change input in first several lines 
		Output: NC file of 2-D parcels' trajectories and whether they encounter deep convection 
	back_traj_encounter.py
		Calculate and save the parcels that has encountered deep convection in the back trajectory model 
		Remember to change out file directory in line 139
	parcel_history.py 
		Calculate parcel' history 
	
	
Plotting
plot/ 
	figures_for_paper_v3.ipynbi
		Calculating the variables in figures and plot them out. 

