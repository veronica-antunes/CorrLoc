import src.CorrLoc as CorrLoc
import src.test_steps as test
from obspy import read, UTCDateTime

# veronica antunes (veronica.antunes@sed.ethz.ch)
# last modification 08/08/2022

######################### SET PARAMETERS #######################################
# paths
inpath='loc_example'
wavs='wavs/L1_drone_shot1.mseed' #'/L1_hammer_shot6.mseed'

# grid parameters [m]
[x0, x1] = [0, 250]
[y0, y1] = [0, 200]
step = 1

# time window length [s]
tw = 3 # in seconds # 2

# velocity values to test
vel1=150
vel2=160
vel_step=10

#filter
fmi = 5 #1
fma = 45 #15

################################################################################

#processing signal suggestion using obspy
filtered_wavs='filtered_wavs.mseed'
w=read(wavs)
w.detrend(type='demean')
w.taper(max_percentage=0.1, type='hann', side='both') #max_lenght=None
w.filter('bandpass', freqmin=fmi, freqmax=fma, corners=2, zerophase='True')
w.write(inpath+'/'+filtered_wavs, format='MSEED')

#plot the waveforms
test.plot_wavs(wav_name=inpath+'/'+filtered_wavs, fig_height=4, fig_width=2.3,
	station_list=None, font_size=10, outfile=inpath+'/wavs.png' )

#step 1: create grid
CorrLoc.grid('stations_file.txt', x0, x1, y0, y1, step, 
	outfile='grid.png', fontsize=13, path=inpath)

#step2: perform cross-correlation
CorrLoc.crosscorr(filtered_wavs, tw, path=inpath)

test.plot_crosscorr(np_pair=inpath+'/pair.npy', np_env=inpath+'/env_struct.npy',
	np_xcorr=inpath+'/xcorr_struct.npy', waveform=inpath+'/'+filtered_wavs, 
	windows='all',pairs=[36,3,31,4,38,32,42,19,30,27], equal_scale= True, fontsize=12,
	outpath=inpath+'/results/plot_crosscorr/', format='png',  pick=0)


#step3 compute diferential times
CorrLoc.difftime(vel1,vel2,vel_step, path=inpath)

test.plot_difftime(arrayA=inpath+'/A_array_0'+str(vel1)+'.npy', 
	np_stat=inpath+'/stat_struct.npy', mesh=inpath+'/mesh_struct.npy',
	np_pair=inpath+'/pair.npy', pairs=[31,36,4,30,32,30],  
	outpath=inpath+'/results/plot_difftime/', format='png', fontsize=11)


#step4 get envelope values for each station pair
CorrLoc.envvalue(filtered_wavs, point=None, indv_norm=False,
			outfile='/results/test_step4/point_plot.png', 
			verbose=True, path=inpath)
			
test.plot_envvalue(arrayA=inpath+'/A_array_0'+str(vel1)+'.npy', 
	np_stat=inpath+'/stat_struct.npy', 
	mesh=inpath+'/mesh_struct.npy', np_pair=inpath+'/pair.npy', 
	arrayB=inpath+'/B_array_0'+str(vel1)+'.npy',
	np_env=inpath+'/env_struct.npy', pairs=[31,36,4,30,32,30], 
	outpath=inpath+'/results/plot_envvalue/',format='png',
	window=0, fontsize=12, plot_point='max', tw=tw)


#step5 stack all individual back-projections
CorrLoc.stackall(final_norm=True, path=inpath)

#step6 plot solution
CorrLoc.plotall('results/Loc', filtered_wavs, source_loc_file='shots_file.txt', 
	source_loc='DRON1', path=inpath, final_loc_file='results/locations.list')


