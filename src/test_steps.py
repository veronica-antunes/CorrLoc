import matplotlib.pyplot as plt
import numpy as np
import matplotlib, os
from obspy import read
import math as m

def plot_wavs(wav_name='filtered_wavs.mseed', fig_height=4, fig_width=2.3,
	station_list=None, font_size=10, outfile='wavs.png', verbose=True ) :
	'''
	Used to plot the waveforms
	Example
	plot_wavs(wav_name='filtered_wavs.mseed', fig_height=4, fig_width=2.3,
	station_list=None, font_size=10, outfile='wavs.png', verbose=True)

	possible to plot specific stations by including a station_list with the names of the 
	stations that we want to plot.
	fig_height and fig_width serve to control the size of the plot
	'''

	from matplotlib.pyplot import figure
	figure(num=None, figsize=(4, 3.5), dpi=300)
	matplotlib.rcParams.update({'font.size': font_size-2})

	#height=4
	#width=2.3
	#wav_name='filtered_wavs.mseed'
	#station_list=['NIR01','NIR03','NIR04','NIR05','NIR06']
	#colours=['b','forestgreen','r','magenta','BlueViolet','SaddleBrown', 'grey', 
	#	'darkorange', 'y', 'lightgreen','cyan','lemongreen','white','k']
	#fs=10

	wavs=read(wav_name)

	if station_list :
		new_wavs=wavs.select(station=station_list[0])
		for n in range(1,len(station_list)):
			new_wavs+=wavs.select(station=station_list[n])
		plot_size=float(len(station_list))
	else:
		new_wavs=wavs
		plot_size=float(len(wavs))

	#new_wavs.trim(starttime=wavs[0].stats.starttime+12*3, endtime=wavs[0].stats.starttime+12*4)
	fig = plt.figure()

	count=0
	for w in new_wavs:
		ax0 = fig.add_subplot(plot_size, 1 ,count+1)

		tr=w.data
		samp_rate=w.stats.sampling_rate

		x=np.linspace(0,len(tr)/samp_rate,len(tr))
		
		#ax0.plot(x,tr, color=colours[int(w.stats.station.split('L1')[1])-1], label=w.stats.station)
		#stat='N'+str(int(w.stats.station[-2::]))
		stat=w.stats.station
		ax0.plot(x, tr, color='b', linewidth=0.5,label=stat)

		ax0.set_xlim(0,int(max(x)))
		
		maxy=max( max(tr), abs(min(tr)) )
		starty, endy = -maxy, maxy
		
		#if maxy < 1:
		#	starty, endy = -1, 1	

		ax0.set_ylim(starty, endy)
		#starty, endy = ax0.get_ylim()
		ax0.yaxis.set_ticks(np.round(np.linspace(starty, endy, 3)[1::]))

		startx, endx = ax0.get_xlim()
		ax0.xaxis.set_ticks(np.linspace(startx, endx, 4))

		ax0.set_title(stat, fontsize=font_size,
			x=0.85, y=0.69-plot_size*0.02)

		if count < plot_size-1 :
			#to remove the labels on top figures
			ax0.axes.get_xaxis().set_ticklabels([]) 
		if count ==  plot_size-1 :
			ax0.set_xlabel('[s]', fontsize=font_size)

		count=count+1

	fig.set_figheight(fig_height)
	fig.set_figwidth(fig_width)
	fig.text(-0.00,0.5, 'Normalizeed Amplitude', ha='center', va='center', 
		rotation=90, fontsize=font_size)
	#plt.suptitle('Window '+str(w))
	#pdf=plt
	#pdf.savefig('wavs.pdf', bbox_inches='tight' )

	plt.savefig(outfile, bbox_inches='tight' )
	plt.close()

def plot_crosscorr(np_pair='pair.npy', np_xcorr='xcorr_struct.npy', np_env='env_struct.npy',
	waveform='filtered_wavs.mseed', windows='all' ,pairs='all', equal_scale= True,
	outpath='2.CrossCorr/', format='png', fontsize=12, pick=0, verbose=True) :
	'''
	Plots the cross-correlation function and envelope of the different station pairs 
	for the different time windows
	Example:
	test_step2 (np_pair='pair.npy', np_xcorr='xcorr_struct.npy', np_env='env_struct.npy',
		waveform='filtered_wavs.mseed', windows='all', pairs='all', equal_scale=True,
		pick=0.0, outpath='test_step2/', format='png', fontsize=12, verbose=True)
	
	windows and pairs should be "all" or a list of numbers (corresponding to the
	respective pairs and window number that we want to plot)
	picks will draw a red line in the respective x-value defined. if None, it will not draw
	the red line
	It is possible to use the same equal scale for all plots if equal_scale is True. If False it
	will draw individual y limits for each plot.
	np_pair, np_xcorr and np_env are outputs from the location code algorithm
	waveform is used to get the sampling rate from the traces (all traces should have the same)
	outpath is the folder where the figures will be saved 
	format is the format of the output figure, ex: 'png', 'pdf'...
	fontize default is 12 and corresponds to the fontsize of the axis. all the other sizes will be
	interpolated from this value (12+2 for labels and 12+3 for title)
	'''
	#from math import sin, cos
	#from __future__ import division # to allow division using /	

	pair=np.load(np_pair)
	xcorr_struct=np.load(np_xcorr)
	env_struct=np.load(np_env)
	amax=env_struct.max()
	#senv_struct=np.load('senv_struct.npy')

	wavs=read(waveform)
	samp_rate=wavs[0].stats.sampling_rate	#sample rate in Hz

	if windows=='all':
		wind_range=range(env_struct.shape[1])
	else:
		try:
			#to test if input is a list
			wind_range=len(windows)
			wind_range=windows
		except:
			raise ValueError('windows must be a list of integers!')
			#print('windows must be a list of integers!')
			#sys.exit()

	if pairs=='all':
		pairs_range=range(env_struct.shape[2])
		plot_size=int(m.ceil(env_struct.shape[2]/2))

	else:
		try:
			#to test if input is a list
			pairs_range=pairs
			plot_size=int(m.ceil(len(pairs)/2))
		except:
			raise ValueError('pairs must be a list of integers!')
			#print('pairs must be a list of integers!')
			#sys.exit()

	#if not (plot_size/2).is_integer():
	#	plot_size=plot_size+1
	
	for w in wind_range:
		plt.close('all')
		fig = plt.figure()

		if verbose:
			print('plotting CC and envelope values for window '+str(w))

		count=0
		for p in pairs_range :
			try :
				xcorr1=xcorr_struct[:,w,p] #326
				env1=env_struct[:,w,p]	    #326
			except:
				break

			ax0 = fig.add_subplot(plot_size,2,count+1)

			new_xlabel=np.linspace(-(len(env1)/2)/samp_rate, (len(env1)/2)/samp_rate, len(env1))

			ax0.plot(new_xlabel,xcorr1, label='p'+str(p)+': '+str(pair[count][1])+'-'+str(pair[count][2]), linewidth=0.5)
			ax0.plot(new_xlabel,env1, '-', linewidth=0.8) #, label='env')
			#ax0.plot(senv1, '-')#, label='env')

			#ax0.set_ylim(-1, 1)
			ax0.set_title('p'+str(p)+': '+str(pair[count][1])+'-'+str(pair[count][2]), fontsize=fontsize,
				x=0.50, y=0.96 )
			#ax0.set_title('pair '+str(p)+', segment 62')

			#starty, endy = ax0.get_ylim()
			if equal_scale:
				ax0.yaxis.set_ticks(np.linspace(-amax, amax, 3))
			
			plt.xticks(np.linspace(-(len(env1)/2)/samp_rate, (len(env1)/2)/samp_rate, 5), size=fontsize)
			#ax0.yaxis.set_ticks(np.arange(-1, 1+0.1, 1 ))
			#ax0.yaxis.set_ticks(np.arange(starty, endy+0.1, 1 ))
			ax0.axes.get_yaxis().set_ticklabels([]) #to remove the ticks

			if pick != None :
				plt.axvline(x=pick, 
					color='r',
					linestyle='-', alpha=0.5) # ,
					#label='yenv= '+str(round(yenv,5)))

			if count == int(env_struct.shape[2])-2 or count == int(env_struct.shape[2])-1 :
				ax0.set_xlabel('time [s]', fontsize=fontsize+2)

			else :
				ax0.axes.get_xaxis().set_ticklabels([]) #to remove the
									  #labels on top figures
			count=count+1
		
		fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.15, hspace=0.4)
		fig.text(0.09,0.5, 'CC value normalized', ha='center', va='center', rotation=90, fontsize=fontsize+2)
		plt.suptitle('window '+str(w), size=fontsize+3)
		#pdf=plt

		if not os.path.exists(outpath):
			os.makedirs(outpath)

		plt.savefig(outpath+'/W'+str(w)+'_P'+str(pairs)[1:-1]+'.'+format, bbox_inches='tight' )
		#pdf.savefig('test_step2/'+str(w)+'.pdf', bbox_inches='tight' )
		plt.close()
		#plt.show()

def plot_difftime(arrayA='A_array_0380.npy', np_stat='stat_struct.npy', mesh='mesh_struct.npy',
	np_pair='pair.npy', pairs='all', fontsize=11, outpath='3.DiffTime/', format='png', 
	verbose=True) :
	'''
	plots the velocity from each grid point to each stat pair
	Example:
	test_step3(arrayA='A_array_0380.npy', np_stat='stat_struct.npy', mesh='mesh_struct.npy',
		np_pair='pair.npy', pairs='all', fontsize=11, outpath='test_step3/', format='png', 
		verbose=True)
	'''

	A=np.load(arrayA)
	stat_net=np.load(np_stat)
	mesh=np.load(mesh)
	pair_list=np.load(np_pair)


	x_mesh=mesh[0]
	y_mesh=mesh[1]
	#plt.plot(x_mesh, y_mesh, '.', color='0.75') 


	x_stat=stat_net[:,1]
	y_stat=stat_net[:,2]
	name_stat=stat_net[:,0]
	#plt.plot(x_stat, y_stat, '^', color='b') 

	#ax.set_xlim([0-50, 1000+50])
	#ax.set_ylim([0-50, 800+50])

	#step=10
	#corner=step/2
	if pairs == 'all':
		pairs_range= range(len(pair_list))
	else:
		try:
			pairs_range=len(pairs)
			pairs_range=pairs
		except:
			raise ValueError('pairs must be a list of integers!')
			#print('pairs must be a list of integers!')
			#sys.exit()		

	for p in pairs_range :
		if verbose: 
			print ('plotting pair... ', p)

		img=A[:,:,p]	#results for pair 2

		sta1=pair_list[p][1]
		sta2=pair_list[p][2]

		ax = plt.subplot(111)

		# for more color-maps check:
		# https://matplotlib.org/stable/tutorials/colors/colormaps.html
		# RdYlBu, RdBu, Seismic, 'Spectral', bwr, coolwarm

		img2=plt.imshow( img.T, cmap='RdBu',
			interpolation='none',	#interpolation: nearest or bicubic
			origin='lower',
			#aspect='auto',
			extent=[0,x_mesh.max(),0,y_mesh.max()] ) 
			#aspect='equal', extent=[0,1000,0,800] ) 
			#, clim=(-0.25, 0.25)) 

		#plt.plot(x_mesh, y_mesh, '.', color='0.75')
		for i in range(len(x_stat)) :
			if pair_list[p][1] == name_stat[i] or pair_list[p][2] == name_stat[i] :
				ax.scatter(float(x_stat[i]), float(y_stat[i]), marker='^', 
					color='k', s=80, edgecolor='w',linewidth=2)
			else :	
				ax.scatter(float(x_stat[i]), float(y_stat[i]), marker='v', 
					color='w', s=50, edgecolor='black',linewidth=0.5)
			
		ax.set_xlim([0, x_mesh.max()])
		ax.set_ylim([0, y_mesh.max()])

		plt.xticks(size=fontsize)
		plt.yticks(size=fontsize)

		plt.xlabel('distance [m]',size=fontsize+2)
		plt.ylabel('distance [m]',size=fontsize+2)

		#plt.colorbar()
		
		colorticks=np.linspace(img.min(),img.max(),7)

		cbar=plt.colorbar(ticks=colorticks)
		cbar.ax.set_ylabel('time [s]',size=fontsize+2)

		plt.title('Diferential time\nStation pair '+
			str(p)+': '+sta1+'-'+sta2, fontsize=fontsize+3)
		#pdf=plt
		
		if not os.path.exists(outpath):
			os.makedirs(outpath)
		
		plt.savefig(outpath+'/pair_'+str(p)+'.'+format) #, bbox_inches='tight' )
		#plt.savefig('test_step3/plot_step10_pair_'+str(p)+'_v250'+'.pdf') #, bbox_inches='tight' )
		plt.close('all')
		plt.close()
		#plt.show()

def plot_CCpairs(arrayA='A_array_0380.npy', arrayB='B_array_0380.npy',
	np_stat='stat_struct.npy', np_env='env_struct.npy', mesh='mesh_struct.npy', 
	np_pair='pair.npy', window=0, fontsize=12, pairs='all', verbose=True, 
	plot_point='max', outpath='4.CCpairs/',format='png', tw=None) :
	'''
	Plots the individual back-projections for each station pair
	Example:
	est_step4 (arrayA='A_array_0380.npy', np_stat='stat_struct.npy', 
	mesh='mesh_struct.npy', np_pair='pair.npy', arrayB='B_array_0380.npy',
	window=0, fontsize=12, pairs='all', outpath='test_step4/',format='png',
	plot_point = 'max', verbose=True, tw=None,
	np_env='env_struct.npy')

	arrayB can be B or D. B plots the individual back-projections normalized
	D plots the individual back-projections (every plot with its maximum value)
	plot_point can be 'max', None or a list of integers with the 4 parameters:
	[n1,n2,n3,n4], where, n1 and n2 is the grid point, n3 is the pair number and
	w is the window number
	wav_name is the waveform file to plot the envelope-amplitude plot or None
	if plot_point is None
	'''
	
	from matplotlib.colors import LinearSegmentedColormap
	matplotlib.rcParams.update({'font.size': fontsize})
	#import location

	#window=1
	#fs=18

	vel=float(arrayA.split('/')[-1].split('_')[2].split('.')[0])
	A=np.load(arrayA)
	stat_net=np.load(np_stat)
	mesh=np.load(mesh)
	pair_list=np.load(np_pair)
	env_struct=np.load(np_env)
	B=np.load(arrayB)
	#B=np.load('D_array_'+str(vel)+'.npy')

	amax=B[:,:,window,:].max() #max for plotting
	omg=B[:,:,window,:]/amax

	x_mesh=mesh[0]
	y_mesh=mesh[1]
	name_stat=stat_net[:,0]

	#plt.plot(x_mesh, y_mesh, '.', color='0.75') 

	x_stat=stat_net[:,1]
	y_stat=stat_net[:,2]
	#plt.plot(x_stat, y_stat, '^', color='b') 

	#ax.set_xlim([0-50, 1000+50])
	#ax.set_ylim([0-50, 800+50])
	#colours2=['w','w','w','w','w','w', 'w', 'w', 'w', 'w','w','w','white','w']

	plt.close('all')
	
	#corner=1
	if pairs == 'all':
		pairs_range= range(len(pair_list))
	else:
		try:
			pairs_range=len(pairs)
			pairs_range=pairs
		except:
			raise ValueError('pairs must be a list of integers!')
			#print('pairs must be a list of integers!')
			#sys.exit()		

	for p in pairs_range :
		if verbose:
			print ('plotting pair... ', p)

		sta1=pair_list[p][1]
		sta2=pair_list[p][2]

		#img=B[:,:,window,p]	#results for pair 2
		img=omg[:,:,p]
		ax = plt.subplot(111)

		vmax=1
		my_cmap = LinearSegmentedColormap.from_list('mycmap', 
			[(0 / vmax, 'darkblue'),
			(0.50 / vmax, 'blue'), 
			(0.75 / vmax, 'cyan'),
			(0.80 / vmax, 'lightgreen'),
			(0.85 / vmax, 'yellow'),
			(0.90 / vmax, 'orange'),
			(0.95 / vmax, 'red'),
			(0.99 / vmax, 'darkred'),
			(1.00 / vmax, 'black')])

		img2=plt.imshow( img.T,
			cmap=my_cmap, clim=(0,1), #clim=(0, 1),
			interpolation='none',	#interpolation: none, nearest or bicubic
			origin='lower',
			#aspect='auto',
			extent=[0,x_mesh.max(),0,y_mesh.max()] ) 
			#aspect='equal', extent=[0,1000,0,800] )
		 
		#plt.plot(x_mesh, y_mesh, '.', color='0.75')
		#plt.plot(x_stat, y_stat, '^', color='b', s=80)

		plt.xlabel('distance [m]')
		plt.ylabel('distance [m]')

		ax.set_xlim([0, x_mesh.max()])
		ax.set_ylim([0, y_mesh.max()])
		ax.set_title('pair '+str(p)+': '+sta1+'-'+sta2+'\nwindow '+
			str(window)+'  vel: '+str(vel))
		colorticks=np.linspace(0,1,6)

		cbar=plt.colorbar(ticks=colorticks)
		cbar.ax.set_ylabel('individual back-projections \n normalized')

		for i in range(len(x_stat)) :
			if pair_list[p][1] == name_stat[i] or pair_list[p][2] == name_stat[i] :
				ax.scatter(float(x_stat[i]), float(y_stat[i]), marker='^', 
					color='k', s=80, edgecolor='w',linewidth=2)
			else :	
				ax.scatter(float(x_stat[i]), float(y_stat[i]), marker='v', 
					color='w', s=50, edgecolor='black',linewidth=0.5)

		if not os.path.exists(outpath):
			os.makedirs(outpath)

		#pdf=plt
		#pdf.savefig('test_step4/stack_pair_'+str(p)+'_yval_'+str(vel)+'.pdf', bbox_inches='tight' )	
		plt.savefig(outpath+'/stk_'+'w'+str(window)+'_p'+str(p)+'_V'+str(vel)+'.'+format, bbox_inches='tight' )
		plt.close('all')
		plt.close()
		#plt.show()

	#plot point of max correlation
	if plot_point == None:
		pass
	elif plot_point == 'max' :
		n1,n2,n3,n4=np.where(B == B.max())
	elif len(plot_point) == 4 :
		try:
			n1,n2,n3,n4=plot_point
		except:
			print('plot_point must be a list of 4 integers, "max", or None')
			raise ValueError('to create the list use: [i,j,p,w] (i,j-position, p-pair, w-window)')
			#print('to create the list use: [i,j,p,w] (i,j-position, p-pair, w-window)')
			#sys.exit()
	else:
		print('plot_point must be a list of 4 integers, "max", or None')
		raise ValueError('to create the list use: [i,j,p,w] (i,j-position, p-pair, w-window)')
		#print('to create the list use: [i,j,p,w] (i,j-position, p-pair, w-window)')
		#sys.exit()

	if n1.any() and tw:
		for i in range(len(n1)):
			#pos=[int(n1[i]),int(n2[i]),int(n4[i]),int(n3[i])]
			pos_name=str(n1[i])+'_'+str(n2[i])+'_'+str(n4[i])+'_'+str(n3[i])
			print('Plotting cross correlation function for grid point ('+str(n1[i])+','+str(n2[i])+
				') window', str(n3[i]), 'and pair', str(n4[i]))

			nslide=int(tw/2)
			env=env_struct[:, n3[i], n4[i]]

			yenv1=B[ n1[i] , n2[i], n3[i] , n4[i]] 
			xenv=(np.abs(env - yenv1)).argmin() 
			yenv=env[xenv]

			plt.close('all')
			new_xlabel=np.linspace(-nslide, nslide, len(env))
			plt.plot(new_xlabel,env)
			new_xenv=new_xlabel[xenv] #(xenv-nslide)/samp_rate

			plt.axvline(x=new_xenv, color='g', linestyle='-', 
				label='xenv='+str(round(new_xenv,3)))
			
			plt.axhline(y=yenv, color='r', linestyle='-' , 
				label='yenv='+'{:.2e}'.format(yenv))

			if new_xenv >= 0 :
				plt.legend(loc="upper left", fontsize=12)
			if new_xenv < 0 :
				plt.legend(loc="upper right", fontsize=12)

			plt.title('Grid point ('+ str(n1[i])+','+
				str(n2[i])+	'), pair '+str(n4[i]) + 
				', window ' + str(n3[i]), size=14)

			plt.xlim(-nslide, nslide)

			plt.xlabel('time [s]', size=13)
			plt.ylabel('Env-CC value', size=13)

			plt.xticks(np.linspace(-nslide,	nslide, 5), size=12)
			plt.yticks(size=12)
			#plt.show()
			outfile=outpath+'/point_plot-'+pos_name+'.'+format
			plt.savefig(outfile, bbox_inches='tight')
				
			plt.close()



			#location.step4_envvalue(wav_name, point=pos, 
			#	outfile=outpath+'/point_plot-'+pos_name+'.'+format, verbose=False)
