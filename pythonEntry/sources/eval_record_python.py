import xgboost as xgb
import numpy as np
import sys


def findDetails(S2):
	#This function finds the first five details of a signal(S2)
	H = np.array([0.125, 0.375, 0.375, 0.125])
	G = np.array([-2, 2])

	lambda_1 = np.array([1.5, 1.12, 1.03, 1.01, 1.00]);

	j = 1;
	numOfStages = 4;
	print "S2 is sdfjlsdflsjdf",S2[2323]
	W=np.zeros((numOfStages,len(S2)))
	#S2 = S;
	while(j<=numOfStages):
		W[j-1,:] = np.convolve(S2, G, mode='same')/lambda_1[j-1]
		S2 = np.convolve(S2, H, mode='same')
		j = j+1
	return W

#load file description matfile
bp_index = []
sampling_freq = 360
with open('../'+sys.argv[1]+'.hea') as fp:
	lineCount = 0;
	for line in fp:
		print line
		if(lineCount == 0):
			splitLine = line.split(" ")
			sampling_freq = int(splitLine[2])
		else :
			if (str.__contains__(line, 'BP') or str.__contains__(line, 'bp') or str.__contains__(line, 'ART') or str.__contains__(line, 'art') or str.__contains__(line, 'ABP') or str.__contains__(line, 'abp')):
				bp_index = lineCount
		lineCount = lineCount + 1
print bp_index
print "sampling freq is ", sampling_freq
#description = sio.loadmat('fileDescription.mat')

#bp_index = description['bp_index']
#it is 1 since we have converted .dat to .txt and python arrays are zero based
ecg_index = 1

testData = np.loadtxt('../'+sys.argv[1]+'.txt')
#testData = np.loadtxt('1003.txt')
#ECG xgboost prediction
result_ecg = []
if(not (np.std(testData[:,ecg_index]) == 0)):
	WINLEN_ecg = 50
	testData[:,ecg_index] = (testData[:,ecg_index] - np.mean(testData[:,ecg_index]))/np.std(testData[:,ecg_index])
	S_ecg = np.concatenate((np.zeros(WINLEN_ecg), testData[:,ecg_index], np.zeros(WINLEN_ecg)))
	W_ecg = findDetails(S_ecg)

	sig_len = len(S_ecg) - WINLEN_ecg*2

	#Setting up ip features
	ip = np.zeros((sig_len, WINLEN_ecg*4+2))

	for j in xrange(1, 2*WINLEN_ecg+2):
		ip[:,j-1] = W_ecg[3, j-1:j-1+sig_len]
		ip[:,j+2*WINLEN_ecg] = W_ecg[2,j-1:j-1+sig_len]
	#load ecg model trained
	model = xgb.Booster(model_file='ecg.model')
	#model = xgb.load('ecg.model')
	#predictions for ecg-xgboost
	xgtest = xgb.DMatrix(ip)
	result_ecg = model.predict(xgtest)
	del W_ecg
	del S_ecg
	np.savetxt('ecg_pred.txt',result_ecg)

#BP xgboost prediction
result_bp = []
if(bool(bp_index) and (np.std(testData[:,bp_index]) != 0)):
	WINLEN_bp = 100
	#bp_index=bp_index[0][0]
	testData[:,bp_index] = (testData[:,bp_index] - np.mean(testData[:,bp_index]))/np.std(testData[:,bp_index]);
	S_bp = np.concatenate((np.zeros(WINLEN_bp), testData[:,bp_index], np.zeros(WINLEN_bp)))
	W_bp = findDetails(S_bp)

	sig_len = len(S_bp) - WINLEN_bp*2

	#Setting up ip features
	ip = np.zeros((sig_len, WINLEN_bp*2+2))

	for j in xrange(1, WINLEN_bp+2):
		ip[:,j-1] = W_bp[3, 2*j-2:2*j+sig_len-2]
		ip[:,j+WINLEN_bp] = W_bp[2,2*j-2:2*j+sig_len-2]
	#load bp model trained
	model = xgb.Booster(model_file='bp.model')
	#predictions for bp-xgboost
	xgtest = xgb.DMatrix(ip)
	result_bp = model.predict(xgtest)
	np.savetxt('bp_pred.txt',result_bp)
###############################################
SHORTWIN = 10
LONGWIN = 30
#RES_TH = median(result_all)
CNT_TH = 5


sig_len = len(result_ecg)

#this is causing some problem - i think it shouldn't be done.
#res_smooth = np.convolve(result_ecg,np.ones(SHORTWIN+1)/(SHORTWIN+1))
res_smooth = result_ecg
isPeak = res_smooth > np.mean(res_smooth)
#isPeak = isPeak.astype(int)

#isPeak = res_smooth > 0.3;
print res_smooth
print res_smooth[0]
print res_smooth[-1]
res_smooth_d = np.concatenate((np.array([res_smooth[0]]), res_smooth, np.array([res_smooth[-1]])))
isPeak  = ( isPeak ) & (res_smooth_d[1:-1] >= res_smooth_d[0:-2]) & ( res_smooth_d[1:-1] > res_smooth_d[2:])

isPeak_back = isPeak
#itemindex = numpy.where(array==item)[0]
peak_locs = np.where(isPeak==True)[0]
#peak_locs = find(isPeak);

for i in xrange(1,len(peak_locs)+1):
	left = max(0,peak_locs[i-1]-LONGWIN+1)
	right = min(sig_len-1,peak_locs[i-1]+LONGWIN-1)
	ind = np.argmax(res_smooth[left:right+1])
	if((ind +left) != peak_locs[i-1]):
		#here there is a confusion as to LHS is ind+left or ind+left-1
		isPeak[peak_locs[i-1]] = False

peaks = np.where(isPeak==True)[0]
peaks_ecg = peaks


peaks_ecg_val = res_smooth[peaks_ecg]
p=np.zeros((len(peaks_ecg)-2,3))
p[:,0] = peaks_ecg_val[0:-2]
p[:,1] = peaks_ecg_val[1:-1]
p[:,2] = peaks_ecg_val[2:]

mask = (abs(p[:,0] - p[:,2]) <= 0.2*p[:,0]) & (p[:,1] <= 0.5*p[:,0])
ind = np.where(mask == False )[0]
peaks_ecg1 = np.concatenate((np.array([peaks_ecg[0]]), peaks_ecg[ind+1], np.array([peaks_ecg[-1]])))

peaks_ecg = peaks_ecg1
peaks_ecg_val = res_smooth[peaks_ecg]
p=np.zeros((len(peaks_ecg)-3,4))
p[:,0] = peaks_ecg_val[0:-3]
p[:,1] = peaks_ecg_val[1:-2]
p[:,2] = peaks_ecg_val[2:-1]
p[:,3] = peaks_ecg_val[3:]

mask = (abs(p[:,0]-p[:,3]) <= 0.2*p[:,0]) & (abs(p[:,1]-p[:,2]) <= 0.2*p[:,0]) & (p[:,1] <= 0.5*p[:,0])
ind = np.where(mask == True)[0]
peaks_ecg[ind+1]=0
peaks_ecg[ind+2]=0
peaks_ecg = peaks_ecg[np.where(peaks_ecg>0)[0]]
peaks_ecg1 = peaks_ecg
#peaks_ecg2 = peaks_ecg


isPeak= np.zeros(len(res_smooth))
isPeak[peaks_ecg] = 1
final = isPeak * res_smooth
LONGWIN = int(round(max(30,0.3*np.median(np.diff(peaks_ecg)))))


for i in xrange(1,len(peaks_ecg)+1):
	left = max(0,peaks_ecg[i-1]-LONGWIN+1)
	right = min(sig_len-1,peaks_ecg[i-1]+LONGWIN-1)
	ind = np.argmax(final[left:right+1])
	if((ind +left) != peaks_ecg[i-1]):
		isPeak[peaks_ecg[i-1]] = 0

peaks = np.where(isPeak == 1)[0]
final = isPeak * res_smooth

temp_sort = np.sort(res_smooth[peaks])
temp_sort_ecg = temp_sort
temp = max(0.1, 0.6*np.mean(temp_sort[0:int(round(0.2*len(temp_sort)))]))
peaks_ecg1 = np.where(final > temp)[0]
peaks_ecg = peaks_ecg1

np.savetxt('peaks_ecg1.txt',peaks_ecg1)
#############################
#new heuristic to remove false peaks
LONGWIN = 300
isPeak= np.zeros(len(res_smooth))
isPeak[peaks_ecg] = 1
final = isPeak * res_smooth


for i in xrange(1,len(res_smooth)-LONGWIN+1):
	left = i-1
	right = left+LONGWIN
	segment = np.copy(final[left:right+1])
	ind1 = np.argmax(segment)
	val1 = segment[ind1]
	segment[ind1]=0
	ind2 = np.argmax(segment)
	val2 = segment[ind2]
	val = (val1+val2)/2.0
	segment = np.copy(final[left:right+1])
	indieces = np.where(segment < 0.5*val)[0]
	indieces_1 = np.where((indieces > min((ind1, ind2))) & (indieces < max((ind1, ind2))))[0]
	segment[indieces[indieces_1]] = 0
	final[left:right+1] = segment

peaks_ecg1 = np.where(final > 0)[0]
peaks_ecg = peaks_ecg1

#######################




#peaks_ecg_val = res_smooth(peaks_ecg)
peaks = np.round(((peaks_ecg1+1)/250.0)*sampling_freq).astype(int)
#disp(size(peaks))
print "saving ecg peaks"
np.savetxt('peaks_ecg.txt',peaks)
#temp_sort = sort(diff(peaks_ecg1))
#LONGWIN = round(mean(temp_sort(1:round(0.2*length(temp_sort))))/2)



#post processing for BP
#bp_index = []
if(result_bp.size != 0):
	LONGWIN = 30
	#RES_TH = median(result_all)
	CNT_TH = 5

	sig_len = len(result_bp)

	#res_smooth = np.convolve(result_bp,np.ones(SHORTWIN+1)/(SHORTWIN+1),mode='same')
	res_smooth = result_bp
	isPeak = res_smooth > np.mean(res_smooth)
	#isPeak = isPeak.astype(int)

	#isPeak = res_smooth > 0.3
	res_smooth_d = np.concatenate((np.array([res_smooth[0]]), res_smooth, np.array([res_smooth[-1]])))
	isPeak  = ( isPeak ) & (res_smooth_d[1:-1] >= res_smooth_d[0:-2]) & ( res_smooth_d[1:-1] > res_smooth_d[2:])

	isPeak_back = isPeak
	#itemindex = np.where(array==item)
	peak_locs = np.where(isPeak==True)[0]
	#peak_locs = find(isPeak)

	for i in xrange(1,len(peak_locs)+1):
		left = max(0,peak_locs[i-1]-LONGWIN+1)
		right = min(sig_len-1,peak_locs[i-1]+LONGWIN-1)
		ind = np.argmax(res_smooth[left:right+1])
		if((ind +left) != peak_locs[i-1]):
			#here there is a confusion as to LHS is ind+left or ind+left-1
			isPeak[peak_locs[i-1]] = False

	peaks = np.where(isPeak==True)[0]
	peaks_bp = peaks


	peaks_bp_val = res_smooth[peaks_bp]
	p=np.zeros((len(peaks_bp)-2,3))
	p[:,0] = peaks_bp_val[0:-2]
	p[:,1] = peaks_bp_val[1:-1]
	p[:,2] = peaks_bp_val[2:]

	mask = (abs(p[:,0] - p[:,2]) <= 0.2*p[:,0]) & (p[:,1] <= 0.5*p[:,0])
	ind = np.where(mask == False )[0]
	peaks_bp1 = np.concatenate((np.array([peaks_bp[0]]),peaks_bp[ind+1],np.array([peaks_bp[-1]])))

	peaks_bp = peaks_bp1
	peaks_bp_val = res_smooth[peaks_bp]
	p=np.zeros((len(peaks_bp)-3,4))
	p[:,0] = peaks_bp_val[0:-3]
	p[:,1] = peaks_bp_val[1:-2]
	p[:,2] = peaks_bp_val[2:-1]
	p[:,3] = peaks_bp_val[3:]

	mask = (abs(p[:,0]-p[:,3]) <= 0.2*p[:,0]) & (abs(p[:,1]-p[:,2]) <= 0.2*p[:,0]) & (p[:,1] <= 0.5*p[:,0])
	ind = np.where(mask == True)[0]
	peaks_bp[ind+1]=0
	peaks_bp[ind+2]=0
	peaks_bp = peaks_bp[np.where(peaks_bp>0)[0]]
	peaks_bp1 = peaks_bp
	#peaks_bp2 = peaks_bp
	
	
	
	isPeak= np.zeros(len(res_smooth))
	isPeak[peaks_bp] = 1
	final = isPeak * res_smooth
	LONGWIN = int(round(max(30,0.3*np.median(np.diff(peaks_bp)))))
	
	
	for i in xrange(1,len(peaks_bp)+1):
		left = max(0,peaks_bp[i-1]-LONGWIN+1)
		right = min(sig_len-1,peaks_bp[i-1]+LONGWIN-1)
		ind = np.argmax(final[left:right+1])
		if((ind +left) != peaks_bp[i-1]):
			isPeak[peaks_bp[i-1]] = 0

	peaks = np.where(isPeak == 1)[0]
	final = isPeak * res_smooth

	temp_sort = np.sort(res_smooth[peaks])
	temp_sort_bp = temp_sort
	temp = max(0.1, 0.6*np.mean(temp_sort[0:int(round(0.2*len(temp_sort)))]))
	peaks_bp1 = np.where(final > temp)[0]


	#############################
	#new heuristic to remove false peaks
	LONGWIN = 300
	isPeak= np.zeros(len(res_smooth))
	isPeak[peaks_bp1] = 1
	final = isPeak * res_smooth


	for i in xrange(1,len(res_smooth)-LONGWIN+1):
		left = i-1
		right = left+LONGWIN
		segment = np.copy(final[left:right+1])
		ind1 = np.argmax(segment)
		val1 = segment[ind1]
		segment[ind1]=0
		ind2 = np.argmax(segment)
		val2 = segment[ind2]
		val = (val1+val2)/2.0
		segment = np.copy(final[left:right+1])
		indieces = np.where(segment < 0.5*val)[0]
		indieces_1 = np.where((indieces > min((ind1, ind2))) & (indieces < max((ind1, ind2))))[0]
		segment[indieces[indieces_1]] = 0
		final[left:right+1] = segment

	peaks_bp1 = np.where(final > 0)[0]
	peaks_bp = peaks_bp1

	

	
	print "peaks_bp1 ",peaks_bp1[0:10]
	np.savetxt('peaks_ecg1.txt',peaks_ecg1.astype(int))
	np.savetxt('peaks_bp1.txt',peaks_bp1.astype(int))
	var_bp = np.var(np.diff(peaks_bp1))
	var_ecg = np.var(np.diff(peaks_ecg1))
	done = 0
	if((abs(len(peaks_ecg)-len(peaks_bp)) < 0.8*len(peaks_ecg1)) and (abs(len(peaks_ecg)-len(peaks_bp)) < 0.8*len(peaks_bp1))):
		if(var_bp > 2*var_ecg):
			peaks = peaks_ecg1;
			done = 1
		elif(var_ecg > 2*var_bp):
			peaks = peaks_bp1;
			done = 1
	if(done == 0):
		interval = 5000
		#final = zeros(1,max(numel(peaks_bp1), numel(peaks_ecg1)))
		final = np.zeros((1,len(testData)))[0]
		i = 0
		shp = np.shape(testData)
		print "shape is",shp
		R = shp[0]
		C = shp[1]
		last = 0
		while(i<R):
			left = i
			right = min(i+interval-1,R-1)
			ecg_loc = np.concatenate((np.array([last]),peaks_ecg1[np.where((peaks_ecg1>left) & (peaks_ecg1<right))[0]]))
			bp_loc = np.concatenate((np.array([last]),peaks_bp1[np.where((peaks_bp1>left) & (peaks_bp1<right))[0]]))
			if ((ecg_loc.size == 0) or (np.diff(ecg_loc).size == 0)):
				final[bp_loc[1:]]=1
				last = bp_loc[-1]
			elif ((bp_loc.size == 0) or (np.diff(bp_loc).size == 0)):
				final[ecg_loc[1:]]=1
				last = ecg_loc[-1]
			elif 0.6*len(ecg_loc) >= len(bp_loc):
				final[ecg_loc[1:]]=1
				last = ecg_loc[-1]
			elif 0.6*len(bp_loc) >= len(ecg_loc):
				final[bp_loc[1:]]=1
				last = bp_loc[-1]
			else:
				if np.var(np.diff(ecg_loc))<=np.var(np.diff(bp_loc)): 
					final[ecg_loc[1:]]=1
					last = ecg_loc[-1]
				else:
					final[bp_loc[1:]]=1
					last = bp_loc[-1]
			i=i+interval
			result = np.where(final==1)[0]

		mask = np.where(np.diff(result)<=60)[0]
		print "mask is ", mask
		if(mask.size != 0):
			result[mask] = []
		peaks  = result
	print peaks[0:10]
	#peaks = peaks_bp1
	peaks = np.round(((peaks+1)/250.0)*sampling_freq).astype(int)

np.savetxt('peaks.txt',peaks)
