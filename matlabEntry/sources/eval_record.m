function [peaks] =  eval_record(file)
%file = 'mitbih_dataset/111';
%file  = 'mgh215.data';
%testData = load(strcat('../data/',num2str(file),'.txt'));

%load('../traindetect/train_params_ecg_nohidneuron.mat');
[isloaded,config] = wfdbloadlib;
if(config.inOctave)
crash_dumps_octave_core(0);
end
fprintf('isloaded = %d\n',isloaded);
annName = 'qrs'; % This annotator name is required by the Challenge rules.
siginfo = wfdbdesc(file); 
description = cell(1,numel(siginfo));
for k=1:numel(siginfo)
    description{k} = siginfo(k).Description;
end

%description = squeeze(struct2cell(siginfo));
%description = description(5,:);
% For more information, run 'help wfdbdesc'.
% Using the helper function defined below, find the ECG and blood pressure
% signals available in the record.
%ecg_index = get_index(description,'MLII');
%ecg_index = unique([1 get_index(description,'ECG')]);

%=== get blood pressure signal
bp_ind1 = cellfun(@any, regexpi(description,'bp','once'));
bp_ind2 = cellfun(@any, regexpi(description,'art','once'));
bp_ind3 = cellfun(@any, regexpi(description,'abp','once'));
bp_ind4 = cellfun(@any, regexpi(description,'ppg','once'));
bp_ind5 = cellfun(@any, regexpi(description,'pleth','once'));
bp_ind6 = cellfun(@any, regexpi(description,'pressure','once'));


if any(bp_ind1)
    bp_index = find(bp_ind1==1,1);
elseif any(bp_ind2)
    bp_index = find(bp_ind2==1,1);
elseif any(bp_ind3)
    bp_index = find(bp_ind3==1,1);
elseif any(bp_ind4)
    bp_index = find(bp_ind4==1,1);
elseif any(bp_ind5)
    bp_index = find(bp_ind5==1,1);
elseif any(bp_ind6)
    bp_index = find(bp_ind6==1,1);
else
    bp_index = [];
end


%bp_index = unique([get_index(description,'BP') get_index(description,'ABP') get_index(description,'ART') ...
%    get_index(description,'art') get_index(description,'PPG') get_index(description,'PLETH')]);

if (~isempty(bp_index))
    bp_index = bp_index(1);
end

%if (~isempty(ecg_index))
%    ecg_index = ecg_index(1);
%else
%	ecg_index = 1;
%end
ecg_index = 1;
%bp_index = [];

[~,testData,sampling_freq] = rdsamp(file,[ecg_index bp_index],[],[],[1],1);  
%testData = load(file);
%testData = testData(1:150000,:);
%sampling_freq = 360;
testData = resample(testData,250,sampling_freq);
if (~isempty(bp_index))
    bp_index = 2;
end

if (~isempty(ecg_index))
    ecg_index = 1;
end

if(~isempty(ecg_index) && ~(std(testData(:,ecg_index)) == 0))
    WINLEN_ecg = 50;
    load('traindetect/train_params_ecg.mat');
    %load('train_params_ecg.mat');
    %load('sources/traindetect/train_params_ecg_nohidneuron.mat');
    testData(:,ecg_index) = (testData(:,ecg_index) - mean(testData(:,ecg_index)))/std(testData(:,ecg_index));
    S_ecg = [zeros(WINLEN_ecg,1); testData(:,ecg_index) ;zeros(WINLEN_ecg,1)];
    [W_ecg] = findDetails(S_ecg);
    [m,n] = size(S_ecg);
    
    sig_len = m-WINLEN_ecg*2;
    
    %Setting up input features
    ip = zeros(sig_len,WINLEN_ecg*4+2);
    
    for j = 1:2*WINLEN_ecg+1
         ip(:,j) = W_ecg(4,j:j+sig_len-1);
         ip(:,j+2*WINLEN_ecg+1) = W_ecg(3,j:j+sig_len-1);
    end
    
    [m,n]= size(ip);
    % Normalization of the data
     inD4 = (ip(:,1:(n/2))-mean4_ecg)./var4_ecg;
     inD3 = (ip(:,n/2+1:n)-mean3_ecg)./var3_ecg;
     ip =  [inD4 inD3];
     ip = normalize(ip, mu_ecg, sigma_ecg);
     nn_ecg.testing = 1;
     nn = nnff(nn_ecg, ip, zeros(size(ip,1),2) );
     P = nn.a{nn.n};
    result_ecg = P(:,1);
    clear W_ecg;
    clear ip;
    clear nn_ecg; 
    clear nn; 
else
    result_ecg = [];
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%    NOW BP DATA 
if(~isempty(bp_index) && ~(std(testData(:,bp_index))==0))
    WINLEN_bp = 100; 
    load('traindetect/train_params_bp.mat');
    %load('train_params_bp.mat');
testData(:,bp_index) = (testData(:,bp_index) - mean(testData(:,bp_index)))/std(testData(:,bp_index));
    %S_bp = [zeros(WINLEN_bp,1); testData(:,3) ;zeros(WINLEN_bp,1)];
    S_bp = [zeros(WINLEN_bp,1); testData(:,bp_index) ;zeros(WINLEN_bp,1)];
    [W_bp] = findDetails(S_bp);
    [m,n] = size(S_bp);
    
    sig_len = m-WINLEN_bp*2;
    
    %Setting up input features
    ip = zeros(sig_len,WINLEN_bp*2+2);
    
    for j = 1:WINLEN_bp+1
             ip(:,j) = W_bp(4,2*j-1:2*j+sig_len-2);
         ip(:,j+WINLEN_bp+1) = W_bp(3,2*j-1:2*j+sig_len-2);
    end
    
    [m,n]= size(ip);
    
    
    % Normalization of the data
     inD4 = (ip(:,1:(n/2))-mean4_bp)./var4_bp;
     inD3 = (ip(:,n/2+1:n)-mean3_bp)./var3_bp;
     ip =  [inD4 inD3];
         ip = normalize(ip, mu_bp, sigma_bp);
     nn_bp.testing = 1;
     nn = nnff(nn_bp, ip, zeros(size(ip,1),2) );
     P = nn.a{nn.n};
    
    result_bp = P(:,1);
    clear W_bp;
    clear ip;
    clear nn_ecg; 
    clear nn; 
else
    result_bp = [];
end;

%if((~isempty(result_ecg)) && (~isempty(result_bp)))
 %   result_all = (result_ecg + result_bp)/2 ;
%elseif (~isempty(result_ecg))
 %   result_all = result_ecg ;
%elseif (~isempty(result_bp))
 %   result_all = result_bp ;
%else 
 %   result_all = zeros(sig_len,1);
%end
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SHORTWIN = 10;
LONGWIN = 60;
%RES_TH = median(result_all);
CNT_TH = 5;


sig_len = length(result_ecg);

res_smooth = conv(result_ecg,ones(SHORTWIN+1,1)/(SHORTWIN+1),'same');
isPeak = res_smooth > mean(res_smooth);

%isPeak = res_smooth > 0.3;
res_smooth_d = [res_smooth(1); res_smooth; res_smooth(end)];
isPeak  = ( isPeak ) & (res_smooth_d(2:end-1) >= res_smooth_d(1:end-2)) & ( res_smooth_d(2:end-1) > res_smooth_d(3:end)) ;

isPeak_back = isPeak;

peak_locs = find(isPeak);

for i = 1:length(peak_locs)
   left = max(1,peak_locs(i)-LONGWIN+1);
   right = min(sig_len,peak_locs(i)+LONGWIN-1);
   [~,ind] = max(res_smooth(left:right));
   if((ind +left-1) ~= peak_locs(i))
        isPeak(peak_locs(i)) = 0;
   end
end

peaks = find(isPeak);
peaks_ecg = peaks;



%peaks = find(isPeak);
%peaks_ecg = peaks;
%peaks_ecg3 = peaks;





peaks_ecg_val = res_smooth(peaks_ecg);
p=zeros(length(peaks_ecg)-2,3);
p(:,1) = peaks_ecg_val(1:end-2);
p(:,2) = peaks_ecg_val(2:end-1);
p(:,3) = peaks_ecg_val(3:end);
    
mask = (abs(p(:,1) - p(:,3))<=0.2*p(:,1)) & (p(:,2)<=0.5*p(:,1));
ind = find(mask == 0 );
peaks_ecg1 = [peaks_ecg(1); peaks_ecg(ind+1); peaks_ecg(end)];

peaks_ecg = peaks_ecg1;
peaks_ecg_val = res_smooth(peaks_ecg);
p=zeros(length(peaks_ecg)-3,4);
p(:,1) = peaks_ecg_val(1:end-3);
p(:,2) = peaks_ecg_val(2:end-2);
p(:,3) = peaks_ecg_val(3:end-1);
p(:,4) = peaks_ecg_val(4:end);

mask = (abs(p(:,1)-p(:,4))<=0.2*p(:,1)) & (abs(p(:,2)-p(:,3))<=0.2*p(:,1)) & (p(:,2)<=0.5*p(:,1));
ind = find(mask==1);
peaks_ecg(ind+1)=0;
peaks_ecg(ind+2)=0;
peaks_ecg = peaks_ecg(find(peaks_ecg>0));
peaks_ecg1 = peaks_ecg;
%peaks_ecg2 = peaks_ecg;



isPeak= zeros(size(res_smooth));
isPeak(peaks_ecg) = 1;
final = isPeak.* res_smooth;
LONGWIN = round(max(60,0.6*median(diff(peaks_ecg))))


for i = 1:length(peaks_ecg)
   left = max(1,peaks_ecg(i)-LONGWIN+1);
   right = min(sig_len,peaks_ecg(i)+LONGWIN-1);
   [~,ind] = max(final(left:right));
   if((ind +left-1) ~= peaks_ecg(i))
        isPeak(peaks_ecg(i)) = 0;
   end
end

final = isPeak.* res_smooth;

temp_sort = sort(res_smooth(peaks));
temp_sort_ecg = temp_sort;
temp = max(0.1,0.6*mean(temp_sort(1:round(0.2*length(temp_sort)))))
peaks_ecg1 = find(final > temp);
peaks_ecg = peaks_ecg1;	
%peaks_ecg_val = res_smooth(peaks_ecg);
peaks = round(peaks_ecg1'/250*sampling_freq);
%disp(size(peaks));

%temp_sort = sort(diff(peaks_ecg1));
%LONGWIN = round(mean(temp_sort(1:round(0.2*length(temp_sort))))/2);

if ~isempty(bp_index)
LONGWIN = 60;
	sig_len = length(result_bp)
	
	res_smooth = conv(result_bp,ones(SHORTWIN+1,1)/(SHORTWIN+1),'same');
	isPeak = res_smooth > mean(res_smooth);
	res_smooth_d = [res_smooth(1); res_smooth; res_smooth(end)];
	isPeak  = ( isPeak ) & (res_smooth_d(2:end-1) >= res_smooth_d(1:end-2)) & ( res_smooth_d(2:end-1) > res_smooth_d(3:end)) ;
	
	isPeak_back = isPeak;
	
	peak_locs = find(isPeak);
	
	for i = 1:length(peak_locs)
	   left = max(1,peak_locs(i)-LONGWIN+1);
	   right = min(sig_len,peak_locs(i)+LONGWIN-1);
	   [~,ind] = max(res_smooth(left:right));
	   if((ind +left-1) ~= peak_locs(i))
	        isPeak(peak_locs(i)) = 0;
	   end
	end
	
	peaks = find(isPeak);
	peaks_bp= peaks;
	peaks_bp1 = peaks;







peaks_bp_val = res_smooth(peaks_bp);
p=zeros(length(peaks_bp)-2,3);
p(:,1) = peaks_bp_val(1:end-2);
p(:,2) = peaks_bp_val(2:end-1);
p(:,3) = peaks_bp_val(3:end);
    
mask = (abs(p(:,1) - p(:,3))<=0.2*p(:,1)) & (p(:,2)<=0.5*p(:,1));
ind = find(mask == 0 );
peaks_bp1 = [peaks_bp(1); peaks_bp(ind+1); peaks_bp(end)];

peaks_bp = peaks_bp1;
peaks_bp_val = res_smooth(peaks_bp);
p=zeros(length(peaks_bp)-3,4);
p(:,1) = peaks_bp_val(1:end-3);
p(:,2) = peaks_bp_val(2:end-2);
p(:,3) = peaks_bp_val(3:end-1);
p(:,4) = peaks_bp_val(4:end);

mask = (abs(p(:,1)-p(:,4))<=0.2*p(:,1)) & (abs(p(:,2)-p(:,3))<=0.2*p(:,1)) & (p(:,2)<=0.5*p(:,1));
ind = find(mask==1);
peaks_bp(ind+1)=0;
peaks_bp(ind+2)=0;
peaks_bp = peaks_bp(find(peaks_bp>0));
peaks_bp1 = peaks_bp;
%peaks_bp2 = peaks_bp;









isPeak= zeros(size(res_smooth));
isPeak(peaks_bp) = 1;
final = isPeak.* res_smooth;
LONGWIN = round(max(60,0.6*median(diff(peaks_bp))))


for i = 1:length(peaks_bp)
   left = max(1,peaks_bp(i)-LONGWIN+1);
   right = min(sig_len,peaks_bp(i)+LONGWIN-1);
   [~,ind] = max(final(left:right));
   if((ind +left-1) ~= peaks_bp(i))
        isPeak(peaks_bp(i)) = 0;
   end
end

peaks = find(isPeak);
%peaks_bp = peaks;
%peaks_bp3 = peaks;


final = isPeak.* res_smooth;

temp_sort = sort(res_smooth(peaks));
temp = max(0.1,0.6*mean(temp_sort(1:round(0.2*length(temp_sort)))))
peaks_bp1 = find(final > temp);


interval = 5000;
%final = zeros(1,max(numel(peaks_bp1), numel(peaks_ecg1)));
final = zeros(1,length(testData));
i = 1;
[R,C] = size(testData);
last = 1;
while(i<=R)
	left = i;
	right = min(i+interval-1,R);
	ecg_loc = [last ;peaks_ecg1(peaks_ecg1>left & peaks_ecg1<right)];
	bp_loc = [last ;peaks_bp1(peaks_bp1>left & peaks_bp1<right)];
	if ( isempty(ecg_loc) || isempty(diff(ecg_loc)) )
		final(bp_loc(2:end))=1;
last = bp_loc(end);
	elseif (isempty(bp_loc) || isempty(diff(bp_loc)))
		final(ecg_loc(2:end))=1;
last = ecg_loc(end);
	elseif 0.6*length(ecg_loc) >= length(bp_loc)
                final(ecg_loc(2:end))=1;
last = ecg_loc(end);
        elseif 0.6*length(bp_loc) >= length(ecg_loc)
                final(bp_loc(2:end))=1;
last = bp_loc(end);
	else
            	
		if var(diff(ecg_loc))<=var(diff(bp_loc)) 
			final(ecg_loc(2:end))=1;
last = ecg_loc(end);
		else
			final(bp_loc(2:end))=1;
last = bp_loc(end);
		end
	end	
        i=i+interval;	
	result = find(final==1);
	
end

mask = find(diff(result)<=60);
result(mask+1) = [];
peaks  = result;






peaks = round(peaks'/250*sampling_freq);

end
%keyboard;
wrann(file,annName,peaks,[],[],[],[]);
%disp(['time for ',file,' : ', num2str(t2)])
end
   

function ind = get_index(description,pattern)
	M = length(description);
	tmp_ind = strfind(description,pattern);
	ind = [];
	for m = 1:M
	    if(~isempty(tmp_ind{m}))
	        ind(end+1) = m;
		break;
	    end
	end
end
