NUM_TR = 150000;
x = [load('../train60_bp') ;  load('../test30_bp')];

input = x(1:NUM_TR,1:end-1);
[mean3_bp, mean4_bp,var3_bp,var4_bp, train_bp] =  shuffle_train_data_wavelet(input);
train_bp = [train_bp,x(1:NUM_TR,end)];

input = x(NUM_TR+1:end,1:end-1);
[m,n]= size(input);
inD4 = (input(:,1:(n/2))-mean4_bp)./var4_bp;
inD3 = (input(:,1+(n/2):n)-mean3_bp)./var3_bp;
test_bp = [inD4,inD3,x(NUM_TR+1:end,end)];

x = [load('../train60_ecg') ;  load('../test30_ecg')];
input = x(1:NUM_TR,1:end-1);
[mean3_ecg, mean4_ecg,var3_ecg,var4_ecg, train_ecg] =  shuffle_train_data_wavelet(input);
train_ecg = [train_ecg,x(1:NUM_TR,end)];

input = x(NUM_TR+1:end,1:end-1);
[m,n]= size(input);
inD4 = (input(:,1:(n/2))-mean4_ecg)./var4_ecg;
inD3 = (input(:,1+(n/2):n)-mean3_ecg)./var3_ecg;
test_ecg = [inD4,inD3,x(NUM_TR+1:end,end)];

save ('-mat','train_test_bp.mat','train_bp','test_bp','mean3_bp','mean4_bp','var3_bp','var4_bp');
save ('-mat','train_test_ecg.mat','train_ecg','test_ecg','mean3_ecg','mean4_ecg','var3_ecg','var4_ecg');
