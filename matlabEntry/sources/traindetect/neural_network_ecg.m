epochs = 150;
epochs = 20;
hidden_neuron = 100;
lambda = 0.0001;
regularize = 1;
dropout =1;


addpath(genpath('/home/hwr/tagme/codes/DeepLearnToolbox-master'));

ECGdata = load('train_test_ecg.mat');
train_x = ECGdata.train_ecg(:,1:end-1);
train_y = ECGdata.train_ecg(:,end);
train_y = [train_y,1-train_y];

test_x  = ECGdata.test_ecg(:,1:end-1);
test_y  = ECGdata.test_ecg(:,end);
test_y=[test_y,1-test_y];

mean3_ecg = ECGdata.mean3_ecg;
mean4_ecg = ECGdata.mean4_ecg;

var3_ecg = ECGdata.var3_ecg;
var4_ecg = ECGdata.var4_ecg;
%ECGdata = load('../train_data_ecg.txt');
%ECGdata = ECGdata(1:150000,:);
%
%cvrun = 1;
%
%ind1 = 145001;
%ind2 = 150000;
%
%test_x = [ECGdata(ind1:ind2,1:202)];
%test_y = [ECGdata(ind1:ind2,203)];
%
%train_x = ECGdata(setdiff(1:size(ECGdata,1),[ind1:ind2]),1:202);
%train_y = ECGdata(setdiff(1:size(ECGdata,1),[ind1:ind2]),203);
%
%train_y = [train_y,1-train_y];
%test_y=[test_y,1-test_y];


clear ECGdata;


%randomize
I = randperm(size(train_x,1));
train_x = train_x(I,:);
train_y = train_y(I,:);

%Normalization
[train_x, mu, sigma] = zscore(train_x);
test_x = normalize(test_x, mu, sigma);

%nn.weightPenaltyL2 = 1e-4;
if(regularize == 1)
    nn.weightPenaltyL2 = lambda;
end

if(dropout == 1)
   nn.dropoutFraction = 0.5;   %  Dropout fraction 
end
nn.output = 'softmax';
nn=nnsetup([size(train_x,2),  2]);
%opts.numepochs =  20;   %  Number of full sweeps through data
opts.numepochs =epochs;
opts.batchsize = 500;  %  Take a mean gradient step over this many samples
[nn, L] = nntrain(nn, train_x, train_y, opts);


[er, bad] = nntest(nn, test_x, test_y);

nn_ecg = rmfield(nn,'a');
mu_ecg = mu;
sigma_ecg = sigma;


disp(sprintf('error = %g',er));

save('-mat', 'train_params_ecg_nohidneuron.mat','nn_ecg','mu_ecg','sigma_ecg','mean3_ecg','mean4_ecg','var3_ecg','var4_ecg');
%=============================================
