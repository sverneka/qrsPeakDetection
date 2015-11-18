%function [er,bad,nn] = neural_network(epochs,hidden_neuron,lambda,regularize,dropout)
epochs = 50;
hidden_neuron = 20;
lambda = 0.0001;
regularize = 1;
dropout =1;


addpath(genpath('/home/hwr/tagme/codes/DeepLearnToolbox-master'));

BPdata = load('train_test_bp.mat');
train_x = BPdata.train_bp(:,1:end-1);
train_y = BPdata.train_bp(:,end);
train_y = [train_y,1-train_y];

test_x  = BPdata.test_bp(:,1:end-1);
test_y  = BPdata.test_bp(:,end);
test_y=[test_y,1-test_y];

mean3_bp = BPdata.mean3_bp;
mean4_bp = BPdata.mean4_bp;

var3_bp = BPdata.var3_bp;
var4_bp = BPdata.var4_bp;
%BPdata = load('../train_data_bp.txt');
%BPdata = BPdata(1:150000,:);
%
%cvrun = 1;
%
%ind1 = 145001;
%ind2 = 150000;
%
%test_x = [BPdata(ind1:ind2,1:202)];
%test_y = [BPdata(ind1:ind2,203)];
%
%train_x = BPdata(setdiff(1:size(BPdata,1),[ind1:ind2]),1:202);
%train_y = BPdata(setdiff(1:size(BPdata,1),[ind1:ind2]),203);
%
%train_y = [train_y,1-train_y];
%test_y=[test_y,1-test_y];


clear BPdata;


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
%nn = nnsetup([202 20 2]);
nn=nnsetup([size(train_x,2), 2]);
%opts.numepochs =  20;   %  Number of full sweeps through data
opts.numepochs =epochs;
opts.batchsize = 500;  %  Take a mean gradient step over this many samples
[nn, L] = nntrain(nn, train_x, train_y, opts);


[er, bad] = nntest(nn, test_x, test_y);

nn_bp = rmfield(nn,'a');

mu_bp = mu;
sigma_bp = sigma;


disp(sprintf('error = %g',er));

save('-mat', 'train_params_bp.mat','nn_bp','mu_bp','sigma_bp','mean3_bp','mean4_bp','var3_bp','var4_bp');
%=============================================
