function test_tagme

addpath(genpath('/home/hwr/tagme/codes/DeepLearnToolbox-master'));

x1 = load('../../../features/Train/alldatahog00');
x2 = load('../../../features/Train/alldatahog01');
x3 = load('../../../features/Train/alldatahog02');
x4 = load('../../../features/Train/alldatahog03');

y1 = load('../../../features/Train/labels00');
y2 = load('../../../features/Train/labels01');
y3 = load('../../../features/Train/labels02');
y4 = load('../../../features/Train/labels03');


%train_x = double(train_x) / 255;
%test_x  = double(test_x)  / 255;
%train_y = double(train_y);
%test_y  = double(test_y);
cverr = 0;

for cvrun = 1:4

train_x = [x1;x2;x3];
train_y = [y1;y2;y3];

test_x = x4;
test_y = y4;

%%  ex1 train a 100 hidden unit RBM and visualize its weights
%rand('state',0)
%dbn.sizes = [100];
%opts.numepochs =  10;
%opts.batchsize = 100;
%opts.momentum  =   0;
%opts.alpha     =   1;
%dbn = dbnsetup(dbn, train_x, opts);
%dbn = dbntrain(dbn, train_x, opts);
%figure; visualize(dbn.rbm{1}.W');   %  Visualize the RBM weights

%%  ex2 train a 100-100 hidden unit DBN and use its weights to initialize a NN
rand('state',0)
%train dbn
dbn.sizes = [200];
opts.numepochs = 10;
opts.batchsize = 25;
opts.momentum  =   0;
opts.alpha     =   1;
dbn = dbnsetup(dbn, train_x, opts);
dbn = dbntrain(dbn, train_x, opts);

%unfold dbn to nn
nn = dbnunfoldtonn(dbn, 5);
nn.activation_function = 'sigm';

%train nn
opts.numepochs = 500;
opts.batchsize = 25;
nn = nntrain(nn, train_x, train_y, opts);
[er, bad] = nntest(nn, test_x, test_y);


% Update the error
cverr = cverr + er; 

% rotate for next cv 
% Not an efficient way but okay for the time being
tx = x1;
x1 = x2;
x2 = x3;
x3 = x4;
x4 = tx;

ty = y1;
y1 = y2;
y2 = y3;
y3 = y4;
y4 = ty;

disp(sprintf('error = %g',er));
end;

disp(sprintf('cv error = %g',cverr/4));
%assert(er < 0.10, 'Too big error');
