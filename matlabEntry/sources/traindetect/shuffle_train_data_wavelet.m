function [mean3, mean4,var3,var4, input] =  shuffle_train_data_wavelet(input)
% [m,n] = size(x);
% 
% input = x(2:m,:);


 [m,n]= size(input);
 inD4 = input(:,1:(n/2));
 inD3 = input(:,1+(n/2):n);
 mean4 = sum(sum(inD4))/(m*n/2);
 mean3 = sum(sum(inD3))/(m*n/2);
 var4 = sqrt(sum(sum((inD4-mean4).^2))/(m*n/2));
 var3 = sqrt(sum(sum((inD3-mean3).^2))/(m*n/2));
 inD4 = (inD4-mean4)./var4;
 inD3 = (inD3-mean3)./var3;
 input = [inD4 inD3];

