function [W] = findDetails(S2)
% This function finds the first five details of a signal(S2)
    H = [0.125 0.375 0.375 0.125];
    G = [-2 2];

    lambda = [1.5 1.12 1.03 1.01 1.00];

    j = 1;
    numOfStages = 4;

    W=zeros(numOfStages,size(S2,1));
   % S2 = S;
    while(j<=numOfStages)
       W(j,:) = conv(S2, G, 'same')/lambda(j);
       S2 = conv(S2, H, 'same');
       j = j+1;
    end
end
