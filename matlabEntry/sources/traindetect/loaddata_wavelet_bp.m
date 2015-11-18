%function loaddata_wavelet(from, to)
y=zeros(60*2000,1);
WINLEN = 100;
count = 1;
%x = zeros(60*2000,4*WINLEN+2);
x = zeros(60*2000,2*WINLEN+2);
%x = zeros(60*2000,2*WINLEN+2);
%for i=from:from
for i=101:101
    testData = load(strcat(num2str(i),'.txt'));
    
S = testData(:,3);
%S = result;
H = [0.125 0.375 0.375 0.125];
G = [-2 2];

lambda = [1.5 1.12 1.03 1.01 1.00];
j=1;
J=4;
W= zeros(J,size(S,1));
S2=S;

while(j<=J)
    W(j, :) = conv(S2,G, 'same')/lambda(j);
    S2 = conv(S2,H,'same');
    j=j+1;
end

D1= W(1,:);
D2= W(2,:);
D3= W(3,:);
D4= W(4,:);

S2=S2';

    %testDataAtr = load(strcat(num2str(i),' (2).txt'));
testDataAtr(:,2) = load(strcat(num2str(i),'_1.txt'));
    [m,n] = size(testData);
    [m_atr,n_atr] = size(testDataAtr);
    
    for j=1:m_atr
    k=testDataAtr(j,2);
    if(k <=m-WINLEN  && k>WINLEN)
        %tempData = (testData(k-WINLEN:k+WINLEN,3))';
        %tempData = [D4(k-WINLEN:1:k+WINLEN) D3(k-WINLEN:1:k+WINLEN)] ;
        tempData = [D4(k-WINLEN:2:k+WINLEN) D3(k-WINLEN:2:k+WINLEN)] ;
        x(count,:) =  tempData;
        y(count) = 1;
        count = count+1;
    end
    end
    
    
     

r_val = zeros(1,1000);
%[m_x, n_x] = size(x);
t=1;
% while(t<=200)
%     peak = 0;
%     r= round(rand * (m-1) + 1);
%     if(r <m-WINLEN  && r>WINLEN)
%        for k=1:m_atr
%        %if(r == testDataAtr(k,2))
%         if(abs(r-testDataAtr(k,2))<10)
%             peak = 1;
%            break;
%        end
%        end
%       %peak = isempty(find(testDataAtr(:,2)==r));
%         if(peak ==0)
%             %tempData = (testData(r-WINLEN:r+WINLEN,3))';
%             tempData = [D4(r-WINLEN:1:r+WINLEN) D3(r-WINLEN:1:r+WINLEN)] ;
%              %tempData = [D4(r-WINLEN:2:r+WINLEN) D3(r-WINLEN:2:r+WINLEN)] ;
%              x(count,:) =  tempData;
%              y(count) = 0;
%             r_val(count) = r;
%             count = count+1;
%             t=t+1;
%             
%         end
%     end
% end

% while(t<=10)
%     peak = 0;
%     r1= round(rand * (m_atr-1) + 1);
%     r= round((rand - 0.5)*20);
%     if(r==0)
%         r = r1+1;
%     else
%         r=r1+r;
%     end
%     if(r <m-WINLEN  && r>WINLEN)
% %        for k=1:m_atr
% %        %if(r == testDataAtr(k,2))
% %         if(abs(r-testDataAtr(k,2))>10)
% %             peak = 1;
% %            break;
% %        end
% %        end
%       %peak = isempty(find(testDataAtr(:,2)==r));
% %         if(peak ==0)
%             %tempData = (testData(r-WINLEN:r+WINLEN,3))';
%            tempData = [D4(r-WINLEN:1:r+WINLEN) D3(r-WINLEN:1:r+WINLEN)] ;
% %              tempData = [D4(r-WINLEN:2:r+WINLEN) D3(r-WINLEN:2:r+WINLEN)] ;
%              x(count,:) =  tempData;
%              y(count) = 0;
%             r_val(count) = r;
%             count = count+1;
%             t=t+1;
%             
%         %end
%     end
% end
% 
% while(t<=210)
%     peak = 0;
%     %r= round(rand * (m-1) + 1);
%      r1= round(rand * (m_atr-1) + 1);
%     r= round((rand - 0.5)*30); 
%     if(r<0)
%         r=r-10+r1;
%     else
%         r=r+10+r1;
%     end
%     if(r <m-WINLEN  && r>WINLEN)
% %        for k=1:m_atr
% %        %if(r == testDataAtr(k,2))
% %         if(~((r-testDataAtr(k,2))<-10 && (r-testDataAtr(k,2))>-25 ))
% %             peak = 1;
% %            break;
% %        end
% %        end
%       %peak = isempty(find(testDataAtr(:,2)==r));
%         %if(peak ==0)
%             %tempData = (testData(r-WINLEN:r+WINLEN,3))';
%             tempData = [D4(r-WINLEN:1:r+WINLEN) D3(r-WINLEN:1:r+WINLEN)] ;
%             %tempData = [D4(r-WINLEN:2:r+WINLEN) D3(r-WINLEN:2:r+WINLEN)] ;
%              x(count,:) =  tempData;
%              y(count) = 0;
%             r_val(count) = r;
%             count = count+1;
%             t=t+1;
%             
%        % end
%     end
% end
% 
% while(t<=1000)
%     peak = 0;
%     r= round(rand * (m-1) + 1);
%     if(r <m-WINLEN  && r>WINLEN)
%        for k=1:m_atr
%       % if(r == testDataAtr(k,2))
%         %if((r-testDataAtr(k,2))<-10 && (r-testDataAtr(k,2))>-30 )
%            if(abs(r-testDataAtr(k,2))<30)
%                peak = 1;
%            break;
%            end
%        end
%       %peak = isempty(find(testDataAtr(:,2)==r));
%         if(peak ==0)
%             %tempData = (testData(r-WINLEN:r+WINLEN,3))';
%            % tempData = [D4(r-WINLEN:1:r+WINLEN) D3(r-WINLEN:1:r+WINLEN)] ;
%             tempData = [D4(r-WINLEN:2:r+WINLEN) D3(r-WINLEN:2:r+WINLEN)] ;
%              x(count,:) =  tempData;
%              y(count) = 0;
%             r_val(count) = r;
%             count = count+1;
%             t=t+1;
%             
%         end
%     end
% end

while(t<=1000)
    peak = 0;
    r= round(rand * (m-1) + 1);
    if(r <m-WINLEN  && r>WINLEN)
       for k=1:m_atr
      if(r == testDataAtr(k,2))
        %if((r-testDataAtr(k,2))<-10 && (r-testDataAtr(k,2))>-30 )
         % if(abs(r-testDataAtr(k,2))<30)
               peak = 1;
           break;
      end
       end
      %peak = isempty(find(testDataAtr(:,2)==r));
        if(peak ==0)
            %tempData = (testData(r-WINLEN:r+WINLEN,3))';
            %tempData = [D4(r-WINLEN:1:r+WINLEN) D3(r-WINLEN:1:r+WINLEN)] ;
             tempData = [D4(r-WINLEN:2:r+WINLEN) D3(r-WINLEN:2:r+WINLEN)] ;
             x(count,:) =  tempData;
             y(count) = 0;
            r_val(count) = r;
            count = count+1;
            t=t+1;
            
        end
    end
end



i=i
end
 
x=x(1:count-1,:);
y=y(1:count-1,1);

