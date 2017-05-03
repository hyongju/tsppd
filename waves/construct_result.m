function [ XX ] = construct_result( XX,v,n,k)
    capa = 0;
%     cont = 0;


    for i = 1:v-1
%        cont = cont+1;
       flag = 1;
       while flag
           [~,index] = max(XX(i,:));
           if 2<=index && index<=1+n, % potential pick-up
               if sum(find(XX(:,index)==1))==0 && capa+1<=k,
                   XX(i,:) = [zeros(1,index-1),1,zeros(1,v-index)];
                   XX(:,index) = [zeros(i-1,1);1;zeros(v-i,1)];
                   flag = 0;
                   capa = capa+1;
               else
                   XX(i,index)=0;
               end
           else          % potential drop-off
               if sum(find(XX(:,index)==1))==0 && sum(find(XX(:,index-n)==1))~=0,
                   XX(i,:) = [zeros(1,index-1),1,zeros(1,v-index)];
                   XX(:,index) = [zeros(i-1,1);1;zeros(v-i,1)];
                   flag = 0;
                   capa = capa-1;
               else
                   XX(i,index)=0;
               end
           end

       end
    end    

end

