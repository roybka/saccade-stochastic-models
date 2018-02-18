function [ train ] = create_train_from_haz(haz ,nsacs,l)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
cutoff=1500;
trace=false(1,l);
        cur=1000;
        newhaz=[haz(1:cutoff) ones(1,10000)*mean(haz(cutoff-200:cutoff))];
        for i=1:nsacs
            sacdone=0;
            sacdone=find(rand(1,length(haz(1:cutoff))+10000)<newhaz,1,'first');
            if sacdone
                
                trace(cur+sacdone)=1;
                gap=sacdone;
                %gaps=[gaps cnt2];
                cur=cur+sacdone;
            else
                cur=cur+1500;
            end
            
            
            cnt2=0;
            %         while ~sacdone & cnt2<1500
            %             cnt2=cnt2+1;
            %             sacdone=rand<haz(cnt2);
            %             if sacdone
            %
            %                 trace(cur+cnt2)=1;
            %                 gap=cnt2;
            %                 %gaps=[gaps cnt2];
            %                 cur=cur+cnt2;
            %                 break
            %             end
            %         end
        end
        
  train=trace;      
end

