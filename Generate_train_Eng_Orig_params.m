function train=Generate_train_Eng_Orig_params(params,N)


epsilon=params(1);
hc=params(2);
lambda1=params(3);
rho1=params(4);
lambda2=params(5);
rho2=params(6);

% mu1=200;sigma1=50;hc=7.83;

%Engbert's model from 2011 paper plus Roy's inhibition between MS addition, works well.
visualize=0;

MSTIM=[];cnt=1;
%Parameters::
L=51;
% epsilon=.001;
% hc=7.83;
lambda=1;
% runs=12;


% mu_sigma(1,:) = [38,115];%[35 107]
% mu_sigma(2,:) = [27,144];%[32,130]
% mu_sigma(3,:) = [34,150];%[42 148]

MSTIM=[];cnt=1;
timesincelast=10000;

delayfac=0;
%tic

    %iternum=1000;

    
    %Initialize maps and variables:
    H=zeros(L);%100,L,L);
    startpos=[randi(L),randi(L)];
    i=startpos(1);
    j=startpos(2);
    i0=round((L-1)/2);
    j0=round((L-1)/2);
    train= false(1,N);
    MS_t_counter=0;
    
    
    potmap=zeros(L);
    for w=1:L
        for ww=1:L
            potmap(w,ww)=lambda*L*(((w-i0)/i0)^2+((ww-j0)/j0)^2);
        end
    end
    j_m1=0;i_m1=0;
%     itrack=zeros(1,N);
%     jtrack=zeros(1,N);
    track_counter=0;
    
    for iter=1:N
        
       if (iter==15000) & ((MS_t_counter==0)|(MS_t_counter>50))
           
           return
       end
        
        neighbours(1,1:2)=[i+1,j];
        neighbours(2,1:2)=[i-1,j];
        neighbours(3,1:2)=[i,j+1];
        neighbours(4,1:2)=[i,j-1];
        neighbours(neighbours<1)=L;
        neighbours(neighbours>L)=1;
        ap_t = ((1+lambda1*exp(-rho1*t_from_last_ms^2))^-1);
        potmap_t = potmap*ap_t;
        
        neighbourvals=[H(neighbours(1,1),neighbours(1,2)),H(neighbours(2,1),neighbours(2,2)),H(neighbours(3,1),neighbours(3,2)),H(neighbours(4,1),neighbours(4,2))];
        neighbourvals=neighbourvals+[potmap_t(neighbours(1,1),neighbours(1,2)),potmap_t(neighbours(2,1),neighbours(2,2)),potmap_t(neighbours(3,1),neighbours(3,2)),potmap_t(neighbours(4,1),neighbours(4,2))];
        mino=101; %this will store the minimum
        indexes=[]; %and it's index
        
        
        [mino, indexes]=min(neighbourvals);
        if sum(neighbourvals==mino)>1
            [~,inds]=find(min(neighbourvals)==neighbourvals);
            indexes=inds(randi(length(inds)));
        end
        
        
        hc_t = hc*(1 + (1/(1 + lambda2*exp(rho2*t_from_last_ms^2))));
        if (H(i,j)*ap_t + potmap_t(i,j))>=hc_t  % MS triggering mechanism
            timesincelast=0;
%            MS_t=[ MS_t iter]; %store the "time"
            MS_t_counter = MS_t_counter+1;
            train(iter) = 1;
           
            
            delayfac=randn*sigma1+mu1;
            %disp('MS!!')
            %Nextt lines help calcuate global min as a destination for the MS.
            K=H+potmap;
            mino=150;
            indice=[];
            [val, i_ind_vec] = min(K);
            [val, j_ind] = min(val);
            i_ind = i_ind_vec(j_ind);
            
            i = i_ind;
            j = j_ind;
            
            %        [minimal_value,[i,j]]=(min(k));
            %        temp=H(i,j)+1;
            %        H=H*(1-epsilon);
            %        H(i,j)=temp;
            %        imagesc(H)
            %        pause(0.01)
            %        continue;
        else %else (no MS ) go to the local minimum we calculated.
            timesincelast=timesincelast+1;
            i=neighbours(indexes,1);
            j=neighbours(indexes,2);
            delayfac=delayfac-1;
        end
        temp=H(i,j)+1; %store new value after stepping
        H=H*(1-epsilon); %update all Matrix
        H(i,j)=temp; %replace last stepped with the tmp.
        if visualize
            imagesc(H)
            hold on
            % plot([1,j_m1],[1,i_m1],'gx')
            plot([1,j],[1,i],'kx')
            pause(0.01);
        end
        %    i_m1=i;
        %    j_m1=j;
%         itrack=[itrack i];
%         jtrack=[jtrack j];
      
%         track_counter = track_counter + 1;
%         itrack(track_counter) = i;
%         jtrack(track_counter) = j;
%         HS(track_counter)=H(i,j);
        
    end
    %MSTIM{cnt}=MS_t;



