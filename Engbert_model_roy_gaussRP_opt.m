clear

%Engbert's model from 2011 paper plus Roy's inhibition between MS addition, works well.
visualize=0;

MSTIM=[];cnt=1;
%Parameters::
L=51;
epsilon=.001;
hc=7.9;
lambda=1;
runs=30;
delayfac=0;
timesincelast=10000;

mu_sigma(1,:) = [35,107];
mu_sigma(2,:) = [32,130];
mu_sigma(3,:) = [43,148];

cont_ind = 2;

sigma1=mu_sigma(cont_ind,1);
mu1=mu_sigma(cont_ind,2);

h = waitbar(0,'Running Engbert Amit...');
run_num=-1;
tic
for iternum=ones(1,runs)*180000
    %iternum=1000;
    run_num = run_num + 1;

    
    %Initialize maps and variables:
    H=zeros(L);%100,L,L);
    startpos=[randi(L),randi(L)];
    i=startpos(1);
    j=startpos(2);
    i0=round((L-1)/2);
    j0=round((L-1)/2);
    MS_t= zeros(1,180000);
    MS_t_counter=0;
    
    
    potmap=zeros(L);
    for w=1:L
        for ww=1:L
            potmap(w,ww)=lambda*L*(((w-i0)/i0)^2+((ww-j0)/j0)^2);
        end
    end
    j_m1=0;i_m1=0;
    itrack=zeros(1,180000);
    jtrack=zeros(1,180000);
    track_counter=0;
    
    for iter=1:iternum
        
        if(mod(iter,1000)==0)
            waitbar((run_num*iternum + iter)/(iternum*runs),h,['Running Engbert Amit %' num2str(round((run_num*iternum + iter)/(iternum*runs)*100))])
        end
        
        neighbours(1,1:2)=[i+1,j];
        neighbours(2,1:2)=[i-1,j];
        neighbours(3,1:2)=[i,j+1];
        neighbours(4,1:2)=[i,j-1];
        neighbours(neighbours<1)=L;
        neighbours(neighbours>L)=1;
        
        neighbourvals=[H(neighbours(1,1),neighbours(1,2)),H(neighbours(2,1),neighbours(2,2)),H(neighbours(3,1),neighbours(3,2)),H(neighbours(4,1),neighbours(4,2))];
        neighbourvals=neighbourvals+[potmap(neighbours(1,1),neighbours(1,2)),potmap(neighbours(2,1),neighbours(2,2)),potmap(neighbours(3,1),neighbours(3,2)),potmap(neighbours(4,1),neighbours(4,2))];
        mino=101; %this will store the minimum
        indexes=[]; %and it's index
        
        for q=1:4 %this loop finds minimum(s) and index(es)
            if neighbourvals(q)<mino
                mino=neighbourvals(q);
                indexes=q;
            elseif neighbourvals(q)==mino
                indexes=[indexes q];
            end
        end
        if length(indexes)>1 %if more than one, fflip a coin
            indexes=indexes(randi( length(indexes)));
        end
        
        if delayfac>1
            hc_c=hc*2;
        else
            hc_c=hc;
        end
        
        if H(i,j)>=hc_c % MS triggering mechanism
            timesincelast=0;
%            MS_t=[ MS_t iter]; %store the "time"
            MS_t_counter = MS_t_counter+1;
            MS_t(MS_t_counter) = iter;
           
            
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
      
        track_counter = track_counter + 1;
        itrack(track_counter) = i;
        jtrack(track_counter) = j;
        
        
    end
    %MSTIM{cnt}=MS_t;
    MSTIM{cnt}=MS_t(1:MS_t_counter);

    Tracks{cnt}=[itrack;jtrack];
    cnt=cnt+1;
end
toc
d=[];
for i=1:runs
    %figure
    d=[d diff(MSTIM{i})];
    
    % hist(diff(MSTIM{i}),100)
end
delete(h);
d2=d(d<1200);
hist(d2,50)
[N, X] =hist(d2,50);
Nint = interp1(X,N/sum(N),0:1200);
Nint(isnan(Nint))=0;
figure;plot(0:1200,Nint/sum(Nint));
hold on

%mu    98.2240  126.5473  148.8615
%sig    33.0756   34.1897   44.7815
%tau   416.9893  446.9567  390.4344

ex_gauss_params = [96.8885  128.3117  149.9466;
                    31.2952   35.3211   44.1432;
                    414.3318  446.3398  391.5703;]

p = exgausspdf(ex_gauss_params(1,cont_ind), ex_gauss_params(2,cont_ind), ex_gauss_params(3,cont_ind), 0:1200);

plot(p);
hold off
params=egfit(d2)
