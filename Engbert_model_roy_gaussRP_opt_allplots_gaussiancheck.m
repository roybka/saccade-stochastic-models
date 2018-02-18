%clear

%Engbert's model from 2011 paper plus Roy's inhibition between MS addition, works well.
visualize=0;

MSTIM=[];cnt=1;
%Parameters::
L=51;
epsilon=.001;
hc=7.83;
lambda=1;
runs=12;


% mu_sigma(1,:) = [38,115];%[35 107]
% mu_sigma(2,:) = [27,144];%[32,130]
% mu_sigma(3,:) = [34,150];%[42 148]

mu_sigma(1,:) = [38,107];%[35 107]
mu_sigma(2,:) = [46,144];%[32,130]
mu_sigma(3,:) = [52,150];%[42 148]

for cont_ind =3
delayfac=0;
MSTIM=[];cnt=1;
timesincelast=10000;
sigma1=mu_sigma(cont_ind,1);
mu1=mu_sigma(cont_ind,2);

h = waitbar(0,'Running Engbert Amit...');
run_num=-1;
tic
for iternum=ones(1,runs)*180000*12
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
        HS(track_counter)=H(i,j);
        
    end
    %MSTIM{cnt}=MS_t;
    MSTIM{cnt}=MS_t(1:MS_t_counter);

    Tracks{cnt}=[itrack;jtrack];
    cnt=cnt+1;
end

MSTIMCOND{cont_ind}=MSTIM;
delete(h);
end


toc
smooth_trace=[]
for cond=3
d=[];n_sacs=[];
MSTIM=MSTIMCOND{cond};
for i=1:runs
    %figure
    
    d=[d diff(MSTIM{i})];
    n_sacs=n_sacs+length(MSTIM{i});
end
    % hist(diff(MSTIM{i}),100)
    d2=d(d<3000);
trace=zeros(1,3000);
for i=d2
    trace(i)=trace(i)+1;
end
smooth_trace(cond,:)=smoothy(trace,100);

end

ex_gauss_params = [96.8885  128.3117  149.9466;
                    31.2952   35.3211   44.1432;
                    414.3318  446.3398  391.5703;]

figure;
subplot(3,1,1)
plot(1:3000,smooth_trace(1,:)./sum(smooth_trace(1,:)),'b')
p = exgausspdf(ex_gauss_params(1,1), ex_gauss_params(2,1), ex_gauss_params(3,1), 0:3000);
hold on;plot(p,'b--')

subplot(3,1,2)
plot(1:3000,smooth_trace(2,:)./sum(smooth_trace(2,:)),'g')
p = exgausspdf(ex_gauss_params(1,2), ex_gauss_params(2,2), ex_gauss_params(3,2), 0:3000);
hold on;plot(p,'g--')

subplot(3,1,3)
plot(1:3000,smooth_trace(3,:)./sum(smooth_trace(3,:)),'b')
p = exgausspdf(ex_gauss_params(1,3), ex_gauss_params(2,3), ex_gauss_params(3,3), 0:3000);
hold on;plot(p,'r--')
legend('model','theory')

cutoff=3000;
d2=d(d<cutoff);
trace=zeros(1,cutoff);
for i=d2
    trace(i)=trace(i)+1;
end
smooth_trace=smoothy(trace,100);
figure;
plot(1:cutoff,smooth_trace./sum(smooth_trace))
p = exgausspdf(ex_gauss_params(1,cont_ind), ex_gauss_params(2,cont_ind), ex_gauss_params(3,cont_ind), 0:cutoff);
hold on;plot(p,'r')

figure
hist(d2,50)
[N, X] =hist(d2,50);
Nint = interp1(X,N/sum(N),0:cutoff);
Nint(isnan(Nint))=0;
figure;plot(0:cutoff,Nint/sum(Nint));
hold on



%mu    98.2240  126.5473  148.8615
%sig    33.0756   34.1897   44.7815
%tau   416.9893  446.9567  390.4344



p = exgausspdf(ex_gauss_params(1,cont_ind), ex_gauss_params(2,cont_ind), ex_gauss_params(3,cont_ind), 0:1200);

plot(p);
hold off
params=egfit(d2)


ACG=zeros(1,2000);
for k=1:length(MSTIM)
for i=1:length(MSTIM{k})-7
    for l=1:7
        if (MSTIM{k}(i+l)-MSTIM{k}(i))<2000
        ACG((MSTIM{k}(i+l)-MSTIM{k}(i)))=ACG((MSTIM{k}(i+l)-MSTIM{k}(i)))+1;
        end
    end
end
end

figure;plot(smoothy(ACG,100))
xlim([0 1500])

j=1;
for i=1:runs

        b=diff(MSTIM{i});
        ns_b(i,j)=length(b);
        %a=a(a<1500);
        egp{i,j}=egfit(b);
        distfit{i,j}=allfitdist(b,'BIC','PDF');
        mus(i,j)=egp{i,j}(1);
        sigs(i,j)=egp{i,j}(2);
        lambdas(i,j)=egp{i,j}(3);
         logL(i,j)=eglike(egp{i,j}, b);
         EGBIC(i,j)=2*logL(i,j)+3*log(length(b));
         
         
         
         for q=1:15
           if strcmp(distfit{i,j}(q).DistName,'exponential')
               expBIC(i,j)=distfit{i,j}(q).BIC;
           elseif  strcmp(distfit{i,j}(q).DistName,'gamma')
               gammaBIC(i,j)=distfit{i,j}(q).BIC;
           elseif  strcmp(distfit{i,j}(q).DistName,'normal')
               gaussBIC(i,j)=distfit{i,j}(q).BIC;
             elseif  strcmp(distfit{i,j}(q).DistName,  'inverse gaussian')
                  invGBIC(i,j)=distfit{i,j}(q).BIC;
             elseif  strcmp(distfit{i,j}(q).DistName,  'birnbaumsaunders')
                  bimGBIC(i,j)=distfit{i,j}(q).BIC;                  
           end
         end
        
         
end
