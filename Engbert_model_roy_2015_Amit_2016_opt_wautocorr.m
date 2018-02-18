%clear all

%Engbert's model from 2015 paper/ original, works well.
visualize=0;

MSTIM=[];cnt=1;
neighbours = zeros(4,2);

%Parameters:
L=50;
epsilon=.00146;
hc=7.55;
lambda=1;

lambda1=1.48*2;%1.48; %after micro-saccade
lambda2=9.07*0.25; %9.07%after micro-saccade
rho1=0.0052*0.5*0.5;%0.0052; %after micro-saccade
rho2=0.00308*0.5*0.5*0.5; %after micro-saccade


runs=4;

h = waitbar(0,'Running Engbert 2015...');
run_num=-1;

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
    MS_t=[];
    potmap=zeros(L);
    Mmap = zeros(L);
    for w=1:L
        for ww=1:L
            potmap(w,ww)=lambda*L*(((w-i0)/i0)^2+((ww-j0)/j0)^2);
        end
    end
    j_m1=0;i_m1=0;
    itrack=[];
    jtrack=[];
    itrack1=zeros(1,180000);
    jtrack1=zeros(1,180000);
    track_counter=0;
    Ms_t1= zeros(1,180000);
    Ms_t_counter=0;
    
    t_from_last_ms = 0;

    for iter=1:iternum
        
        if(mod(iter,1000)==0)
            waitbar((run_num*iternum + iter)/(iternum*runs),h,['Running Engbert 2015 %' num2str(round((run_num*iternum + iter)/(iternum*runs)*100))])
        end
        
        % update the potmap according to post sccadic behavior function engbert 2015 (8)  
        ap_t = ((1+lambda1*exp(-rho1*t_from_last_ms^2))^-1);
        potmap_t = potmap*ap_t;
        neighbours(1,1:2)=[i+1,j];
        neighbours(2,1:2)=[i-1,j];
        neighbours(3,1:2)=[i,j+1];
        neighbours(4,1:2)=[i,j-1];
        neighbours(neighbours<1)=L;
        neighbours(neighbours>L)=1;
        
        neighbourvals=[H(neighbours(1,1),neighbours(1,2)),H(neighbours(2,1),neighbours(2,2)),H(neighbours(3,1),neighbours(3,2)),H(neighbours(4,1),neighbours(4,2))];
        neighbourvals=neighbourvals+[potmap_t(neighbours(1,1),neighbours(1,2)),potmap_t(neighbours(2,1),neighbours(2,2)),potmap_t(neighbours(3,1),neighbours(3,2)),potmap_t(neighbours(4,1),neighbours(4,2))];
        
        vmin = min(neighbourvals);
        indexes = find(neighbourvals == vmin);
        
        if length(indexes)>1 %if more than one, fflip a coin
            indexes=indexes(randi( length(indexes)));
        end
        
        %b_crit from the engbert 2015 (9)
        hc_t = hc*(1 + (1/(1 + lambda2*exp(rho2*t_from_last_ms^2))));
        if (H(i,j)*ap_t + potmap_t(i,j))>=hc_t % MS triggering mechanism
            
            
            Ms_t1(Ms_t_counter+1) = iter;
            Ms_t_counter = Ms_t_counter+1;
                
            
            %disp('MS!!')
            %Nextt lines help calcuate global min as a destination for the MS.
            
            K=H*ap_t+potmap_t+Mmap*ap_t;
            
            [val, i_ind_vec] = min(K);
            [val, j_ind] = min(val);
            i_ind = i_ind_vec(j_ind);
            
            i = i_ind;
            j = j_ind;
                       
            t_from_last_ms = 0; %reset time
        else %else (no MS ) go to the local minimum we calculated.        
            i=neighbours(indexes,1);
            j=neighbours(indexes,2);
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
        
        track_counter = track_counter + 1;
        itrack1(track_counter) = i;
        jtrack1(track_counter) = j;
        
        t_from_last_ms = t_from_last_ms + 1;    
    end
    
    MSTIM{cnt}=Ms_t1(1:Ms_t_counter);
    Tracks{cnt}=[itrack1;jtrack1];
    
    cnt=cnt+1;    
end

d=[];
for i=1:runs
    %figure
    d=[d diff(MSTIM{i})];
    
    % hist(diff(MSTIM{i}),100)
end
delete(h);
d2=d(d<3000);
figure
hist(d2,50)
[N, X] =hist(d2,50);
Nint = interp1(X,N/sum(N),0:1200);
Nint(isnan(Nint))=0;
figure;%plot(0:1200,Nint/sum(Nint));
hold on
p = exgausspdf(141.9822, 44.3091, 384.5009, 0:1200);
plot(p);
hold off


trace=zeros(1,3000);
for i=d2
    trace(i)=trace(i)+1;
end
smooth_trace=smoothy(trace,100);
hold on
plot(smooth_trace./sum(smooth_trace),'k')
%% autocorrelation
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