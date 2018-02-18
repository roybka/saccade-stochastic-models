clear

%Engbert's model from 2011 paper/ original, works well.
visualize=0;
figure

MSTIM=[];cnt=1;
%Parameters::
L=51;
epsilon=.001;
hc=7.9;
lambda=1;
runs=2;
delayfac=0;
timesincelast=10000;
k=.999;
mu1=250;
for iternum=ones(1,runs)*180000
    %iternum=1000;
    
    
    %Initialize maps and variables:
    H=zeros(L);%100,L,L);
    startpos=[randi(L),randi(L)];
    i=startpos(1);
    j=startpos(2);
    i0=round((L-1)/2);
    j0=round((L-1)/2);
    MS_t=[];
    potmap=zeros(L);
    for w=1:L
        for ww=1:L
            potmap(w,ww)=lambda*L*(((w-i0)/i0)^2+((ww-j0)/j0)^2);
        end
    end
    j_m1=0;i_m1=0;
    itrack=[];
    jtrack=[];
    
    for iter=1:iternum
        
        neighbours(1,1:2)=[i+1,j];
        neighbours(2,1:2)=[i-1,j];
        neighbours(3,1:2)=[i,j+1];
        neighbours(4,1:2)=[i,j-1];
        for poi=1:4
            if neighbours(poi,1)<1
                neighbours(poi,1)=L;
            elseif neighbours(poi,1)>L
                neighbours(poi,1)=1;
            end
            if neighbours(poi,2)<1
                neighbours(poi,2)=L;
            elseif neighbours(poi,2)>L
                neighbours(poi,2)=1;
            end
        end
        
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
       
        if timesincelast>mu1
            activity=H(i,j);
        else
            activity=H(i,j)*k^(mu1+1-timesincelast);
        end
           
        if activity>=hc % MS triggering mechanism
            timesincelast=0;
            MS_t=[ MS_t iter]; %store the "time"
           % delayfac=randn*sigma1+mu1;
            %disp('MS!!')
            %Nextt lines help calcuate global min as a destination for the MS.
            K=H+potmap;
            mino=150;
            indice=[];
            for p=1:L
                for pp=1:L
                    if K(p,pp)<mino
                        mino=K(p,pp);
                        indice=[p,pp];
                    end
                end
            end
            i=indice(1);
            j=indice(2);
            
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
           % delayfac=delayfac-1;
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
        itrack=[itrack i];
        jtrack=[jtrack j];
        
    end
    MSTIM{cnt}=MS_t;
    Tracks{cnt}=[itrack;jtrack];
    cnt=cnt+1;
end

d=[];
for i=1:runs
    %figure
    d=[d diff(MSTIM{i})];
    
    % hist(diff(MSTIM{i}),100)
end
%delete(h);
d2=d(d<1200);
figure
hist(d2,50)
[N, X] =hist(d2,50);
Nint = interp1(X,N/sum(N),0:1200);
Nint(isnan(Nint))=0;
figure;plot(0:1200,Nint/sum(Nint));
hold on
p = exgausspdf(141.9822, 44.3091, 384.5009, 0:1200);
plot(p);
hold off
