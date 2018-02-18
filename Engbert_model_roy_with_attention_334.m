function MSTS=Engbert_model_roy_with_attention_334();


%Engbert's model from 2011 paper/ 
visualize=0;


MSTIM=[];cnt=1;
%Parameters::
L=51;
epsilon=.001;
hc=7.9;
lambda=1;

lambda1=.2;
lambda2=.7;
ro1=0.0002;
ro2=0.02;
beta=.3;
kappa=1;
tauA=30;
tauP=150;
iterations=10000;

runs=2000;
% appre=zeros(1,iterations);
% for yy=1:iterations
% appre(yy)=ap(yy);
% end

for iternum=ones(1,runs)*iterations;
%iternum=1000;


%Initialize maps and variables:
H=zeros(L);%100,L,L);
startpos=[randi(L),randi(L)];
i=startpos(1);
j=startpos(2);
i0=round((L-1)/2);
j0=round((L-1)/2);
MS_t=[];
potmaporig=zeros(L);

j_m1=0;i_m1=0;
   itrack=[];
   jtrack=[];
iter=0;timesincelast=0;
for w=1:L
    for ww=1:L
        potmaporig(w,ww)=lambda*L*(((w-i0)/i0)^2+((ww-j0)/j0)^2);
    end
end

while iter<iternum
    iter=iter+1;
    timesincelast=timesincelast+1;
    
%  if mod(iter,5000)==0
%      disp('-')
%  end
 iterap=ap(timesincelast);
potmap=potmaporig*iterap;



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
   thresh=hc*Hc(timesincelast);
   tholds(iter)=thresh;
   if H(i,j)>=thresh % MS triggering mechanism
       MS_t=[ MS_t iter]; %store the "time"
       %disp('MS!!')
       %Nextt lines help calcuate global min as a destination for the MS.
       if iter==7000
           
       timesincelast=0;
       end
       K=H+potmap;
       mino=100;
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
  pause(0.001);
   end
%    i_m1=i;
%    j_m1=j;
itrack=[itrack i];
jtrack=[jtrack j];

end
MSTIM{cnt}=MS_t(MS_t>5000);
Tracks{cnt}=[itrack;jtrack];
cnt=cnt+1;
end


    function out=Cfunc(in)
        out=lambda1*exp(-ro1*in^2);
    end

    function out=ap(t)
        out=1/(1+Cfunc(t-tauP));
    end
    function out=aA(t)
        if t<0
            out=1;
        end
        out=1/(1+Dfunc(t-(tauP+tauA)));
    end
    function out=Dfunc(t)
        out=lambda2*((ro2^kappa)/factorial(kappa+1))*(t^kappa)*exp(-ro2*t);%WTF??? pg. 8036 eq.6
    end

    function out=Hc(t)
        out=1/(1+beta*((1-f(t))+(1-aA(t))));
    end
    function out=f(t)
        out=ap(t-tauP);
    end
MSTS=MSTIM;
save('try1')
end