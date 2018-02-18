L=51;
% epsilon=.001;
% hc=7.83;
    i0=round((L-1)/2);
    j0=round((L-1)/2);
cnt=0;

for lambda=[.4:.4:4]
  cnt=cnt+1;
potmap=zeros(L);
    for w=1:L
        for ww=1:L
            potmap(w,ww)=lambda*L*(((w-i0)/i0)^2+((ww-j0)/j0)^2);
        end
    end
    
    subplot(5,2,cnt)
    imagesc(potmap);
    a=gca;a.CLim=[0 150]
    title(num2str(lambda))
end


cnt=0;gaps={};mus=[];amps={};
lambdas=repmat([.95:.01:1.15],1,1)
lambdas=[1,1,1,1,1,1];
ms=[ 0.022 0.022 0.022 .024 .024 .024]
parfor k=1:length(lambdas)
  
N=1800*1000
params=[.001,7.9,150,40,lambdas(k)]
[train{k}, isok(k), amps{k}]=Generate_train_Eng_Roy_params_new2(params,N,ms(k));
a=diff(find(train{k}));
gaps{k}=a;
meanamp(k)=mean(amps{k});
    b=a(a<3000);
    egparams=egfit(b);
    mus(k)=egparams(1);
end

[lambdas;mus;meanamp]

r=[];p=[];
for k=1:length(lambdas)
    if isok(k)
    [r(k),p(k)]=corr(amps{k}(2:end)',amps{k}(1:end-1)');
    [r2(k),p2(k)]=corr(amps{k}(1:end-1)',gaps{k}');
    mins(k)=min(amps{k});
     maxz(k)=max(amps{k});
    end
end