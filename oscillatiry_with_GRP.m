clear

%% model parameters
T=2000;% time in secs. notice this is important according to bargad 2015
m=.4; %the strength of the oscilatory modulation
P=1.8/1000; %saccades per milisecond (excluding blinks, etc. )
f0=4; % frequency of the carrier wave

mu=200; sig=20;
%% run
x=false(1,T*1000); %exp. time
lastsactime=-1000;
cutoff=5000;
time=0.001:0.001:T;
tau=200;
for t=2:length(x)
    timesincelast=t-lastsactime;
    
    p=P*(1+m*cos(2*pi*f0*time(t)));
            if timesincelast>tau
            x(t)=rand<p;
            if x(t)
                lastsactime=t;
                
                 tau=randn*sig+mu;
%                 taus=[taus tau];
            end
            else
                x(t)=0;
            end
end


%% plot
a=find(x);
b=diff(a);
figure;set(gcf, 'Position', get(0,'Screensize'));
subplot(2,2,1)
hist(b(b<1000),100)
% title(['ISIs ' 'm=' num2str(m) '  T=' num2str(T) '  rate=' num2str(P)])
% %title(['sig=' num2str(sigm) 'tau=' num2str(tau2) 'rate=' num2str(p*1000)])
% xlim([0 1000])
% % figure;
% % autocorr(x,700)
subplot(2,2,2)
ar=myautocorr(x,1000,1);
title('autocorr')
%title(['sig=' num2str(sigm) 'tau=' num2str(tau2)])
% a=zeros(length(taus),1000);
% for i=1:length(taus)
%     a(i,1:round(taus(i)))=1;
% end
pd=zeros(1,cutoff);
for i=1:length(b)
    if b(i)<cutoff
        pd(b(i))=pd(b(i))+1;
    end
end
pd=smoothy(pd,50);
pd=pd/sum(pd);
%subplot(2,2,3)
plot(pd(1:1200))
title('pdf')
cd=cumsum(pd);
su=1-cd;
haz=pd./su;
subplot(2,2,4)
plot(haz(1:2000))
title('HAZARD BARGAD')
  