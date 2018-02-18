function autocorr_from_pdf(mu,sig,tau)
%this function generates an cotocorrelogram from exgaussian pdf. 

Acor=zeros(1,10000);
l=5000;
%mu=213;sig=35;tau=230;
ppdf=exgausspdf(mu,sig,tau,1:l);
Acor(1:l)=ppdf;
figure;plot(ppdf)
for i=1:length(ppdf)
    Acor(i:i+l-1)= Acor(i:i+l-1)+Acor(i)*ppdf;
    
end
figure;
plot(Acor(1:1500));
hold on
plot(ppdf,'r')
xlim([0 1500])



end