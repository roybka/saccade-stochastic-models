function Acor=autocorr_from_pdf_data(ppdf,varargin)
%this function generates an cotocorrelogram from exgaussian pdf. 

Acor=zeros(1,10000);
l=5000;%length(ppdf);
if ~isempty(varargin)
    l=varargin{1};
end
ppdf=ppdf(1:l);
%ppdf=exgausspdf(mu,sig,tau,1:l);
Acor(1:l)=ppdf(1:l);
%figure;plot(ppdf)
for i=1:length(ppdf)
    Acor(i:i+l-1)= Acor(i:i+l-1)+Acor(i)*ppdf;
    
end
%figure;
%plot(smoothy(Acor(1:1500),50));

end