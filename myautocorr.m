function ar=myautocorr(sig,lim,toplot)
fac=10;
a=find(sig);
ar=zeros(1,lim);
for i=11:length(a)
    for k=1:fac
    gap=a(i)-a(i-k);
    if gap<lim
        ar(gap)=ar(gap)+1;
    end
    end
end
ar(1)=0;
if toplot
    plot(smooth(ar/(length(a)-fac),50))
end
end