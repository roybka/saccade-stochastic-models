subjs=[1,3:4,7:9, 11:13 ,15,16]%15,16];
cnt=0;
for i=subjs
    cnt=cnt+1;
    a=res{cnt};
    for j=1:6
        [val,inds]=sort(a(:,5+j));
        bestparams{cnt,j}=a(inds(1),2:5);
        mus(cnt,j)=a(inds(1),4);
        sigs(cnt,j)=a(inds(1),5);
         hcs(cnt,j)=a(inds(1),3);
    end
end