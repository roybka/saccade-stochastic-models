function [ train,pdf,probs ] = generate_train_with_probs( mu,sig,B,B_dur,r,T )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

x=false(1,T);
tau=randn*sig+mu;
lastsactime=0;
r=r/1000;
for t=2:T
    p_new=r;
    
    timesincelast=t-lastsactime;
    if timesincelast<(B_dur)
        p_new=r*B;
    end
    if timesincelast>tau
    x(t)=rand<p_new;
            if x(t)
                lastsactime=t;
                tau=randn*sig+mu;
                %taus=[taus tau];
            end
    else
        x(t)=0;
    end
    if timesincelast<tau
        p_new=0;
    end
    probs(t)=p_new;
end

train=x;
pdf=diff(find(x));

end

