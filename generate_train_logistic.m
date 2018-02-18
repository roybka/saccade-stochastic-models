function [ train,pdf,logf ] = generate_train_logistic( mu,slope,B,B_dur,r,T )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

x=false(1,T);
%tau=randn*sig+mu;
lastsactime=0;
r=r/1000;
function_vals=1./(1+exp(-slope*([1:5000]-mu)));
logf=function_vals;
for t=2:T
    p_new=r;
    %probs(t)=p_new;
    timesincelast=t-lastsactime;
    if timesincelast<(B_dur)
        p_new=r*B;
    end
        p_new=p_new*function_vals(timesincelast);
        x(t)=rand<p_new;
            if x(t)
                lastsactime=t;
%                 tau=randn*sig+mu;
                %taus=[taus tau];
            end

end

train=x;
pdf=diff(find(x));

end

