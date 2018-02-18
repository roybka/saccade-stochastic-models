function [ train,pdf ] = generate_train_RP_comp( mu,sig,B,B_dur,r,T )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
N=r*T/1000;
RP_time=N*mu;
clean_time=T-RP_time;
new_r=(N/(clean_time/1000));
r=new_r;


x=false(1,T);
tau=randn*sig+mu;
lastsactime=0;
r=r/1000;
for t=2:T
    p_new=r;
    %probs(t)=p_new;
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
end

train=x;
pdf=diff(find(x));

end

