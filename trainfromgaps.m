function [ train ] = trainfromgaps( gaps )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
train=false(1,sum(gaps)+1000);
cur=500;
train(cur)=1;
for i=gaps
    cur=cur+i;
    train(cur)=1;
end

end

