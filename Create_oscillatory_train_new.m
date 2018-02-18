function [ train ] = Create_oscillatory_train_new(f,l,m,r0 )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
Fs=1000;
t=1:l;

    p=(r0*(m*cos(2*pi*t*(f/Fs))+1))/Fs;

train=false(1,l);

a=rand(1,l);
train( (a)<p)=1;
train=single(train);
end

