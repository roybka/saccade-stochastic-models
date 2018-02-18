i=[1]%:12 ,14,16]
xs={};cnt=0;
 sigm=[50];
    tau2=[ 200]; % mean ref. period in ms
cnt=cnt+1;
x=false(1,length(sp_trains{i})); %exp. time
taus=[];
mvec=.001:.001:(length(x)/1000);
freq=3.5;
m=.4;
mvec=m*cos(2*pi*freq*mvec);
lastsactime=-1000;
p=2/1000; %saccades per second (excluding blinks, etc. )

tau=randn*sigm+tau2;
for t=2:length(x)
    
    timesincelast=t-lastsactime;
    if timesincelast>tau
       p_new=p+p*mvec(t);
       x(t)=rand<p_new;
            if x(t)
                lastsactime=t;
                tau=randn*sigm+tau2;
                taus=[taus tau];
            end
    else
        x(t)=0;
    end
end



spikeTrain=double(x);

totalTime = length(spikeTrain)

a=find(spikeTrain);
b=diff(a);
pd=zeros(1,10000);
for i=b
    pd(i)=pd(i)+1;
end
pd=smoothy(pd,50);
pd=pd./sum(pd)
cd=cumsum(pd);
su=1-cd;
haz=pd./su;

figure;plot(haz)











freq = 4;
fs = 1000;
time=[0:1/fs:totalTime/fs];
%baseRate = 0.03;
%modulation = 0.25;
refPeriod = 0;
winLen = 10000;
% rateFunc = (modulation*baseRate) * (cos(2*pi*freq*time))+baseRate;
% spikeTrain = generatePoissonTrain(totalTime, rateFunc, refPeriod);
firingRate = sum(spikeTrain)/(totalTime/1000);
[spectrum, freqRange,snr, peakPower, peakFreq] = powerSpectrum(spikeTrain, fs, freq-2, freq+2, winLen);
modulationIndex = getModulationIndex(peakPower,firingRate/fs, totalTime , winLen, fs);

% disp('Rate function modulation:');
% disp(modulation);
disp('Extracted modulation index:');
disp(modulationIndex);
m_minds3(i)=modulationIndex;
m_mfreqs3(i)=peakFreq;

