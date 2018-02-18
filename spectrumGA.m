figure;
cnt=0;

for snum=[1:12]
spikeTrain=double(c_trains{snum});
totalTime = length(spikeTrain)
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
[spectrum{snum}, freqRange{snum},snr, peakPower, peakFreq] = powerSpectrum(spikeTrain, fs, freq-2, freq+2, winLen);
modulationIndex = getModulationIndex(peakPower,firingRate/fs, totalTime , winLen, fs);
cnt=cnt+1;
%subplot(3,1,cnt)
%plot(freqRange(3:end),spectrum(3:end))
%xlim([0 20])
% a=gca;
% a.Box='Off'
end

for snum=1:12
    sp(snum,:)=spectrum{snum};
end