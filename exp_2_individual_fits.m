%% simulation of exp2 (with burst and gaussian refp)


i=9
params{7}.sigm=25;
params{7}.tau2=140;
params{7}.burst=4;
params{7}.pfac=.8;
params{7}.refpend=300;

params{1}.sigm=75;
params{1}.tau2=250;
params{1}.burst=1;
params{1}.pfac=1.3;
params{1}.refpend=300;

params{9}.sigm=20;
params{9}.tau2=130;
params{9}.burst=5.1;
params{9}.pfac=.9;
params{9}.refpend=240;




burst=params{i}.burst;%bursts(i);
xs={};cnt=0;
sigm=[params{i}.sigm];
tau2=[ params{i}.tau2]; % mean ref. period in ms
cnt=cnt+1;
x=false(1,max(data{i}.c_allsacps(:,7)+1000)); %exp. time
nsacs=length(data{i}.c_allsacps);
rate=nsacs/(length(x)/1000)
taus=[];

lastsactime=-1000;
p=params{i}.pfac*rate/1000; %saccades per second (excluding blinks, etc. )


tau=randn*sigm+tau2;
for t=2:length(x)
    p_new=p;
    
    timesincelast=t-lastsactime;
    if timesincelast<(params{i}.refpend)
        p_new=p*burst;
    end
    if timesincelast>tau
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
% m_minds(i)=modulationIndex;
% m_mfreqs(i)=peakFreq;
% 
% sp_trains{i}=spikeTrain;
a=find(spikeTrain);
b=diff(a);
b=b(b<1500);
pd=smoothy(gaps2flat(b),50);
pd=pd./sum(pd);
figure;
subplot(3,1,1)
plot(pd)

realpd=smoothy(gaps2flat(collapsed_ms_gaps{i}),50);
realpd=realpd./sum(realpd)
hold on
plot(realpd,'r')
title(['PDF  s=' num2str(i) 'm=' num2str(modulationIndex)]) 
legend('model','real')
realac=autocorr_from_pdf_data(realpd);
realcd=cumsum(realpd);
realsu=1-realcd;
realhaz=realpd./realsu;


ac=autocorr_from_pdf_data(pd);
cd=cumsum(pd);
su=1-cd;
haz=pd./su;


subplot(3,1,2);plot(haz(1:1000));
hold on
plot(realhaz(1:1000),'r')
title('hazard')
subplot(3,1,3)

plot(ac(1:1500))
hold on
plot(realac(1:1500),'r')
title('autocorr')
