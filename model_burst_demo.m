restab=zeros(2500,7);
%restab=single(restab);
jumps=18;
rates=[1.2 ,1.8 ,2.4];
s_mus=170
s_sigmas=20;
s_Bs=[1:.2:5];
s_Blens=[230:30:330]
n_iters=length(s_mus)*length(s_sigmas)*length(s_Bs)*length(s_Blens)*length(rates)*12;
runtime=(n_iters/1000)*14/60/60;
disp(['estimated runtime is  ' num2str(runtime) '   Hours'])
cnt3=0;

for k=1:length(rates)
    rate=rates(k);
    for mus=s_mus
        for sigmas=s_sigmas
            for Bs=s_Bs
                for Blens=mus:30:330
                    [ train,gaps ] = generate_train( mus,sigmas,Bs,Blens,rate,(2000000) );
                    Nsacs=sum(train);
                    cnt3=cnt3+1;
                    
                    
                    [n,c]=hist([0,gaps(gaps<1000),1000],jumps);
                     n(1)=n(1)-1;
                     n(end)=n(end)-1;
                     n=n./(sum(n));
%                     xi=((n- norig{cnt}).^2)./norig{cnt};
                     restab(cnt3,:)=[mus,sigmas,Bs,Blens,rate,0,Nsacs];
                        
                    
                    
                    if mod(cnt3,1000)==0
                        
                       
                 
                    end
                end
            end
        end
    end
end


inds1=find(((restab(:,5)==1.8)&(restab(:,4)==230)))
partial_tab=restab(inds1,:);
for i=1:length(partial_tab)
     [ train,gaps ] = generate_train( partial_tab(i,1),partial_tab(i,2),partial_tab(i,3),partial_tab(i,4),partial_tab(i,5),(2000000) );
     
     spikeTrain=double(train);

totalTime = length(spikeTrain);

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
ms(1,i)=modulationIndex;
end

figure;plot(partial_tab(:,3),ms(1,:))

inds2=find(((restab(:,5)==1.2)&(restab(:,4)==230)))
partial_tab=restab(inds2,:);
for i=1:length(partial_tab)
     [ train,gaps ] = generate_train( partial_tab(i,1),partial_tab(i,2),partial_tab(i,3),partial_tab(i,4),partial_tab(i,5),(2000000) );
     
     spikeTrain=double(train);

totalTime = length(spikeTrain);

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
ms(2,i)=modulationIndex;
end
hold on
plot(partial_tab(:,3),ms(2,:))


inds3=find(((restab(:,5)==2.4)&(restab(:,4)==230)))
partial_tab=restab(inds3,:);
for i=1:length(partial_tab)
     [ train,gaps ] = generate_train( partial_tab(i,1),partial_tab(i,2),partial_tab(i,3),partial_tab(i,4),partial_tab(i,5),(2000000) );
     
     spikeTrain=double(train);

totalTime = length(spikeTrain);

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
ms(3,i)=modulationIndex;
end
hold on
plot(partial_tab(:,3),ms(3,:))