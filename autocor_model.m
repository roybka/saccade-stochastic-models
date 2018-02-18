function [trace]=autocor_model (MSTIM)

trace=zeros(length(MSTIM),2000);
for i=1:length(MSTIM)
    for j=1:length(MSTIM{i})
        gap=0;cnt=1;
        while gap < 2000 & ((j+cnt)<numel(MSTIM{i}))
            gap=MSTIM{i}(j+cnt)-MSTIM{i}(j);
            if gap<2000
            cnt=cnt+1;
            trace(i,gap)= trace(i,gap)+1;
            end
        end
    end
end
figure;
if size(trace,1)>1
plot(smoothy(mean(trace),100))
else
plot(smoothy((trace),100))
end    