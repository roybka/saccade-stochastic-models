function []=plotrate(MSTS)
a=zeros(1,10000);
for i=1:length(MSTS)
    for j=1:length(MSTS{i})
        a(MSTS{i}(j))=a(MSTS{i}(j))+1;
    end
end
disp((sum(a)))
a=smoothy(a,100);
plot(a);