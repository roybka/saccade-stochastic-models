
function [resmat]= model_runner_forfit(snum)
st=GetSecs;
subjs=[1,3:4,7:9, 11:13 ,15,16]%15,16];
load('C:\Users\Owner\Google Drive\Roy- experiments\checkerMicSac\model\finals\expdata_forfit\exp1')

epsilons=.001;%+(-.00002*5:.00002:.00002*5);
hcs=7.93+(-.005*5:.005:.005*5)
s_mus=120:5:220;
s_sigmas=20:10:120;


resmat=zeros(length(snum)*length(hcs)*length(s_mus)*length(s_sigmas),12);
    
    cnt=0;

for i=snum
    
    
    for k=1:length(hcs)
        for u=1:length(s_mus)
            for w=1:length(s_sigmas)
                
                train=Generate_train_Eng_Roy_params([epsilons,hcs(k),s_mus(u),s_sigmas(w)],180000);
                gaps=diff(find(train));
                for j=1:6
                    pop=all_gaps{i,j};
                    chis(j)=test_fit2(gaps,pop);
                end
                cnt=cnt+1;
                resmat(cnt,:)=[i,epsilons,hcs(k),s_mus(u),s_sigmas(w),chis,sum(train)];
            end
        end
    end
    

end
resmat(sum(resmat,2)==0,:)=[];
