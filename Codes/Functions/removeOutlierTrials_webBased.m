function ha = removeOutlierTrials_webBased(ha_org,i_newblock)

spikeThre=20;
nS=size(ha_org,1);
nT=size(ha_org,2);
for s=1:nS
    ha_withOutlierTrials=ha_org(s,:);
    
    trial_out1=find(ha_withOutlierTrials>70 | ha_withOutlierTrials<-50); % with the reasoning that people should adapt ~20deg in the positive direction
    if abs(ha_withOutlierTrials(1))>spikeThre % when the first trial has too large error
        trial_out1=[1 trial_out1];
    end
    
    dha_s=diff(ha_withOutlierTrials); % base on sudden change in hand angle
    diffLarge=find(abs(dha_s)>spikeThre)+1; % find how errors change from trial to trial
    iSpike=find(diff(diffLarge)==1); % a "spike" is when there is an oulier trial in which error was suddenly very different from the trial before and the trial after
    trial_out2=diffLarge(iSpike);
    trial_out2(trial_out2==i_newblock)=[];
    
    trial_out=[trial_out1 trial_out2];
    ha_withOutlierTrials(trial_out)=nan;
    
%     figure
%     plot(1:nT,ha_org(s,:),trial_out,ha_org(s,trial_out),'x')
%     title(['Subject #' num2str(s)],'fontsize',16)
    
    ha(s,:)=ha_withOutlierTrials;
end

end

