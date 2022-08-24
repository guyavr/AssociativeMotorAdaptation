function plotDiffHandAngle_Differential_sim(a)
% a- the data of all דןצוךשאקג participants

set(groot, 'DefaultAxesXColor', [0,0,0], ...
    'DefaultAxesYColor', [0,0,0], ...
    'DefaultAxesZColor', [0,0,0])

co_afHit= [121,161,197]/255; % blue
co_afErr= 0.6*co_afHit; % darker color

nS=size(a,1); % number of participants;
nC=size(a,2); % number of conditions;

a_s=a;
% mean of each condition (each cell)
ma_s=nanmean(a_s);

bw=.25;
blocx=[1-.2 1+.2 2-.2 2+.2];

lw=2;
afs=24;
tfs=18;

xLim=[.5 2.5];
yLim=[-1.7 2.5];

yTick=-6:1:6;

yLabel='\DeltaHeading Angle (deg)';

yBase=mean(ma_s);

figure('position',[50 100 274 400])
hold on
plot(xLim,yBase*[1 1],':','color',[1 1 1]*.2,'linewidth',2)
for c=1:nC
    if (c==1 || c==2)
        co=co_afErr;
    else
        co=co_afHit;
    end
    
    if mod(c,2)
        bar(blocx(c),ma_s(c),'basevalue',yBase,'barwidth',bw,'linewidth',lw, 'edgecolor',co,'facecolor',co);
    else
        bar(blocx(c),ma_s(c),'basevalue',yBase,'barwidth',bw,'linewidth',lw, 'edgecolor',co,'facecolor','none');
    end
    
end

set(gca,'xtick',blocx,'xticklabel',{'CS+','CS-','CS+','CS-'},'ytick',yTick,'fontsize',tfs)
xtickangle(-60)
ylabel(yLabel,'fontsize',afs)
xlim(xLim);
ylim(yLim);

% subtract according to the error in the previous trial (showing
% the change that is contibuted by the cue)
da_s=nan(nS,2);
da_s(:,1)=a_s(:,1)-a_s(:,2);
da_s(:,2)=a_s(:,3)-a_s(:,4);

end

