function plotDiffHandAngle_CompoundToneLight_sim(a)
% a- the data of all simulated participants

set(groot, 'DefaultAxesXColor', [0,0,0], ...
    'DefaultAxesYColor', [0,0,0], ...
    'DefaultAxesZColor', [0,0,0])

co_afHit= [174,89,116]/255; % moderate pink


a_s=a;

nC=size(a_s,2); % number of conditions;

% mean of each condition (each cell)
ma_s=nanmean(a_s);

bw=.5;
blocx=1:3;

lw=2;
afs=24;
tfs=18;

xLim=[.5 3.5];

yTick=-0.8:.4:0.8;
yLabel='\DeltaHeading Angle (deg)';

yLim=[-0.56 0.84];

figure('position',[50 100 274 400])
hold on
plot(xLim,[0 0],':','color',[1 1 1]*.2,'linewidth',2)

for c=1:nC
    co=co_afHit;
    
    if c==1
        bar(blocx(c),ma_s(c),'basevalue',0,'barwidth',bw,'linewidth',lw, 'edgecolor',co,'facecolor',co);
    else
        bar(blocx(c),ma_s(c),'basevalue',0,'barwidth',bw,'linewidth',lw, 'edgecolor',co,'facecolor','none');
    end
    
end

set(gca,'xtick',blocx,'xticklabel',{'Comp','Tone','Light'},'ytick',yTick,'fontsize',tfs,'tickdir','out')
xtickangle(-60)
ylabel(yLabel,'fontsize',afs)
xlim(xLim);
ylim(yLim);

end

