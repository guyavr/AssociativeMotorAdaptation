function plotBarsMeanDiffHandAngle_CompoundToneLight(a,tr_c)
% a- the data of all participants
% tr_c- trials for each condition

set(groot, 'DefaultAxesXColor', [0,0,0], ...
           'DefaultAxesYColor', [0,0,0], ...
           'DefaultAxesZColor', [0,0,0])
       
co= [174,89,116]/255; % moderate pink

a_s=a.indiv;

nC=size(tr_c,2); % number of conditions;

% mean and standard error of each condition
ma_s=a.m;
sea_s=a.se;

vara_s=[ma_s-sea_s;ma_s+sea_s];

bw=.5;
blocx=1:3;

lw=2;
afs=24;
tfs=18;

xLim=[.5 3.5];
yLim=[-0.56 0.84];

yTick=-.8:.4:.8;
yLabel='\DeltaHeading Angle (deg)';

figure('position',[50 100 300 400])
hold on
plot(xLim,[0 0],':','color',[1 1 1]*.2,'linewidth',2)
for c=1:nC
    
    if c==1
        bar(blocx(c),ma_s(c),'basevalue',0,'barwidth',bw,'linewidth',lw, 'edgecolor',co,'facecolor',co);
    else
        bar(blocx(c),ma_s(c),'basevalue',0,'barwidth',bw,'linewidth',lw, 'edgecolor',co,'facecolor','none');
    end
    
    plot([blocx(c) blocx(c)],vara_s(:,c),'-k','linewidth',2)
    
end

% add p values
txtfs=16;
xtxtloc=xLim(1)+.5*(xLim(2)-xLim(1));
ydisttxt=.1*(yLim(2)-yLim(1));
ydistline=.06*(yLim(2)-yLim(1));
ytxt_loc1=yLim(2);

% For Friedman
pVal_CompFrame=a.mltcmp(1,6);
pVal_CompTone=a.mltcmp(2,6);

if pVal_CompFrame<0.05 %compound vs frame
    if pVal_CompFrame<0.001
        text(xtxtloc,ytxt_loc1,'p<0.001','fontsize',txtfs,'horizontalalignment','center','fontweight','bold')
    else
        text(xtxtloc,ytxt_loc1,sprintf('p=%.3f',pVal_CompFrame),'fontsize',txtfs,'horizontalalignment','center','fontweight','bold')
    end
else
    text(xtxtloc,ytxt_loc1,sprintf('p=%.3f',pVal_CompFrame),'fontsize',txtfs,'horizontalalignment','center')
end
plot([1 3],ones(1,2)*ytxt_loc1-ydistline,'-k','linewidth',3)

if pVal_CompTone<0.05 %compound vs frame
    if pVal_CompTone<0.001
        text(xtxtloc,ytxt_loc1,'p<0.001','fontsize',txtfs,'horizontalalignment','center','fontweight','bold')
    else
        text(1.5,ytxt_loc1-1.3*ydisttxt,sprintf('p=%.3f',pVal_CompTone),'fontsize',txtfs,'horizontalalignment','center','fontweight','bold')
    end
else
    text(xtxtloc,ytxt_loc1-1.3*ydisttxt,sprintf('p=%.3f',pVal_CompTone),'fontsize',txtfs,'horizontalalignment','center')
end
plot([1 2],ones(1,2)*ytxt_loc1-1.3*ydisttxt-ydistline,'-k','linewidth',3)

set(gca,'xtick',blocx,'xticklabel',{'Comp','Tone','Light'},'ytick',yTick,'fontsize',tfs,'tickdir','out')
xtickangle(-60)
ylabel(yLabel,'fontsize',afs)
xlim(xLim);
ylim(yLim);

% Plot correlation
a_tl=a_s(:,2:3);
xLim_corr=[-1.5 1];
yLim_corr=xLim_corr;
figure('position',[50 100 400 400])
hold on
plot([-1.5 1],[1 -1.5],':','color',[1 1 1]*.4,'linewidth',3)
pc=scatter(a_tl(:,1),a_tl(:,2),'o','markeredgecolor','none','markerfacecolor',co,'sizedata',200);
alpha 0.7
uistack(pc,'top')
xlabel({'\DeltaHeading Angle (deg)','Tone Alone'},'fontsize',20)
ylabel({'\DeltaHeading Angle (deg)','Light Alone'},'fontsize',20)
set(gca,'xlim',xLim_corr,'ylim',yLim_corr,'xtick',-1:1:1,'ytick',-1:1:1,'fontsize',18,'tickdir','out')

end

