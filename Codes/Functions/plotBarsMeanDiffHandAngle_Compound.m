function plotBarsMeanDiffHandAngle_Compound(a,tr_c,col)
% a- the data of all participants
% tr_c- trials for each condition
% col- colors

set(groot, 'DefaultAxesXColor', [0,0,0], ...
           'DefaultAxesYColor', [0,0,0], ...
           'DefaultAxesZColor', [0,0,0])
       
co_afHit= col(1,:);
co_afErr= col(2,:);

nC=size(tr_c,2); % number of conditions;

% mean and standard error of each condition
ma_s=a.m;
sea_s=a.se;

vara_s=[ma_s-sea_s;ma_s+sea_s];

bw=.25;
blocx=[1-.2 1+.2 2-.2 2+.2];

lw=2;
afs=24;
tfs=18;

xLim=[.5 2.5];
yLim=[-0.6 0.8];

yTick=-6:.4:6;
yLabel='\DeltaHeading Angle (deg)';

figure('position',[50 100 300 400])
hold on
plot(xLim,[0 0],':','color',[1 1 1]*.2,'linewidth',2)
for c=1:nC
    if (c==1 || c==2)
        co=co_afErr;
    else
        co=co_afHit;
    end
    
    if mod(c,2)
        bar(blocx(c),ma_s(c),'basevalue',0,'barwidth',bw,'linewidth',lw, 'edgecolor',co,'facecolor',co);
    else
        bar(blocx(c),ma_s(c),'basevalue',0,'barwidth',bw,'linewidth',lw, 'edgecolor',co,'facecolor','none');
    end
    
    plot([blocx(c) blocx(c)],vara_s(:,c),'-k','linewidth',lw)
    
end

% add p values
txtfs=16;
xtxtloc=xLim(1)+.5*(xLim(2)-xLim(1));
ydisttxt=.1*(yLim(2)-yLim(1));
ydistline=.05*(yLim(2)-yLim(1));
ytxt_loc1=yLim(2);

if a.pPrev<0.05
    if a.pPrev<0.001
        text(xtxtloc,ytxt_loc1,'p<0.001','fontsize',txtfs,'horizontalalignment','center','fontweight','bold')
    else
        text(xtxtloc,ytxt_loc1,sprintf('p=%.3f',a.pPrev),'fontsize',txtfs,'horizontalalignment','center','fontweight','bold')
    end
else
    text(xtxtloc,ytxt_loc1,sprintf('p=%.3f',a.pPrev),'fontsize',txtfs,'horizontalalignment','center')
end
plot([1 2],ones(1,2)*ytxt_loc1-ydistline,'-k','linewidth',3)
for pl=1:2
    plot([pl-bw pl+bw],ones(1,2)*ytxt_loc1-2*ydistline,'-k','linewidth',1)
    plot([pl pl],[ytxt_loc1-2*ydistline ytxt_loc1-ydistline],':k','linewidth',2)
end

if a.pCurr<0.05
    if a.pCurr<0.001
        text(xtxtloc,ytxt_loc1-1.5*ydisttxt,'p<0.001','fontsize',txtfs,'horizontalalignment','center','fontweight','bold')
    else
        text(xtxtloc,ytxt_loc1-1.5*ydisttxt,sprintf('p=%.3f',a.pCurr),'fontsize',txtfs,'horizontalalignment','center','fontweight','bold')
    end
else
    text(xtxtloc,ytxt_loc1-1.5*ydisttxt,sprintf('p=%.3f',a.pCurr),'fontsize',txtfs,'horizontalalignment','center')
end
plot([(1-bw/2+diff([1-bw/2 2-bw/2])/2) (1+bw/2+diff([1+bw/2 2+bw/2])/2)],ones(1,2)*ytxt_loc1-1.5*ydisttxt-ydistline,'-k','linewidth',3)
plot([1-bw/2 2-bw/2],ones(1,2)*ytxt_loc1-1.5*ydisttxt-2*ydistline,'-k','linewidth',1)
plot(ones(1,2)*1-bw/2+diff([1-bw/2 2-bw/2])/2,[ytxt_loc1-1.5*ydisttxt-2*ydistline ytxt_loc1-1.5*ydisttxt-ydistline],':k','linewidth',2)
plot([1+bw/2 2+bw/2],ones(1,2)*ytxt_loc1-1.5*ydisttxt-3*ydistline,'-k','linewidth',1)
plot(ones(1,2)*1+bw/2+diff([1+bw/2 2+bw/2])/2,[ytxt_loc1-1.5*ydisttxt-3*ydistline ytxt_loc1-1.5*ydisttxt-ydistline],':k','linewidth',2)

set(gca,'xtick',blocx,'xticklabel',{'Comp','Single','Comp','Single'},'ytick',yTick,'fontsize',tfs)
xtickangle(-60)
ylabel(yLabel,'fontsize',afs)
xlim(xLim);
ylim(yLim);

end

