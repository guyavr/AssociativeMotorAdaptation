function plotBarsMeanDiffHandAngle_CompoundToneLight_indiv(a,col)
% a- the data of all participants
% col- colors

set(groot, 'DefaultAxesXColor', [0,0,0], ...
           'DefaultAxesYColor', [0,0,0], ...
           'DefaultAxesZColor', [0,0,0])
       
a_s=a.indiv;

nS=size(a_s,1); % number of participants;

blocx=1:3;

tfs=18;

xLim=[.5 3.5];
yLim=[-1.5 2.25];

yTick=-6:1:6;

figure('position',[50 100 252 400])
hold on
plot(xLim,[0 0],':','color',[1 1 1]*.2,'linewidth',2)

col_i=nan(nS,3);
for rgb=1:3
    col_i(:,rgb)=linspace(col(1,rgb),col(2,rgb),nS);
end
% individuals' data
for s=1:nS
    plot(blocx, a_s(s,:) ,'-','color',col_i(s,:),'linewidth',1);
end

set(gca,'xtick',blocx,'xticklabel',{'Comp','Tone','Light'},'ytick',yTick,'fontsize',tfs,'tickdir','out')
xtickangle(-60)
xlim(xLim);
ylim(yLim);


end

