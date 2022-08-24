function plotBarsMeanDiffHandAngle_Differential_indiv(a,col)
% a- the data of all participants
% col- colors

set(groot, 'DefaultAxesXColor', [0,0,0], ...
           'DefaultAxesYColor', [0,0,0], ...
           'DefaultAxesZColor', [0,0,0])
       
a_s=a.indiv;

nS=size(a_s,1); % number of participants;

blocx=[1-.2 1+.2 2-.2 2+.2];

lw=2;
tfs=18;

xLim=[.5 2.5];
yLim=[-3.4 5];
yTick=-6:2:6;

figure('position',[50 100 252 400])
hold on
plot(xLim,[0 0],':','color',[1 1 1]*.2,'linewidth',lw)

col_i=nan(nS,3);
for rgb=1:3
    col_i(:,rgb)=linspace(col(1,rgb),col(2,rgb),nS);
end
% individuals' data
for s=1:nS
    plot(blocx, a_s(s,:) ,'-','color',col_i(s,:),'linewidth',1);
end

set(gca,'xtick',blocx,'xticklabel',{'CS+','CS-','CS+','CS-'},'ytick',yTick,'fontsize',tfs,'tickdir','out')
xtickangle(-60)
xlim(xLim);
ylim(yLim);

end

