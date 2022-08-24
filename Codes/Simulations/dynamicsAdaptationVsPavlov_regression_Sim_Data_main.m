%% Dynamics of Pavlovian and adaptation effects (regression analysis)- Simulation and Data- Fig. 5
clear;close all;clc;

addpath('../DataAnalysis')
%% load data
load data_diffCond
load Differential;
diff_cond.ha=Differential.timeCourse_indiv_org;
[nsubs,ntrials] = size(diff_cond.ha);
nlearning = 600; % learning trials
nprobe = ntrials-nlearning; % washout

load('RW_sim_ha') % From the following simulation code: rw_sim_diff_csScheduleOfExp1.m
sim_ha_RW=V;
load('SS_sim_ha') % From the following simulation code: ss_sim_diff_csScheduleOfExp1.m)

sim_ha_SS=x;

%%%%%%%%%%%%%%%%%%%%%
%% RESCORLA-WAGNER %%
%%%%%%%%%%%%%%%%%%%%%

block_l = 50; % trials per phase
Ps = 1:block_l:600+block_l; % lerning blocks
Ps_w = 601:block_l:ntrials+block_l; % testing blocks
Ps_w(end) = 800;

for j = 1:nsubs
    ha = [nan diff(diff_cond.ha(j,1:end))]; % we're just fitting on learning
    
    us = diff_cond.cs_p(j,1:end); % cs schedule
    
    X_learn = [zscore(us(1:600))' zscore(us(2:601))' zscore(us(1:600))'.*zscore(us(2:601))'];
    B_learn(j,:) = regress(ha(2:601)',[X_learn ones(600,1)]);
    
    X_test = [zscore(us(600:799))' zscore(us(601:800))' zscore(us(600:799))'.*zscore(us(601:800))'];
    B_test(j,:) = regress(ha(601:800)',[X_test ones(200,1)]);
        
    %%%%%%%%%%%%%%%%%%%%%%%%
    %% now do it by phase %%
    %%%%%%%%%%%%%%%%%%%%%%%%
    %% DATA %%
    % learning
    for k = 1:nlearning/block_l
        X = [us(Ps(k):Ps(k+1)-1)' us(Ps(k)+1:Ps(k+1))' us(Ps(k):Ps(k+1)-1)'.*us(Ps(k)+1:Ps(k+1))'];
        hatmp = ha(Ps(k)+1:Ps(k+1));
        tmp = regress(hatmp',[X ones(block_l,1)]);
        Bs_prev(j,k) = tmp(1);
        Bs_cur(j,k) = tmp(2);
    end
    
    % testing/washout
    for k = 1:nprobe/block_l
        X = [us(Ps_w(k)+1:Ps_w(k+1))']; % only pavlovian
        hatmp = ha(Ps_w(k)+1:Ps_w(k+1));
        bl = block_l;
        if k == nprobe/block_l
            bl = block_l-1;
        end
        tmp = regress(hatmp',[X ones(bl,1)]);
        Bs_cur_w(j,k) = tmp(1);
    end
    
    %% SIMS %%
    % RW
    ha_sim_RW = [nan diff(sim_ha_RW(j,:))];
    
    % learning
    for k = 1:nlearning/block_l
        X = [us(Ps(k):Ps(k+1)-1)' us(Ps(k)+1:Ps(k+1))' us(Ps(k):Ps(k+1)-1)'.*us(Ps(k)+1:Ps(k+1))'];
        hatmp = ha_sim_RW(Ps(k)+1:Ps(k+1));
        tmp = regress(hatmp',[X ones(block_l,1)]);
        Bs_prev_sim_RW(j,k) = tmp(1);
        Bs_cur_sim_RW(j,k) = tmp(2);
        
    end
    % testing/washout
    for k = 1:nprobe/block_l
        X = [us(Ps_w(k)+1:Ps_w(k+1))'];
        hatmp = ha_sim_RW(Ps_w(k)+1:Ps_w(k+1));
        bl = block_l;
        if k == nprobe/block_l
            bl = block_l-1;
        end
        tmp = regress(hatmp',[X ones(bl,1)]);
        Bs_cur_sim_w_RW(j,k) = tmp(1);
    end
    
    
    % SS
    ha_sim_SS = [nan diff(sim_ha_SS(j,:))];
    
    % learning
    for k = 1:nlearning/block_l
        X = [us(Ps(k):Ps(k+1)-1)' us(Ps(k)+1:Ps(k+1))' us(Ps(k):Ps(k+1)-1)'.*us(Ps(k)+1:Ps(k+1))'];
        hatmp = ha_sim_SS(Ps(k)+1:Ps(k+1));
        tmp = regress(hatmp',[X ones(block_l,1)]);
        Bs_prev_sim_SS(j,k) = tmp(1);
        Bs_cur_sim_SS(j,k) = tmp(2);
        
    end
    % testing/washout
    for k = 1:nprobe/block_l
        X = [us(Ps_w(k):Ps_w(k+1)-1)' us(Ps_w(k)+1:Ps_w(k+1))' us(Ps_w(k):Ps_w(k+1)-1)'.*us(Ps_w(k)+1:Ps_w(k+1))'];
        hatmp = ha_sim_SS(Ps_w(k)+1:Ps_w(k+1));
        bl = block_l;
        if k == nprobe/block_l
            bl = block_l-1;
        end
        tmp = regress(hatmp',[X ones(bl,1)]);
        Bs_prev_sim_w_SS(j,k) = tmp(1);
        Bs_cur_sim_w_SS(j,k) = tmp(2);
    end
    
end

%% Plot Data
afs=24;
tfs=18;
ms=10;

xLim=[0.8 12.2];
yLim=[-0.5 4.5];
co = [[1 1 1]*180;97 129 158]/255;
co_indiv = [[1 1 1]*200;121,161,197]/255;

lw = 1;
figSize=[50 100 450 400];

bins=1:nlearning/block_l;

figure('position',figSize);
hold on
for p=1:2
    if p==1
        Bs=Bs_prev;
    else
        Bs=Bs_cur;
    end
    
    mBs=nanmean(Bs);
    seBs=nanstd(Bs)/sqrt(nsubs);
    
    varBs=[mBs-seBs;mBs+seBs];
    
    fill([bins flip(bins)],[varBs(1,bins) flip(varBs(2,bins))]',co(p,:),'linestyle','none','facealpha',0.3);
    plot(bins,mBs,'o','markerfacecolor',co(p,:),'markeredgecolor',co(p,:),'markersize',ms,'linewidth',lw);
    
end
h=lsline;
h(1).Color=.9*co(2,:);
h(2).Color=.9*co(1,:);
h(1).LineWidth=5;
h(2).LineWidth=5;
xlim(xLim)
ylim(yLim)
xlabel('Bin (50 trials)','fontsize',afs);
ylabel('Regression \beta','fontsize',afs);
set(gca,'xtick',2:2:12,'ytick',0:1:4,'fontsize',tfs);
box off;

% organize data for linear mixed model regression
nbins=nlearning/block_l;
Bs_prev_vec=reshape(Bs_prev,[],1);
Bs_cur_vec=reshape(Bs_cur,[],1);
Betas=[Bs_prev_vec;Bs_cur_vec];
SN=repmat((1:nsubs)',nbins*2,1);
BIN_singleType=reshape(repmat(1:nbins,nsubs,1),[],1);
BIN=repmat(BIN_singleType,2,1);
Type=[zeros(nsubs*nbins,1);ones(nsubs*nbins,1)];

T = table(SN,BIN,Type,Betas);
save('AdaptPavlovEffects_table','T');
writetable(T,'AdaptPavlovEffects_table.csv');

lme = fitlme(T,'Betas ~ BIN * Type + (1 | SN)','FitMethod','REML');
[beta,betanames,stats_Fixed] = fixedEffects(lme,'DFMethod','satterthwaite');
[B,Bnames,stats_Random] = randomEffects(lme,'DFMethod','satterthwaite');
stats = anova(lme,'DFMethod','satterthwaite');
pVal = coefTest(lme);

%% PLOT SIM

% SS (plotted on top ot the RW sim)
bins=1:nlearning/block_l;
xLim=[0.8 12.2];
yLim=[-0.5 4.5];

figure('position',figSize);
hold on
for p=1:2
    if p==1
        Bs_SS=Bs_prev_sim_SS;
        Bs_RW=Bs_prev_sim_RW;
    else
        Bs_SS=Bs_cur_sim_SS;
        Bs_RW=Bs_cur_sim_RW;
    end
    
    mBs_SS=nanmean(Bs_SS);
    seBs_SS=nanstd(Bs_SS)/sqrt(nsubs);
    
    varBs_SS=[mBs_SS-seBs_SS;mBs_SS+seBs_SS];
    
    mBs_RW=nanmean(Bs_RW);
    seBs_RW=nanstd(Bs_RW)/sqrt(nsubs);
    
    varBs_RW=[mBs_RW-seBs_RW;mBs_RW+seBs_RW];
    
    plot(bins,mBs_SS,':','color',.9*co(p,:),'linewidth',3);
    plot(bins,mBs_RW,'-','color',.9*co(p,:),'linewidth',5);

end
xlim(xLim)
ylim(yLim)
xlabel('Bin (50 trials)','fontsize',afs);
ylabel('Regression \beta','fontsize',afs);
set(gca,'xtick',2:2:12,'ytick',0:1:4,'fontsize',tfs);
box off;
