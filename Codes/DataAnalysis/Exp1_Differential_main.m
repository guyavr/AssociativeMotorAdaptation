%% Associative adaptation - Data analysis - Exp 1 - Differential conditioning (lab-based)
clc; clear; close all

addpath('../Functions')

load('../../Data/AssociativeAdaptation_Exp1_Differential_trials');

nS=max(T.SN);
nT=max(T.TN);
nB=max(T.BN);

clamp_ccw=reshape(T.CCW,nT,[])';
clamp_ccw_s=clamp_ccw(:,1);

block_mat=reshape(T.BN,nT,[]);
block=block_mat(:,1);

cs_p_type_mat=reshape(T.cond,nT,[])';
cs_p_type=cs_p_type_mat(:,1);

clamp_i=find(block==1,1);
wash_i=find(block==2,1);

tr_clamp=clamp_i:(wash_i-1);
tr_wash=wash_i:nT;

ha_all=T.hand_theta_maxv;

% create a matrix of mea (nS X nT)
ha_ccw_cw=reshape(ha_all,nT,[])';

% switch hand angle direction for ccw clamp conditions
ha_cw=ha_ccw_cw;
ha_cw((clamp_ccw_s==1),:)=-ha_cw((clamp_ccw_s==1),:);

% Remove outlier trials
ha = removeOutlierTrials_labBased(ha_cw,wash_i);
ha_timecourse_org = ha; % only for ploting the time course (before removing large RT and MT trials)

per_outliers_allTrials=100*sum(sum(isnan(ha_timecourse_org)))/(nT*nS);
per_outliers_perSubj=100*sum(isnan(ha_timecourse_org),2)/nT;

% Remove bias
for s=1:nS
    if isfinite(ha(s,1))
        ha_timecourse(s,:) = ha_timecourse_org(s,:)-ha_timecourse_org(s,1);
    else
        ha_timecourse(s,:) = ha_timecourse_org(s,:);
        ha_timecourse(s,1) = 0;
    end
end

% Reaction Time
rt=reshape(T.RT,nT,[])';
rt_thresh=0.4; % "start faster" message in the experiment

rt_above_thresh_mat=zeros(nS,nT);
rt_above_thresh_mat(rt>rt_thresh)=1;

tr_high_rt=cell(nS,2); % one column for adaptation and a second for washout
nTr_high_rt=nan(nS,2); % one column for adaptation and a second for washout
for s=1:nS
    tr_high_rt{s,1}=find(rt(s,tr_clamp)>rt_thresh); 
    tr_high_rt{s,2}=find(rt(s,tr_wash)>rt_thresh); 
    nTr_high_rt(s,1)=length(tr_high_rt{s,1});
    nTr_high_rt(s,2)=length(tr_high_rt{s,2});
end

% Movement Time
mt=reshape(T.MT,nT,[])';
mt_thresh=0.3; % "move faster" message in the experiment

mt_above_thresh_mat=zeros(nS,nT);
mt_above_thresh_mat(mt>mt_thresh)=1;

tr_high_mt=cell(nS,2); % one column for adaptation and a second for washout
nTr_high_mt=nan(nS,2); % one column for adaptation and a second for washout
for s=1:nS
    tr_high_mt{s,1}=find(mt(s,tr_clamp)>mt_thresh); 
    tr_high_mt{s,2}=find(mt(s,tr_wash)>mt_thresh); 
    nTr_high_mt(s,1)=length(tr_high_mt{s,1});
    nTr_high_mt(s,2)=length(tr_high_mt{s,2});
end

% remove trials with high reaction time
ha(find(rt_above_thresh_mat))=nan;

% remove trials with high movement time
ha(find(mt_above_thresh_mat))=nan;

rt_mt_above_thresh_comb_mat=rt_above_thresh_mat;
rt_mt_above_thresh_comb_mat(find(mt_above_thresh_mat))=1;

nT_outiers_rt_mt_subj=nan(nS,2); % separated for adaptation (first column) and washout (second)
nT_outiers_rt_mt_subj(:,1)=sum(rt_mt_above_thresh_comb_mat(:,tr_clamp),2);
nT_outiers_rt_mt_subj(:,2)=sum(rt_mt_above_thresh_comb_mat(:,tr_wash),2);

nT_outiers_rt_mt=sum(nT_outiers_rt_mt_subj);
per_outiers_rt_mt=100*[nT_outiers_rt_mt(1)/(nS*length(tr_clamp)) nT_outiers_rt_mt(2)/(nS*length(tr_wash))];

% trial by trial *change* in hand angle (the nan column is for consistency
% of the trial number- the change in hand angle according to the previous trial)
dha=[NaN(nS,1) diff(ha,1,2)]; 

rot=reshape(T.ri,nT,[])';
% in this design, the rotation schedule is different across all participants
cs_p=rot;
cs_p(rot~=0)=1;
diff_cond.ha=ha;
diff_cond.dha=dha;
diff_cond.cs_p=cs_p;
diff_cond.rt=rt;

save('data_diffCond','diff_cond');

% Dissociate trial type. Absolute error trials are according to 'rotation'
% - it is defined also for no feedback trials
t_n1=cell(nS,4); % number of conditions
for s=1:nS
    t_n_err=find(rot(s,:))'; % CS of err trials n
    t_n_hit=setdiff(1:nT,t_n_err)'; % CS of hit trials n
    
    dt_n_err=diff(t_n_err);
    dt_n_hit=diff(t_n_hit);
    
    cs_err=zeros(nT,1);
    cs_err(t_n_err)=1;
    cs_hit=zeros(nT,1);
    cs_hit(t_n_hit)=1;
    
    t_n1_errAerr=t_n_err(find(dt_n_err==1)+1);
    t_n1_errAhit=find(diff(cs_err)==1)+1;
    t_n1_hitAhit=t_n_hit(find(dt_n_hit==1)+1);
    t_n1_hitAerr=find(diff(cs_hit)==1)+1;
    
    % combine the above to a single table
    t_n1(s,:)={t_n1_errAerr,t_n1_hitAerr,t_n1_errAhit,t_n1_hitAhit};
end

t_adapt=tr_clamp(1:tr_clamp(end));
t_wash=tr_wash(2:end); % remove the first trial since it might be contaminated by the feedback on the last clamp trial

% mean for each condition in each block
dha_mC_adapt=mean_cond_differential(dha,t_n1,t_adapt);
dha_mC_wash=mean_cond_differential(dha,t_n1,t_wash);

% summary analysis
Conditioning.dHA.adapt=summaryAna_differential(dha_mC_adapt);
Conditioning.dHA.wash=summaryAna_differential(dha_mC_wash);

% analyze movement time
mt_afExclude=mt;
mt_afExclude(find(mt_above_thresh_mat))=nan;
mMT_adapt=nanmean(mt_afExclude(:,t_adapt),2);
mmMT_adapt=nanmean(mMT_adapt);
semMT_adapt=nanstd(mMT_adapt)/sqrt(nS);

% analyze reaction time
rt_afExclude=rt;
rt_afExclude(find(rt_above_thresh_mat))=nan;
mRT=nanmean(rt_afExclude,'all');
stdRT=nanstd(reshape(rt_afExclude,1,[]));

Differential=Conditioning.dHA;
Differential.timeCourse_indiv=ha;
Differential.timeCourse_indiv_org=ha_timecourse_org;
save('Differential','Differential');

% awaerness: check if participants that identified the pattern showed
% higher pavlovian effect
s_aware=[2,5,6,8,12,13,14];
s_unaware=setdiff(1:nS,s_aware);
adapt_aware_PavEffect=Differential.adapt.da_sCurr(s_aware);
adapt_unaware_PavEffect=Differential.adapt.da_sCurr(s_unaware);
wash_aware_PavEffect=Differential.wash.da_sCurr(s_aware);
wash_unaware_PavEffect=Differential.wash.da_sCurr(s_unaware);

% BayesFactor: A Matlab package for Bayes Factor statistical analysis- https://zenodo.org/record/7006300
% Cohen's D: Ruggero G. Bettinardi (2022). computeCohen_d(x1, x2, varargin) (https://www.mathworks.com/matlabcentral/fileexchange/62957-computecohen_d-x1-x2-varargin), MATLAB Central File Exchange. Retrieved August 19, 2022.
PavEff_awareness_adapt_diffMeans=mean(adapt_unaware_PavEffect)-mean(adapt_aware_PavEffect);
[PavEff_awareness_adapt_h,PavEff_awareness_adapt_p,PavEff_awareness_adapt_ci,PavEff_awareness_adapt_stats] = ttest2(adapt_unaware_PavEffect,adapt_aware_PavEffect);
[PavEff_awareness_adapt_bf10,~] = bf.ttest2(adapt_unaware_PavEffect,adapt_aware_PavEffect);
PavEff_awareness_adapt_cohend = computeCohen_d(adapt_unaware_PavEffect, adapt_aware_PavEffect, 'independent');
PavEff_awareness_wash_diffMeans=mean(wash_unaware_PavEffect)-mean(wash_aware_PavEffect);
[PavEff_awareness_wash_h,PavEff_awareness_wash_p,PavEff_awareness_wash_ci,PavEff_awareness_wash_stats] = ttest2(wash_unaware_PavEffect,wash_aware_PavEffect);
[PavEff_awareness_wash_bf10,~] = bf.ttest2(wash_unaware_PavEffect,wash_aware_PavEffect);
PavEff_awareness_wash_cohend = computeCohen_d(wash_unaware_PavEffect, wash_aware_PavEffect, 'independent');

% Power analysis for Exp 4
dha_CSp_adapt=[Differential.adapt.indiv(:,1); Differential.adapt.indiv(:,3)];
dha_CSm_adapt=[Differential.adapt.indiv(:,2); Differential.adapt.indiv(:,4)];
PavEff_adapt_cohend = computeCohen_d(dha_CSp_adapt, dha_CSm_adapt, 'paired'); 

dha_CSp_wash=[Differential.wash.indiv(:,1); Differential.wash.indiv(:,3)];
dha_CSm_wash=[Differential.wash.indiv(:,2); Differential.wash.indiv(:,4)];
PavEff_wash_cohend = computeCohen_d(dha_CSp_wash, dha_CSm_wash, 'paired'); % for power analysis for main exp

%% Plots

% colors
co_afCSm= [121,161,197]/255; % blue
co_afCSp= 0.6*co_afCSm; % darker color
col_timeCourse= 0.8*co_afCSm;
col=[co_afCSm;co_afCSp];

% Plot mean time courses across participants (Fig. 2B)
plotMeanTimeCourse_Differential(ha_timecourse,t_wash,col_timeCourse)

% Bar plots- summary analysis - Change in Hand Angle (Figs. 2C, 2D)
plotBarsMeanDiffHandAngle_Differential_Violin(Conditioning.dHA.adapt,t_n1,col)
plotBarsMeanDiffHandAngle_Differential_Violin(Conditioning.dHA.wash,t_n1,col)

plotBarsMeanDiffHandAngle_Differential_indiv(Conditioning.dHA.adapt,col)
plotBarsMeanDiffHandAngle_Differential_indiv(Conditioning.dHA.wash,col)
