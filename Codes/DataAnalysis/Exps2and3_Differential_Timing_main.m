%% Associative adaptation - Data analysis - Exps 2 & 3 - CS-US timing in Differential conditioning (web-based)
clc; clear; close all

addpath('../Functions')  

co_afCSm_delay= [255,191,0]/255; % gold
co_afCSm_control= [94,183,123]/255; % green

delayGroup=2; % 1- Delay; 2- Control (matched ITI)
if delayGroup==1
    load('../../Data/AssociativeAdaptation_Exp2_DifferentialTiming_Delay_trials');
else
    load('../../Data/AssociativeAdaptation_Exp3_DifferentialTiming_Control_trials');
end

nS=max(T.SN);
nT=max(T.TN);
BN=T.BN;

clamp_ccw=reshape(T.CCW,nT,[])';
clamp_ccw_s=clamp_ccw(:,1);

block_mat=reshape(BN,nT,[]);
block=block_mat(:,1);

clamp_i=find(block==1,1);
wash_i=find(block==2,1);

tr_clamp=clamp_i:(wash_i-1);
tr_wash=wash_i:nT;

ha_all=T.hand_theta;

% create a matrix of mea (nS X nT)
ha_ccw_cw=reshape(ha_all,nT,[])';

% switch hand angle direction for ccw clamp conditions
ha_cw=ha_ccw_cw;
ha_cw((clamp_ccw_s==1),:)=-ha_cw((clamp_ccw_s==1),:);

% Remove outlier trials
ha = removeOutlierTrials_webBased(ha_cw,wash_i);
ha_timecourse_org = ha; % only for ploting the time course (before removing large RT trials)

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
rt=reshape(T.RT/1000,nT,[])';
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

rt_below_thresh_vals=nan(nS,nT);
rt_below_thresh_vals(rt_above_thresh_mat==0)=rt(rt_above_thresh_mat==0);
mRT_below_thresh=nanmean(rt_below_thresh_vals,'all');
stdRT_below_thresh=nanstd(rt_below_thresh_vals,0,'all');

% Movement Time
mt=reshape(T.MT/1000,nT,[])';
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

nT_outiers_rt_mt=sum(nT_outiers_rt_mt_subj, 1);
per_outiers_rt_mt=100*[nT_outiers_rt_mt(1)/(nS*length(tr_clamp)) nT_outiers_rt_mt(2)/(nS*length(tr_wash))];

% trial by trial *change* in hand angle (the nan column is for consistency
% of the trial number - the change in hand angle according to the previous trial)
dha=[NaN(nS,1) diff(ha,1,2)]; 

% Dissociate trial type. Absolute error trials are according to 'rotation'
% - it is defined also for no feedback trials
rot=reshape(T.ri,nT,[])';
% in this design, the rotation schedule is different across all participants
cs_p=rot;
cs_p(rot~=0)=1;
diff_cond.ha=ha;
diff_cond.dha=dha;
diff_cond.cs_p=cs_p;
diff_cond.rt=rt;

if delayGroup==1
    save('data_diffCond_Timing_VaryingDelay_Comb','diff_cond');
else
    save('data_diffCond_Timing_LongITI_Comb','diff_cond');
end

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
% if delayGroup==1 % for testing stats after removal of extreme participant
%     dha_mC_adapt(64,:)=[];
%     dha_mC_wash(64,:)=[];
% end

% summary analysis
Conditioning.dHA.adapt=summaryAna_differential(dha_mC_adapt);
Conditioning.dHA.wash=summaryAna_differential(dha_mC_wash);

Differential=Conditioning.dHA;
Differential.timeCourse_indiv=ha;
Differential.timeCourse_indiv_org=ha_timecourse_org;

if delayGroup==1
    save('Differential_Timing_Delay','Differential');
else
    save('Differential_Timing_Control','Differential');
end

%% Plots

if delayGroup==1
    co_afCSm= co_afCSm_delay;
else
    co_afCSm= co_afCSm_control;
end
co_afCSp= 0.6*co_afCSm; % darker color
col_timeCourse= 0.8*co_afCSm;
col=[co_afCSm;co_afCSp]; % bar graphs

% Plot mean time courses across participants (Fig. S1A)
plotMeanTimeCourse_Differential(ha_timecourse,t_wash,col_timeCourse)

if delayGroup==1
    plotBarsMeanDiffHandAngle_Differential_Timing_Violin_axisBreak(Conditioning.dHA.adapt,t_n1,col)
else
    plotBarsMeanDiffHandAngle_Differential_Timing_Violin(Conditioning.dHA.adapt,t_n1,col)
end
plotBarsMeanDiffHandAngle_Differential_Timing_Violin(Conditioning.dHA.wash,t_n1,col)

plotBarsMeanDiffHandAngle_Differential_indiv(Conditioning.dHA.adapt,col)
plotBarsMeanDiffHandAngle_Differential_indiv(Conditioning.dHA.wash,col)
