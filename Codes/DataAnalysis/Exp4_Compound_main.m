%% Associative adaptation - Data analysis - Exp 4 - Compound conditioning (lab-based)
clc; clear; close all

addpath('../Functions')  

load('../../Data/AssociativeAdaptation_Exp4_Compound_trials');

nS=max(T.SN);
nT=max(T.TN);
nB=max(T.BN);

clamp_ccw=reshape(T.CCW,nT,[])';
clamp_ccw_s=clamp_ccw(:,1);

block_mat=reshape(T.BN,nT,[]);
block=block_mat(:,1);

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

tone=reshape(T.cs_tone,nT,[])';
frame=reshape(T.cs_frame,nT,[])';

comp.ha=ha;
comp.dha=dha;
comp.cs_tone=tone;
comp.cs_light=light;
comp.rt=rt;

save('data_compound','comp');

t_n1_all=cell(nS,3); % number of conditions
t_n1_SepCurrCS=cell(nS,3); % number of conditions
t_n1=cell(nS,4); % number of conditions

for s=1:nS
    t_n_comp_all=find(tone(s,:)==frame(s,:))';
    t_n_comp=setdiff(t_n_comp_all(wash_i:end),wash_i); % Removal of trial 601 from all conditions
    t_n_tone=setdiff(find(tone(s,:)), [find(frame(s,:)) wash_i] )';
    t_n_frame=setdiff(find(frame(s,:)), [find(tone(s,:)) wash_i] )';
    
    % combine tone and frame trials to a single CS condition;
    t_n_single=union(t_n_tone,t_n_frame);
    
    dt_n_comp=diff(t_n_comp);
    dt_n_single=diff(t_n_single);
    
    cs_comp=zeros(nT,1);
    cs_comp(t_n_comp_all)=1;
    cs_single=zeros(nT,1);
    cs_single(t_n_single)=1;
    
    cs_comp(tr_clamp)=nan;
    cs_single(tr_clamp)=nan;
    
    t_n1_cAc=t_n_comp(find(dt_n_comp==1)+1);
    t_n1_cAs=find(diff(cs_comp)==1)+1;
    t_n1_sAs=t_n_single(find(dt_n_single==1)+1);
    t_n1_sAc=find(diff(cs_single)==1)+1;
    
    % compound after single (t_n1_cAs), tone after comp, light after comp
    t_n1_toneAc=intersect(t_n_tone,t_n1_sAc);
    t_n1_frameAc=intersect(t_n_frame,t_n1_sAc);
    
    % compound after single (t_n1_cAs), tone after comp/frame, frame after
    % comp/tone
    t_minus_f=tone(s,:)-frame(s,:);
    d_t_minus_f=[nan diff(t_minus_f)];
    t_n1_toneAframe=find(d_t_minus_f==2);
    t_n1_frameAtone=find(d_t_minus_f==-2);
    
    t_n1_toneAother=union(t_n1_toneAc,t_n1_toneAframe);
    t_n1_frameAother=union(t_n1_frameAc,t_n1_frameAtone);
    
    % combine the above to single tables
    t_n1_all(s,:)={t_n_comp_all,t_n_tone,t_n_frame};
    t_n1_SepCurrCS(s,:)={t_n_comp,t_n_tone,t_n_frame};
    t_n1(s,:)={t_n1_cAc,t_n1_sAc,t_n1_cAs,t_n1_sAs};   
    
end

% mean for each condition in each block
dha_mC_SepCurrCS=mean_cond_compound(dha,t_n1_SepCurrCS);
dha_mC=mean_cond_compound(dha,t_n1);

% summary analysis
Conditioning.dHA=summaryAna_compound(dha_mC);
Conditioning.dHA_SepCS=summaryAna_Compound_SepCS(dha_mC_SepCurrCS);
close all

Compound=Conditioning.dHA;
Compound.timeCourse_indiv=ha;
Compound.timeCourse_indiv_org=ha_timecourse_org;
Compound.CompToneLight=Conditioning.dHA_SepCS;
save('Compound','Compound');

%% Plots

co_afHit= [203,115,141]/255; % light pink
co_afErr= [144,62,90]/255; % pink
col=[co_afHit;co_afErr];
col_timeCourse= mean(col);

% Plot mean time courses across participants (Fig. 8B)
plotMeanTimeCourse_Compound(ha_timecourse,tr_wash,col_timeCourse)

% Bar plots- summary analysis - Change in Hand Angle (Figs. 8C, 8E)
plotBarsMeanDiffHandAngle_CompoundToneLight(Conditioning.dHA_SepCS,t_n1_SepCurrCS)
plotBarsMeanDiffHandAngle_CompoundToneLight_indiv(Conditioning.dHA_SepCS,col)

% (Figs. 8F)
plotBarsMeanDiffHandAngle_Compound(Conditioning.dHA,t_n1,col)
plotBarsMeanDiffHandAngle_Compound_indiv(Conditioning.dHA,col)


