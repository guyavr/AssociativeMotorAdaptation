%% simulate rescorla-wagner model
clear;close all;clc;

addpath('../DataAnalysis')
load('data_diffCond')

%% differential conditioning paradigm
nsims = size(diff_cond.ha,1); % number of simulations - for the number of subj
ntrials = size(diff_cond.ha,2); % total number of trials
nlearning = 600; % learning trials
nwash = ntrials-nlearning; % washout

% init associative values (and set first trial to 0)
V_plan = nan(nsims,ntrials+1,1);  V_plan(:,1) = 0;
V_tone = nan(nsims,ntrials+1,1);  V_tone(:,1) = 0;
V_light = nan(nsims,ntrials+1,1);  V_light(:,1) = 0;
V = nan(nsims,ntrials+1,1);  V(:,1) = 0;

% model parameters - good try for delta HA- 
% Now updated based on the values in sim_SS_RW_diff_comp_current
alpha_plan = 0.99; % plan salience
alpha_tone = .002; % tone salience
alpha_light = .002; % light salience
beta = .12; % learning rate
lambda = 15; % maximum associative strength over all stimuli

% experiment specifications of all subjs
us_m=diff_cond.cs_p;
us_m(:,(nlearning+1):ntrials)=0;

% sim loop
for s = 1:nsims
    
    % trial protocol
    us = us_m(s,:); % clamp trials (after shuffling, with 200 washouts)
    cs_plus = diff_cond.cs_p(s,:); % init tone cs
    
    % trial loop
    for t = 1:ntrials
        % clamp trial?
        if us(t)
            err = lambda;
        else
            err = 0;
        end
        
        % which cs?
        if cs_plus(t)
            V(s,t) = V_plan(s,t) + V_tone(s,t);
            % update tone only
            delta_tone = alpha_tone * beta * (err - V(s,t));
            V_tone(s,t+1) = V_tone(s,t) + delta_tone;
            V_light(s,t+1) = V_light(s,t);
            % update plan
            delta_plan = alpha_plan * beta * (err - V(s,t));
            V_plan(s,t+1) = V_plan(s,t) + delta_plan; 
        else
            V(s,t) = V_plan(s,t) + V_light(s,t);
            % update light only
            delta_light = alpha_light * beta * (err - V(s,t));
            V_light(s,t+1) = V_light(s,t) + delta_light;
            V_tone(s,t+1) = V_tone(s,t);
            % update plan
            delta_plan = alpha_plan * beta * (err - V(s,t));
            V_plan(s,t+1) = V_plan(s,t) + delta_plan;
        end
    end
    
end

save('RW_sim_ha','V')
