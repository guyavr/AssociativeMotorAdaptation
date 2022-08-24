%% simulate state-space model
clear;close all;clc;

addpath('../DataAnalysis')
load('data_diffCond')

%% differential conditioning paradigm
nsims = size(diff_cond.ha,1); % number of simulations - for the number of subj
ntrials = size(diff_cond.ha,2); % total number of trials
nlearning = 600; % learning trials
nwash = ntrials-nlearning; % washout

% init state (and set first trial to 0)
x = nan(nsims,ntrials+1,1);  x(:,1) = 0;

% model parameters 
lambda = 15; % error size
A = 0.9;
B= 0.12;

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
        
        x(s,t+1) = A * x(s,t) + B * (err - x(s,t));           
        
    end
    
end

save('SS_sim_ha','x')
