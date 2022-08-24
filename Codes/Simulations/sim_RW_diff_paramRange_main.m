%% Simulate rescorla wagner model- test parameters ranges (Fig. 3)
clear;close all;clc;

%% differential conditioning paradigm
nsims = 100; % number of simulations
ntrials = 800; % total number of trials
nlearning = 600; % learning trials
nwash = ntrials-nlearning; % washout

% init associative values (and set first trial to 0)
V_plan = nan(nsims,ntrials+1,1);  V_plan(:,1) = 0;
V_tone = nan(nsims,ntrials+1,1);  V_tone(:,1) = 0;
V_light = nan(nsims,ntrials+1,1);  V_light(:,1) = 0;
V = nan(nsims,ntrials+1,1);  V(:,1) = 0;

A = 1;
                    
beta_const=0.2;
lambda_const=20;
beta_range=0.15:0.15:0.45;
lambda_range=5:15:35;

nAlphaRange=40;

alpha_plan_range=linspace(0.05,4,nAlphaRange);
alpha_cs_range=linspace(0.0001,0.001,nAlphaRange);

pavEffect_adapt=cell(2,3); % 2- constant beta/lambda, 3- range beta/lambda
pavEffect_wash=cell(2,3); % 2- constant beta/lambda, 3- range beta/lambda

newSet=0;

if newSet
    
    for b_l=1:2 % 1- beta_const; 2- lambda_const
        
        for b_l_range=1:3
            % model parameters
            if b_l==1
                beta=beta_const;
                lambda=lambda_range(b_l_range);
            else
                beta=beta_range(b_l_range);
                lambda=lambda_const;
            end
            
            pavEffect_adapt_cond=nan(nAlphaRange);
            pavEffect_wash_cond=nan(nAlphaRange);
            
            for a_p=1:nAlphaRange
                alpha_plan = alpha_plan_range(a_p); % plan salience
                for a_cs=1:nAlphaRange
                    alpha_tone = alpha_cs_range(a_cs); % tone salience
                    alpha_light = alpha_cs_range(a_cs); % light salience
                    
                    % sim loop
                    for sim = 1:nsims
                        
                        % trial protocol
                        tmp_us = rectpulse([0,1],nlearning/2); % 50 percent clamp schedule
                        us = [tmp_us(randperm(length(tmp_us))) zeros(1,ntrials-length(tmp_us))]; % clamp trials (after shuffling, with 200 washouts)
                        cs_plan = ones(1,ntrials); % plan CS is assumed to happen EVERY TRIAL
                        cs_tone = zeros(1,ntrials); % init tone cs
                        cs_light = zeros(1,ntrials); % init light cs
                        idx_on = find(us(1:length(tmp_us))==1); % paired cs-us
                        idx_off = find(us(1:length(tmp_us))==0); % unpaired cs-us
                        cs_tone(idx_on) = 1; % set tone
                        cs_light(idx_off) = 1; % set light
                        tmp_wash = rectpulse([0,1],nwash/2); % washout, intermixed cs trial
                        wash = tmp_wash(randperm(length(tmp_wash))); % shuffle
                        cs_tone(length(tmp_us)+1:end) = wash; % set tone
                        cs_light(length(tmp_us)+1:end) = 1-wash; % set light
                        
                        % trial loop
                        for t = 1:ntrials
                            % clamp trial
                            if us(t)
                                err = lambda;
                            else
                                err = 0;
                            end
                            
                            % which cs?
                            if cs_tone(t)
                                V(sim,t) = V_plan(sim,t) + V_tone(sim,t);
                                % update tone only
                                delta_tone = alpha_tone * beta * (err - V(sim,t));
                                V_tone(sim,t+1) = V_tone(sim,t) + delta_tone;
                                V_light(sim,t+1) = V_light(sim,t);
                                % update plan
                                delta_plan = alpha_plan * beta * (err - V(sim,t));
                                V_plan(sim,t+1) = A*V_plan(sim,t) + delta_plan; % A is for the retention factor of the state space
                            else
                                V(sim,t) = V_plan(sim,t) + V_light(sim,t);
                                % update light only
                                delta_light = alpha_light * beta * (err - V(sim,t));
                                V_light(sim,t+1) = V_light(sim,t) + delta_light;
                                V_tone(sim,t+1) = V_tone(sim,t);
                                % update plan
                                delta_plan = alpha_plan * beta * (err - V(sim,t));
                                V_plan(sim,t+1) = A*V_plan(sim,t) + delta_plan;
                            end
                        end
                        
                        %% compute changes in hand angles in each of four conditions
                        % learning
                        tpt=[];tpl=[];lpt=[];lpl=[];
                        for n = 2:nlearning
                            dv = V(sim,n) - V(sim,n-1);
                            if cs_tone(n) == 1 && cs_tone(n-1) == 1
                                tpt = [tpt;dv];
                            elseif cs_tone(n) == 1 && cs_tone(n-1) == 0
                                tpl = [tpl;dv];
                            elseif cs_light(n) == 1 && cs_light(n-1) == 1
                                lpl = [lpl;dv];
                            elseif cs_light(n) == 1 && cs_light(n-1) == 0
                                lpt = [lpt;dv];
                            end
                        end
                        tone_post_tone(sim,1) = nanmean(tpt);
                        tone_post_light(sim,1) = nanmean(tpl);
                        light_post_tone(sim,1) = nanmean(lpt);
                        light_post_light(sim,1) = nanmean(lpl);
                        % washout
                        tpt=[];tpl=[];lpt=[];lpl=[];
                        for n = nlearning+1:ntrials
                            dv = V(sim,n) - V(sim,n-1);
                            if cs_tone(n) == 1 && cs_tone(n-1) == 1
                                tpt = [tpt;dv];
                            elseif cs_tone(n) == 1 && cs_tone(n-1) == 0
                                tpl = [tpl;dv];
                            elseif cs_light(n) == 1 && cs_light(n-1) == 1
                                lpl = [lpl;dv];
                            elseif cs_light(n) == 1 && cs_light(n-1) == 0
                                lpt = [lpt;dv];
                            end
                        end
                        tone_post_tone(sim,2) = nanmean(tpt);
                        tone_post_light(sim,2) = nanmean(tpl);
                        light_post_tone(sim,2) = nanmean(lpt);
                        light_post_light(sim,2) = nanmean(lpl);
                    end
                    
                    dha_adapt=[tone_post_tone(:,1),light_post_tone(:,1), tone_post_light(:,1),light_post_light(:,1)]; %tone-err; light-hit
                    dha_wash=[tone_post_tone(:,2),light_post_tone(:,2), tone_post_light(:,2),light_post_light(:,2)]; %tone-err; light-hit
                    
                    pavEffect_adapt_cond_loop=mean([mean(dha_adapt(:,1)-dha_adapt(:,2)), mean(dha_adapt(:,3)-dha_adapt(:,4))]);
                    pavEffect_wash_cond_loop=mean([mean(dha_wash(:,1)-dha_wash(:,2)), mean(dha_wash(:,3)-dha_wash(:,4))]);
                    
                    pavEffect_adapt_cond(a_cs,a_p)=pavEffect_adapt_cond_loop; % alpha_cs - raws; alpha_plan - columns
                    pavEffect_wash_cond(a_cs,a_p)=pavEffect_wash_cond_loop;
                    
                    V=V(:,1:ntrials);
                    V_plan=V_plan(:,1:ntrials);
                    V_tone=V_tone(:,1:ntrials);
                    V_light=V_light(:,1:ntrials);
                    
                    V_all={V;V_plan;V_tone;V_light};
                    
                end
            end
            
            pavEffect_adapt{b_l,b_l_range}=pavEffect_adapt_cond;
            pavEffect_wash{b_l,b_l_range}=pavEffect_wash_cond;
        end
    end
    save('diff_sim_set_paramRange_rw','pavEffect_adapt','pavEffect_wash')
    
else
    load('diff_sim_set_paramRange_rw')
end

min_pavEffect_adapt=min(min(cell2mat(pavEffect_adapt)));
max_pavEffect_adapt=max(max(cell2mat(pavEffect_adapt)));

min_pavEffect_wash=min(min(cell2mat(pavEffect_wash)));
max_pavEffect_wash=max(max(cell2mat(pavEffect_wash)));

min_pavEffect_rw=min(min_pavEffect_adapt,min_pavEffect_wash);
max_pavEffect_rw=max(max_pavEffect_adapt,max_pavEffect_wash);

save('extreme_pavEffect_RW_sim','min_pavEffect_rw','max_pavEffect_rw')
load('extreme_pavEffect_SS_sim') % load the conditioning effect range from the SS simulation
% scale both state space and rescorla-wagner to the same color scale
min_pavEffect=min(min_pavEffect_rw,min_pavEffect_ss);
max_pavEffect=max(max_pavEffect_rw,max_pavEffect_ss);

%%
clc
close all

set(groot, 'DefaultAxesXColor', [0,0,0], ...
           'DefaultAxesYColor', [0,0,0], ...
           'DefaultAxesZColor', [0,0,0])
       
% for the contour plots. load the actual data
load ('../DataAnalysis/Differential.mat')

for block=1:2 % 1-adapt; 2- wash
    if block==1
        pavEffect=pavEffect_adapt;
        pavEffect_rangeData=Differential.adapt.seda_s_rangeCurr;
    else
        pavEffect=pavEffect_wash;
        pavEffect_rangeData=Differential.wash.seda_s_rangeCurr;
    end
    
    for b_l=1:2 % 1- beta_const; 2- lambda_const 
        figure('position',[50,50,1200,320])
        hold on
        for b_l_range=1:3
            subplot(1,3,b_l_range)
            
            h=imagesc([alpha_plan_range(1) alpha_plan_range(end)],[alpha_cs_range(1) alpha_cs_range(end)],(pavEffect{b_l,b_l_range}),'CDataMapping','scaled',[0 max_pavEffect]);
            set(gca,'YDir','normal','xtick',[alpha_plan_range(1) alpha_plan_range(end)],'ytick',[alpha_cs_range(1) alpha_cs_range(end)],'fontsize',14)
            colormap(brewermap([],'*Spectral'))
            
            hold on
            contour(repmat(alpha_plan_range,nAlphaRange,1),repmat(alpha_cs_range',1,nAlphaRange),pavEffect{b_l,b_l_range},pavEffect_rangeData(1)*[1 1],'color',0*[1 1 1],'linewidth',2)
            contour(repmat(alpha_plan_range,nAlphaRange,1),repmat(alpha_cs_range',1,nAlphaRange),pavEffect{b_l,b_l_range},pavEffect_rangeData(2)*[1 1],'k','linewidth',4)
        end
    end
        
end
