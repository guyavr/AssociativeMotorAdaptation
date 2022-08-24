%% Model comparison (Fig. 4)
clc; clear; close all

addpath('../Functions')

% differential conditioning paradigm
ntrials = 800; % total number of trials
nlearning = 600; % learning trials
nprobe = ntrials-nlearning; % washout

% load data
load data_diffCond;
[nsubs,ntrials] = size(diff_cond.ha);

niter = 200; % fitting iterations
    
newFit_RW=0;
newFit_SS=0;

if newFit_RW
    %%%%%%%%%%%%%%%%%%%%%
    %% RESCORLA-WAGNER %%
    %%%%%%%%%%%%%%%%%%%%%
    
    for s = 1:nsubs
        disp(['now fitting subject ',num2str(s),' rescorla-wagner']);
        ha = diff_cond.ha(s,:);
        us = diff_cond.cs_p(s,:); % cs schedule
        
        for iter = 1:niter
            
            % RW model parameter init
            alpha_plan = rand; % plan salience
            alpha_tone = rand; % tone salience
            alpha_light = rand; % light salience, set to zero because no US.
            beta = rand; % learning rate
            lambda = unifrnd(-30,60); % maximum associative strength over all stimuli
            
            params = [alpha_plan,alpha_tone,alpha_light,beta,lambda];
            options=optimset('display','off');
            LB = [0 0 0 0 -30];
            UB = [1  1 1 1 60];
            [ps, sse] = fmincon(@func_RW,params,[],[],[],[],LB,UB,[],options,ha,us);
            mod.ps(iter,:) = ps;
            mod.sse(iter) = sse;
        end
        [fitRW.sse(s),best] = min(mod.sse);
        fitRW.alpha_plan(s) = mod.ps(best,1);
        fitRW.alpha_tone(s) = mod.ps(best,2);
        fitRW.alpha_light(s) = mod.ps(best,3);
        fitRW.beta(s) = mod.ps(best,4);
        fitRW.lambda(s) = mod.ps(best,5);
        fitRW.aic(s) = 2*length(params) + ntrials*log(fitRW.sse(s)/ntrials);
        
        % test the model
        V = 0;
        V_plan = 0;
        V_tone = 0;
        V_light = 0;

        for t = 1:ntrials
            SIM(t) = V;
            if ~isnan(ha(t))
                % clamp trial?
                if us(t)
                    err = fitRW.lambda(s); % error
                else
                    err = 0;
                end
                
                if us(t)
                    V = V_plan + V_tone;
                    % update tone only
                    delta_tone = fitRW.alpha_tone(s) * fitRW.beta(s) * (err - V);
                    V_tone = V_tone + delta_tone;
                else
                    V = V_plan + V_light;
                    % update light only
                    delta_light = fitRW.alpha_light(s) * fitRW.beta(s) * (err - V);
                    V_light = V_light + delta_light;
                end
                % update plan (every trial)
                delta_plan = fitRW.alpha_plan(s) * fitRW.beta(s) * (err - V);
                V_plan = V_plan + delta_plan;
            end
        end
        
        % plot the best fit
        figure
        hold on
        plot(1:ntrials,ha,'or')
        plot(1:ntrials,SIM,'-k')
    end
    save fitRW fitRW;
    clear mod;
        
    % plot params
    figure;
    subplot(1,5,1);hold on;
    histogram(fitRW.alpha_plan);title('\alpha - plan');
    ylabel('Freq');
    subplot(1,5,2);hold on;
    histogram(fitRW.alpha_tone);title('\alpha - tone');
    subplot(1,5,3);hold on;
    histogram(fitRW.alpha_light);title('\alpha - light');
    subplot(1,5,4);hold on;
    histogram(fitRW.beta);title('\beta');
    subplot(1,5,5);hold on;
    histogram(fitRW.lambda);title('\lambda');
    
else
    load fitRW
end

if newFit_SS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STATE-SPACE (1 state) %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for s = 1:nsubs
        disp(['now fitting subject ',num2str(s),' state-space']);
        ha = diff_cond.ha(s,:); % 1:600 for just learning
        us = diff_cond.cs_p(s,:); % cs schedule
            
        for iter = 1:niter
            
            % RW model parameter init
            A = rand; % retention
            B = rand; % lr
            err = unifrnd(-30,60); % maximum asymtote
        
            params = [A,B,err];
            options=optimset('display','off');
            LB = [0 0 -30];
            UB = [1 1 60];
            [ps, sse] = fmincon(@func_SS,params,[],[],[],[],LB,UB,[],options,ha,us);
            mod2.ps(iter,:) = ps;
            mod2.sse(iter) = sse;
        end
        [fitSS.sse(s),best] = min(mod2.sse);
        fitSS.A(s) = mod2.ps(best,1);
        fitSS.B(s) = mod2.ps(best,2);
        fitSS.err(s) = mod2.ps(best,3);
        fitSS.aic(s) = 2*length(params) + ntrials*log(fitSS.sse(s)/ntrials);
        
        % test the model
        X = zeros(1,ntrials);
        clamp_size = 15;
        for t = 1:(ntrials-1)
            if ~isnan(ha(t))
                if t <= nlearning
                    spe = fitSS.err(s)*us(t) - X(t);
                else
                    spe = 0 - X(t);
                end
                
                X(t+1) = fitSS.A(s)*X(t) + fitSS.B(s)*spe;
            end
        end
        
        % plot the best fit
        figure
        hold on
        plot(1:ntrials,ha,'o')
        plot(1:ntrials,X,'-k')
        
        
    end
    save fitSS fitSS;
    clear mod;
    
    % plot params
    figure;
    subplot(1,2,1);hold on;
    histogram(fitSS.A);title('Retention');
    ylabel('Freq');
    subplot(1,2,2);hold on;
    histogram(fitSS.B);title('Learning Rate');

else
    load fitSS
end

%%%%%%%%%%%%%%%%%
%% Comparisons %%
%%%%%%%%%%%%%%%%%

% Fig. 4
set(groot, 'DefaultAxesXColor', [0,0,0], ...
           'DefaultAxesYColor', [0,0,0], ...
           'DefaultAxesZColor', [0,0,0])
       
col = [255,180,100]/255; % orange
figure('position',[50 100 1000 400]);

afs=24;
tfs=18;

subplot(1,2,1);
[sortedSSE,order] = sort(fitRW.sse-fitSS.sse);
bar(sortedSSE,'facecolor',col,'edgecolor',col,'linewidth',1);
ylabel('\Delta SSR (RW-SS)','fontsize',afs);
xlabel('Participant','fontsize',afs);
box off
set(gca,'xtick',0:5:20,'ytick',-2e5:1e5:1e5,'fontsize',tfs);
xlim([0 17])
ylim([-2.9e5 0.5e5])

subplot(1,2,2);
bar(fitRW.aic(order)-fitSS.aic(order),'facecolor',col,'edgecolor',col,'linewidth',1);
ylabel('\Delta AIC (RW-SS)','fontsize',afs);
xlabel('Participant','fontsize',afs);
box off
set(gca,'xtick',0:5:20,'ytick',-800:400:0,'fontsize',tfs);
xlim([0 17])
ylim([-1100 200])

SSE.mDiff=mean(fitRW.sse-fitSS.sse);
[SSE.h,SSE.p,SSE.ci,SSE.stats] = ttest(fitRW.sse,fitSS.sse);
SSE.cohend = computeCohen_d(fitRW.sse, fitSS.sse, 'paired');

AIC.mDiff=mean(fitRW.aic-fitSS.aic);
[AIC.h,AIC.p,AIC.ci,AIC.stats] = ttest(fitRW.aic,fitSS.aic);
AIC.cohend = computeCohen_d(fitRW.aic, fitSS.aic, 'paired');
