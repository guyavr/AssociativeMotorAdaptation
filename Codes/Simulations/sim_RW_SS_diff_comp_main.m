%% simulations for state space and Rescorla-Wagner models for differential
%% and compound conditioning experiments. (Figs. 2E, 2F, and 8D)
clc; close all; clear

addpath('../Functions')  

nSims = 20000; % number of simulations SHOULD BE 20000
nAcquis = 600; % acquisition/learning trials

% models parameters
% state space
A = 0.9; % retention
B = 0.15; % learning rate
err = 15; % sensory prediction error (it's a clamp, but we treat it like rotation, i.e., an internal SPE that decreases as adaptation proceeds)

newSet_exp1 = 0;
newSet_exp4 = 0;

nT_Ac=nan(nSims,2); % after compound: column 1 for tone, column 2 for light

for exp= [1,4] % 1- differential; 4- compound
    
    if exp == 1
        
        if ~newSet_exp1
            continue
        end
        
        paramVal_RW=[0.99 0.002 0.002 0.12 15];

        nT = 800; % number of trials  
        nProbe = nT - nAcquis;  
    
        dha_acq_SS = nan(nSims,4);
        dha_probe_SS = nan(nSims,4);
        dha_acq_RW = nan(nSims,4);
        dha_probe_RW = nan(nSims,4);
        
    else
        
        if ~newSet_exp4
            continue
        end
        
        paramVal_RW=[0.99 0.1 0.1 0.02 15];

        nT = 900; % number of trials
        nProbe = nT - nAcquis;  
    
        dha_SS = nan(nSims,3);
        dha_RW = nan(nSims,3);
        
    end
        
    alpha_goal = paramVal_RW(1); % goal salience
    alpha_tone = unifrnd(paramVal_RW(2),paramVal_RW(3),nSims,1); % tone salience
    alpha_light = alpha_tone; % for varying CS across participants but having the same learning rate fot tone and light
    beta = paramVal_RW(4); % learning rate of US
    lambda = paramVal_RW(5); % maximum associative strength over all stimuli (acquivalent to initial SPE)
    
    for sim = 1:nSims
        
        % init state for state space model(and set first trial to 0)
        x_ss = nan(1,nT);  x_ss(1) = 0;
        
        % init associative values (and set first trial to 0)
        V_goal = nan(1,nT,1);  V_goal(1) = 0;
        V_tone = nan(1,nT,1);  V_tone(1) = 0;
        V_light = nan(1,nT,1);  V_light(1) = 0;
        V = nan(1,nT,1); % initializes below as the sum of all presented CSs
        
        % init states for multi-context state space model(and set first trial to 0)
        x_goal = nan(1,nT,1);  x_goal(1) = 0;
        x_tone = nan(1,nT,1);  x_tone(1) = 0;
        x_light = nan(1,nT,1);  x_light(1) = 0;
        x_comp = nan(1,nT,1);  x_comp(1) = 0;

        % trial protocol
        if exp == 1
            tmp_us = rectpulse([0,1],nAcquis/2);
            us = [tmp_us(randperm(length(tmp_us))) zeros(1,nProbe)];
            cs_tone = zeros(1,nT); % init tone cs
            cs_light = zeros(1,nT); % init light cs
            idx_on = find(us(1:length(tmp_us))==1); % cs+
            idx_off = find(us(1:length(tmp_us))==0); % cs-
            cs_tone(idx_on) = 1; % set tone (cs+)
            cs_light(idx_off) = 1; % set light (cs-)
            tmp_probe = rectpulse([0,1],nProbe/2); % probe, intermixed cs trial
            probe_cs = tmp_probe(randperm(length(tmp_probe))); % shuffle
            cs_tone(nAcquis+1:end) = probe_cs; % set tone
            cs_light(nAcquis+1:end) = 1-probe_cs; % set light

        else
            us = [ones(1,nAcquis) zeros(1,nProbe)];
            cs_tone = zeros(1,nT); % init tone cs
            cs_light = zeros(1,nT); % init light cs
            cs_tone(1:nAcquis) = 1; % set tone
            cs_light(1:nAcquis) = 1; % set light
            tmp_probe = rectpulse([0,1,2],nProbe/3); % probe, intermixed compund and 1-cs trials
            probe_cs = tmp_probe(randperm(length(tmp_probe))); % shuffle
            idx_tone = probe_cs~=2;
            cs_tone(nAcquis+1:end) = idx_tone; % set tone
            idx_light = probe_cs~=1;
            cs_light(nAcquis+1:end) = idx_light; % set tone
        end
                        
        for m = 1:2 % 1- state space; 2- Rescorla-Wagner
            
            for t = 2:nT
                
                if m==1
                    
                    x_ss(t) = A * x_ss(t-1) + B * (err * us(t-1) - x_ss(t-1));
                    
                elseif m==2
                    
                    if cs_tone(t-1) && cs_light(t-1) % update both tone & light
                        
                        V(t-1) = V_goal(t-1) + V_tone(t-1) + V_light(t-1);
                        
                        V_tone(t) = V_tone(t-1) + alpha_tone(sim) * beta * (lambda * us(t-1) - V(t-1));
                        V_light(t) = V_light(t-1) + alpha_light(sim) * beta * (lambda * us(t-1) - V(t-1));
                                                
                    elseif cs_tone(t-1) && ~cs_light(t-1) % update tone
                        V(t-1) = V_goal(t-1) + V_tone(t-1);
                        
                        V_tone(t) = V_tone(t-1) + alpha_tone(sim) * beta * (lambda * us(t-1) - V(t-1));
                        V_light(t) = V_light(t-1);
                        
                    elseif ~cs_tone(t-1) && cs_light(t-1) % update light
                        V(t-1) = V_goal(t-1) + V_light(t-1);
                        
                        V_tone(t) = V_tone(t-1);
                        V_light(t) = V_light(t-1) + alpha_light(sim) * beta * (lambda * us(t-1) - V(t-1));
                        
                    end
                    
                    V_goal(t) = V_goal(t-1) + alpha_goal * beta * (lambda * us(t-1) - V(t-1)); % V_goal is updated on every trial
                    
                end
                
            end
            
            if m==1
                X = x_ss;
            elseif m==2
                X = V;
            end
            
            % Compute changes in hand angle
            dX = [nan diff(X)];
            
            if exp==1
                itone = find(cs_tone); % tone (cs+) trials
                ilight = find(cs_light); % light (cs-) trials
                
                ditone=diff(itone);
                dilight=diff(ilight);
                
                itpt = itone(find(ditone==1)+1); % trials for tone (cs+) after tone
                ilpt = find(diff(cs_light)==1)+1;% trials for light (cs-) after tone
                itpl = find(diff(cs_tone)==1)+1;% trials for tone after light
                ilpl = ilight(find(dilight==1)+1); % trials for light after light
                
                dX_tpt = [nanmean(dX(intersect(1:nAcquis,itpt))) nanmean(dX(intersect((nAcquis+1):nT,itpt)))]; % acquisition and probe
                dX_lpt = [nanmean(dX(intersect(1:nAcquis,ilpt))) nanmean(dX(intersect((nAcquis+1):nT,ilpt)))]; % acquisition and probe
                dX_tpl = [nanmean(dX(intersect(1:nAcquis,itpl))) nanmean(dX(intersect((nAcquis+1):nT,itpl)))]; % acquisition and probe
                dX_lpl = [nanmean(dX(intersect(1:nAcquis,ilpl))) nanmean(dX(intersect((nAcquis+1):nT,ilpl)))]; % acquisition and probe
                
                dha_acq = [dX_tpt(1), dX_lpt(1), dX_tpl(1), dX_lpl(1)];
                dha_probe = [dX_tpt(2), dX_lpt(2), dX_tpl(2), dX_lpl(2)];
                
                if m==1
                    dha_acq_SS(sim,:) = dha_acq;
                    dha_probe_SS(sim,:) = dha_probe;
                elseif m==2
                    dha_acq_RW(sim,:) = dha_acq;
                    dha_probe_RW(sim,:) = dha_probe;
                else
                    dha_acq_SScntxt(sim,:) = dha_acq;
                    dha_probe_SScntxt(sim,:) = dha_probe;
                end
                
            else
                
                icomp = find(cs_tone & cs_light);
                itone = find(cs_tone & ~cs_light);
                ilight = find(~cs_tone & cs_light);
                
                isingle=union(itone,ilight);
                
                dt_n_comp=diff(icomp);
                dt_n_single=diff(isingle);
    
                cs_comp=zeros(nT,1);
                cs_comp(icomp)=1;
                cs_single=zeros(nT,1);
                cs_single(isingle)=1;
                
                t_n1_cAc=icomp(find(dt_n_comp==1)+1);
                t_n1_cAs=find(diff(cs_comp)==1)+1;
                t_n1_sAs=isingle(find(dt_n_single==1)+1);
                t_n1_sAc=find(diff(cs_single)==1)+1;
                
                t_n1_toneAc=intersect(itone,t_n1_sAc);
                t_n1_lightAc=intersect(ilight,t_n1_sAc);
    
                dX_comp = nanmean(dX(intersect((nAcquis+1):nT,icomp))); % here only interested in the probe phase
                dX_tone = nanmean(dX(itone)); % single CSs only appear in the probe phase
                dX_light = nanmean(dX(ilight)); % single CSs only appear in the probe phase
                
                dha = [dX_comp dX_tone dX_light];
                
                if m==1
                    dha_SS(sim,:) = dha;
                elseif m==2
                    dha_RW(sim,:) = dha;
                end
            end
            
        end
        
    end

end
        
for exp=[1,4]
    if exp==1
        
        if newSet_exp1
            save('simSet_Exp1_Diff','dha_acq_SS','dha_probe_SS','dha_acq_RW','dha_probe_RW')
        else
            load('simSet_Exp1_Diff')
        end
        
        plotDiffHandAngle_Differential_sim(dha_acq_RW)
        plotDiffHandAngle_Differential_sim(dha_probe_RW)
        
    else
        
        if newSet_exp4
            save('simSet_Exp4_Comp','dha_SS','dha_RW')
        else
            load('simSet_Exp4_Comp')
        end
        
        plotDiffHandAngle_CompoundToneLight_sim(dha_RW)
        
    end
end