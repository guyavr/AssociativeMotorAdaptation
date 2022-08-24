%% objective function fro fitting Rescorla-Wagner model
function [error] = func_RW(params,ha,us)

alpha_plan = params(1);
alpha_tone = params(2);
alpha_light = params(3);
beta = params(4);
lambda = params(5);

ntrials = length(us);
% separate learning from probe
cs = us;
us(601:end) = 0; % probe/washout

V = 0;
V_plan = 0;
V_tone = 0;
V_light = 0;

% trial loop
for t = 1:ntrials
    SIM(t) = V;
    if ~isnan(ha(t))
        % clamp trial?
        if us(t)
            err = lambda; % error
        else
            err = 0;
        end
        
        if cs(t)
            V = V_plan + V_tone;
            % update tone only
            delta_tone = alpha_tone * beta * (err - V);
            V_tone = V_tone + delta_tone;
        else
            V = V_plan + V_light;
            % update light only
            delta_light = alpha_light * beta * (err - V);
            V_light = V_light + delta_light;
        end
        % update plan (every trial)
        delta_plan = alpha_plan * beta * (err - V);
        V_plan = V_plan + delta_plan;
    end
end

valid = find(~isnan(ha));
error = nansum((ha(valid)-SIM(valid)).^2);
