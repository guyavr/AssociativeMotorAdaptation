%% objective function for fitting state-space model
function [error] = func_SS(params,ha,us)

A = params(1);
B = params(2);
err = params(3);

ntrials = length(us);
% clamp_size = 15;
X = zeros(1,ntrials+1);
us(601:end) = 0; % probe/washout

for t = 1:ntrials
    if ~isnan(ha(t))
        if t <= 600
            spe = err*us(t) - X(t);
        else
            spe = 0 - X(t);
        end
        X(t+1) = A*X(t) + B*spe;
    end
end

valid = find(~isnan(ha));
error = nansum((ha(valid)-X(valid)).^2);
