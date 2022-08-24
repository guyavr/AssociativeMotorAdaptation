function mea_mC_mB = mean_cond_compound(mea,tr_c)
% return the mean of each condition in each block.
% mea- data
% tr_c- cell array of trials of each condition

nC=size(tr_c,2); % number of conditions;
nS=size(mea,1); % number of participants

mea_mC_mB=nan(nS,nC);

for s=1:nS
    for c=1:nC
        mea_mC_mB(s,c)=nanmean(mea(s,tr_c{s,c}));
    end
end

end

