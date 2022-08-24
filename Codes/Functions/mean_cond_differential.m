function mea_mC_mB = mean_cond_differential(mea,tr_c,tr_b)
% return the mean of each condition in each block.
% tr_c- cell array of trials of each condition
% tr_b- trials of the block

nC=size(tr_c,2); % number of conditions;
nS=size(mea,1); % number of participants

mea_mC_mB=nan(nS,nC);

for s=1:nS
    for c=1:nC
        tr=intersect(tr_c{s,c},tr_b);
        mea_mC_mB(s,c)=nanmean(mea(s,tr));
    end
end

end

