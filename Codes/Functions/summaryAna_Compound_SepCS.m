function ana = summaryAna_Compound_SepCS(a_s)
% a- the data of all participants

nS=size(a_s,1);
nC=size(a_s,2);

% bootstrap for non-parametric confidence interval of the mean
nboot=5000;
bootfun=@(x)nanmean(x);

lilh=nan(1,nC);
lilp=nan(1,nC);

if nS>=4
    for c=1:nC
        [lilh(c),lilp(c)] = lillietest(a_s(:,c));
    end
end

m=nanmean(a_s);
se=std(a_s)/sqrt(nS); % standard error of the mean

med=nanmedian(a_s);
ci=bootci(nboot,bootfun,a_s);
interqr=iqr(a_s);
prctile25_75=prctile(a_s,[25 75],1);

% Since the the data of from the compound condition does not pass normality test,
% use nonparametric Friedman's test
[pCond,tbl,stats] = friedman(a_s);
mltcmp=multcompare(stats,'ctype','bonferroni');
X2=tbl{2, 5};
kendallsW=X2/(nS*(nC-1)); %Kendall's W Value effect size W = X2/N(K-1)

% correlation
a_tl=a_s(:,2:3);
[rho_pear,pval_pear] = corr(a_tl);
corrPear_r = rho_pear(1,2);
corrPear_p = pval_pear(1,2);
[corrPear_bf10,~,~] = bf.corr(a_tl(:,1),a_tl(:,2));

[rho_spear,pval_spear] = corr(a_tl,'Type','Spearman');
corrSpear_r = rho_spear(1,2);
corrSpear_p = pval_spear(1,2);

ana.indiv=a_s;
ana.m=m;
ana.med=med;
ana.ci=ci;
ana.se=se;
ana.interqr=interqr;
ana.prctile25_75=prctile25_75;
ana.friedman_tbl=tbl;
ana.friedman_stats=stats;
ana.kendallsW=kendallsW;
ana.pCond=pCond;
ana.mltcmp=mltcmp;
ana.lilh=lilh;
ana.lilp=lilp;

ana.corrPear_r=corrPear_r;
ana.corrPear_p=corrPear_p;
ana.corrPear_bf10=corrPear_bf10;
ana.corrSpear_r=corrSpear_r;
ana.corrSpear_p=corrSpear_p;

end

