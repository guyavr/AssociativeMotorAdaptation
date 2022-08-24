function ana = summaryAna_differential(a_s)
% a_s- the data of all participants

nS=size(a_s,1);
nC=size(a_s,2);

% bootstrap for non-parametric confidence interval of the mean
nboot=1000;
bootfun=@(x)nanmean(x);

lilh=nan(1,nC);
lilp=nan(1,nC);

if nS>=4
    for c=1:nC
        [lilh(c),lilp(c)] = lillietest(a_s(:,c));
    end
end
            
m=nanmean(a_s);
ci=bootci(nboot,bootfun,a_s);
se=std(a_s)/sqrt(nS); % standard error of the mean

% summary statistics of main effects
a_sPrev=[nanmean(a_s(:,1:2),2) nanmean(a_s(:,3:4),2)];
a_sCurr=[nanmean(a_s(:,[1 3]),2) nanmean(a_s(:,[2 4]),2)];

mPrev=mean(a_sPrev);
mCurr=mean(a_sCurr);
sePrev=std(a_sPrev)/sqrt(nS); % standard error of the mean
seCurr=std(a_sCurr)/sqrt(nS); % standard error of the mean
ciPrev=1.96*sePrev;
ciCurr=1.96*seCurr;

% difference scores (adaptation effect)
da_sPrev=a_sPrev(:,1)-a_sPrev(:,2);
mda_sPrev=mean(da_sPrev);
seda_sPrev=std(da_sPrev)/sqrt(nS); % standard error of the mean
seda_s_rangePrev=[mda_sPrev-seda_sPrev;mda_sPrev+seda_sPrev];
cida_sPrev=1.96*seda_sPrev; % CI
cida_s_rangePrev=[mda_sPrev-cida_sPrev;mda_sPrev+cida_sPrev];

% difference scores (Pavlovian effect)
da_sCurr=a_sCurr(:,1)-a_sCurr(:,2);
mda_sCurr=mean(da_sCurr);
seda_sCurr=std(da_sCurr)/sqrt(nS); % standard error of the mean
seda_s_rangeCurr=[mda_sCurr-seda_sCurr;mda_sCurr+seda_sCurr];
cida_sCurr=1.96*seda_sCurr; % CI
cida_s_rangeCurr=[mda_sCurr-cida_sCurr;mda_sCurr+cida_sCurr];

% ANOVA
y=a_s;
tbl=table(y(:,1),y(:,2),y(:,3),y(:,4),'VariableNames',{'eAe','hAe','eAh','hAh'});
conds = table(categorical({'e' 'e' 'h' 'h'}'),categorical({'e' 'h' 'e' 'h'}'),'VariableNames',{'Prev','Curr'});
rm = fitrm(tbl,'eAe-hAh~1','WithinDesign',conds);
[ranovatbl,~,C,~ ]= ranova(rm,'WithinModel','Prev*Curr');

pPrev=table2array(ranovatbl('(Intercept):Prev','pValue'));
pCurr=table2array(ranovatbl('(Intercept):Curr','pValue'));
pPrevByCurr=table2array(ranovatbl('(Intercept):Prev:Curr','pValue'));

compPrev=multcompare(rm,'Prev');
compCurr=multcompare(rm,'Curr');
compPrevByCurr=multcompare(rm,'Prev','By','Curr','ComparisonType','bonferroni');
compCurrByPrev=multcompare(rm,'Curr','By','Prev','ComparisonType','bonferroni');

% effect size: partial eta squared
SS_Prev=table2array(ranovatbl('(Intercept):Prev','SumSq'));
SS_PrevErr=table2array(ranovatbl('Error(Prev)','SumSq'));
SS_Curr=table2array(ranovatbl('(Intercept):Curr','SumSq'));
SS_CurrErr=table2array(ranovatbl('Error(Curr)','SumSq'));
SS_PrevByCurr=table2array(ranovatbl('(Intercept):Prev:Curr','SumSq'));
SS_PrevByCurrErr=table2array(ranovatbl('Error(Prev:Curr)','SumSq'));
pEtaSq_Prev=SS_Prev/(SS_Prev+SS_PrevErr);
pEtaSq_Curr=SS_Curr/(SS_Curr+SS_CurrErr);
pEtaSq_PrevByCurr=SS_PrevByCurr/(SS_PrevByCurr+SS_PrevByCurrErr);

% Bayes Factor
F_Prev=table2array(ranovatbl('(Intercept):Prev','F'));
F_Curr=table2array(ranovatbl('(Intercept):Curr','F'));
F_PrevByCurr=table2array(ranovatbl('(Intercept):Prev:Curr','F'));
df1=table2array(ranovatbl('(Intercept)','DF'));
df2=table2array(ranovatbl('Error','DF'));
bf10_Prev = bf.bfFromF(F_Prev,df1,df2,nS);
bf10_Curr = bf.bfFromF(F_Curr,df1,df2,nS);
bf10_PrevByCurr = bf.bfFromF(F_PrevByCurr,df1,df2,nS);

ana.indiv=a_s;
ana.m=m;
ana.ci=ci;
ana.se=se;
ana.mPrev=mPrev;
ana.sePrev=sePrev;
ana.ciPrev=ciPrev;
ana.mCurr=mCurr;
ana.seCurr=seCurr;
ana.ciCurr=ciCurr;

ana.da_sPrev=da_sPrev;
ana.mda_sPrev=mda_sPrev;
ana.seda_s_rangePrev=seda_s_rangePrev;
ana.cida_s_rangePrev=cida_s_rangePrev;
ana.da_sCurr=da_sCurr;
ana.mda_sCurr=mda_sCurr;
ana.seda_s_rangeCurr=seda_s_rangeCurr;
ana.cida_s_rangeCurr=cida_s_rangeCurr;

ana.ranovatbl=ranovatbl;
ana.pPrev=pPrev;
ana.pCurr=pCurr;
ana.pPrevByCurr=pPrevByCurr;
ana.compPrevByCurr=compPrevByCurr;
ana.compCurrByPrev=compCurrByPrev;
ana.pEtaSq_Prev=pEtaSq_Prev;
ana.pEtaSq_Curr=pEtaSq_Curr;
ana.pEtaSq_PrevByCurr=pEtaSq_PrevByCurr;
ana.lilh=lilh;
ana.lilp=lilp;

ana.bf10_Prev=bf10_Prev;
ana.bf10_Curr=bf10_Curr;
ana.bf10_PrevByCurr=bf10_PrevByCurr;

end

