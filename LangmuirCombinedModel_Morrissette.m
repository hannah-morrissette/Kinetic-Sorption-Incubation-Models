% %% FULL DATA SET
% 
% load('JBF fit')
% jbtest=fitnlm(tin,JBF,@fitsat,[2.64 0.77 4.2]);
% 
% figure()
% scatter(tin(1:24),JBF(1:24))
% hold on
% scatter(tin(1:24), jbtest.Fitted(1:24))
% title('JBF low - full data set')
% xlabel('Time (hrs)')
% ylabel('DOCw (mg)')
% legend('obs','model')
% 
% figure()
% scatter(tin(25:48), JBF(25:48))
% hold on
% scatter(tin(25:48), jbtest.Fitted(25:48))
% title('JBF high - full data set')
% xlabel('Time (hrs)')
% ylabel('DOCw (mg)')
% legend('obs','model')

%% data
global y0z qm modout

char = xlsread('qmax.xlsx');
qmax = char(:,1);
y0all = xlsread('y0_analytical.xlsx');
rates = xlsread('analytical_output.xlsx');

rownames = {'kads';'kdes';'adjRsq';'AIC';'AICc';'RMSE'};
colnames = {'JBF' 'JBS' 'TAF' 'TAS' 'WCF' 'WCS' 'WIF' 'WIS' 'WMF' 'WMS' 'PCF' 'PCS' 'PIF' 'PIS' 'PMF' 'PMS'};

%% REDUCED DATA SET - JBF
JBFred = xlsread('obs_analytical_JB.xlsx');
JBFred2 = JBFred(1:30,1);
tinred = JBFred(1:30,2);
y0z = [[y0all(3,2) y0all(1,2)];[y0all(3,1) y0all(1,1)]] ;
beta = [rates(1,1) rates(2,1)];
qm = qmax(1);
jbftest=fitnlm(tinred,JBFred2,@fitsatredfixS0,beta);
jbfmodout = modout;

figure()
scatter(tinred(1:15),JBFred2(1:15))
hold on
scatter(tinred(1:15), jbftest.Fitted(1:15))
hold on
scatter(tinred(16:30), JBFred2(16:30))
hold on
scatter(tinred(16:30), jbftest.Fitted(16:30))
title('JBF saturated')
xlabel('Time (hrs)')
ylabel('DOCw (mg)')
legend('low obs','low model','high obs','high model')
saveas(gcf,['JBFsat'],'png');

%% JBS
JBSred = xlsread('obs_analytical_JB.xlsx');
JBSred2 = JBSred(31:60,1);
tinred = JBSred(31:60,2);
y0z = [[y0all(4,2) y0all(2,2)];[y0all(4,1) y0all(2,1)]];
beta = [rates(1,2) rates(2,2)];
qm = qmax(2);
jbstest=fitnlm(tinred,JBSred2,@fitsatredfixS0,beta);
jbsmodout = modout;

figure()
scatter(tinred(1:15),JBSred2(1:15))
hold on
scatter(tinred(1:15), jbstest.Fitted(1:15))
hold on
scatter(tinred(16:30), JBSred2(16:30))
hold on
scatter(tinred(16:30), jbstest.Fitted(16:30))
title('JBS saturated')
xlabel('Time (hrs)')
ylabel('DOCw (mg)')
legend('low obs','low model','high obs','high model')
saveas(gcf,['JBSsat'],'png');
%% REDUCED DATA SET - TAF
TAFred = xlsread('obs_analytical_TA.xlsx');
TAFred2 = TAFred(1:30,1);
tinred = TAFred(1:30,2);
y0z = [[y0all(19,2) y0all(17,2)];[y0all(19,1) y0all(17,1)]] ;
beta = [rates(1,3) rates(2,3)];
qm = qmax(5);

taftest=fitnlm(tinred,TAFred2,@fitsatredfixS0,beta);
tafmodout = modout;

figure()
scatter(tinred(1:15),TAFred2(1:15))
hold on
scatter(tinred(1:15), taftest.Fitted(1:15))
hold on
scatter(tinred(16:30), TAFred2(16:30))
hold on
scatter(tinred(16:30), taftest.Fitted(16:30))
title('TAF saturated')
xlabel('Time (hrs)')
ylabel('DOCw (mg)')
legend('low obs','low model','high obs','high model')
saveas(gcf,['TAFsat'],'png');

%% TAS
TASred = xlsread('obs_analytical_TA.xlsx');
TASred2 = TASred(31:60,1);
tinred = TASred(31:60,2);
y0z = [[y0all(20,2) y0all(18,2)];[y0all(20,1) y0all(18,1)]];
beta = [0.1 rates(2,4)];
qm = qmax(6);

tastest=fitnlm(tinred,TASred2,@fitsatredfixS0,beta);
tasmodout = modout;

figure()
scatter(tinred(1:15),TASred2(1:15))
hold on
scatter(tinred(1:15), tastest.Fitted(1:15))
hold on
scatter(tinred(16:30), TASred2(16:30))
hold on
scatter(tinred(16:30), tastest.Fitted(16:30))
title('TAS saturated')
xlabel('Time (hrs)')
ylabel('DOCw (mg)')
legend('low obs','low model','high obs','high model')
saveas(gcf,['TASsat'],'png');

%% REDUCED DATA SET - WCF
WCFred = xlsread('obs_analytical_WC.xlsx');
WCFred2 = WCFred(1:30,1);
tinred = WCFred(1:30,2);
y0z = [[y0all(23,2) y0all(21,2)];[y0all(23,1) y0all(21,1)]];
%beta = [rates(1,5) rates(2,5)];
beta = [0.01 0.01];
qm = qmax(11);

wcftest=fitnlm(tinred,WCFred2,@fitsatredfixS0,beta);
wcfmodout = modout;

figure()
scatter(tinred(1:15),WCFred2(1:15))
hold on
scatter(tinred(1:15), wcftest.Fitted(1:15))
hold on
scatter(tinred(16:30), WCFred2(16:30))
hold on
scatter(tinred(16:30), wcftest.Fitted(16:30))
title('WCF saturated')
xlabel('Time (hrs)')
ylabel('DOCw (mg)')
legend('low obs','low model','high obs','high model')
saveas(gcf,['WCFsat'],'png');

%% WCS
WCSred = xlsread('obs_analytical_WC.xlsx');
WCSred2 = WCSred(31:57,1);
tinred = WCSred(31:57,2);
y0z = [[y0all(24,2) y0all(22,2)];[y0all(24,1) y0all(22,1)]];
%beta = [rates(1,6) rates(2,6)];
beta = [0.1 0.1];
qm = qmax(10);

wcstest=fitnlm(tinred,WCSred2,@fitsatredfixS0,beta);
wcsmodout = modout;

figure()
scatter(tinred(1:12),WCSred2(1:12))
hold on
scatter(tinred(1:12), wcstest.Fitted(1:12))
hold on
scatter(tinred(13:27), WCSred2(13:27))
hold on
scatter(tinred(13:27), wcstest.Fitted(13:27))
title('WCS saturated')
xlabel('Time (hrs)')
ylabel('DOCw (mg)')
legend('low obs','low model','high obs','high model')
saveas(gcf,['WCSsat'],'png');

%% REDUCED DATA SET - WIF
WIFred = xlsread('obs_analytical_WI.xlsx');
WIFred2 = WIFred(1:33,1);
tinred = WIFred(1:33,2);
y0z = [[y0all(27,2) y0all(25,2)];[y0all(27,1) y0all(25,1)]];
%beta = [rates(1,7) rates(2,7)];
beta = [0.1 0.1];
qm = qmax(13);

wiftest=fitnlm(tinred,WIFred2,@fitsatredfixS0,beta);
wifmodout = modout;

figure()
scatter(tinred(1:15),WIFred2(1:15))
hold on
scatter(tinred(1:15), wiftest.Fitted(1:15))
hold on
scatter(tinred(16:33), WIFred2(16:33))
hold on
scatter(tinred(16:33), wiftest.Fitted(16:33))
title('WIF saturated')
xlabel('Time (hrs)')
ylabel('DOCw (mg)')
legend('low obs','low model','high obs','high model')
saveas(gcf,['WIFsat'],'png');

%% WIS
WISred = xlsread('obs_analytical_WI.xlsx');
WISred2 = WISred(34:63,1);
tinred = WISred(34:63,2);
y0z = [[y0all(28,2) y0all(26,2)];[y0all(28,1) y0all(26,1)]];
%beta = [rates(1,8) rates(2,8)];
beta = [0.1 0.1];
qm = qmax(14);

wistest=fitnlm(tinred,WISred2,@fitsatredfixS0,beta);
wismodout = modout;

figure()
scatter(tinred(1:15),WISred2(1:15))
hold on
scatter(tinred(1:15), wistest.Fitted(1:15))
hold on
scatter(tinred(16:30), WISred2(16:30))
hold on
scatter(tinred(16:30), wistest.Fitted(16:30))
title('WIS saturated')
xlabel('Time (hrs)')
ylabel('DOCw (mg)')
legend('low obs','low model','high obs','high model')
saveas(gcf,['WISsat'],'png');

%% REDUCED DATA SET - WMF
WMFred = xlsread('obs_analytical_WM.xlsx');
WMFred2 = WMFred(1:30,1);
tinred = WMFred(1:30,2);
y0z = [[y0all(31,2) y0all(29,2)];[y0all(31,1) y0all(29,1)]];
beta = [rates(1,9) rates(2,9)];
qm = qmax(17);

wmftest=fitnlm(tinred,WMFred2,@fitsatredfixS0,beta);
wmfmodout = modout;

figure()
scatter(tinred(1:15),WMFred2(1:15))
hold on
scatter(tinred(1:15), wmftest.Fitted(1:15))
hold on
scatter(tinred(16:30), WMFred2(16:30))
hold on
scatter(tinred(16:30), wmftest.Fitted(16:30))
title('WMF saturated')
xlabel('Time (hrs)')
ylabel('DOCw (mg)')
legend('low obs','low model','high obs','high model')
saveas(gcf,['WMFsat'],'png');

%% WMS
WMSred = xlsread('obs_analytical_WM.xlsx');
WMSred2 = WMSred(31:60,1);
tinred = WMSred(31:60,2);
y0z = [[y0all(32,2) y0all(30,2)];[y0all(32,1) y0all(30,1)]];
%beta = [rates(1,10) rates(2,10)];
beta = [0.1 0.1];
qm = qmax(18);

wmstest=fitnlm(tinred,WMSred2,@fitsatredfixS0,beta);
wmsmodout = modout;

figure()
scatter(tinred(1:15),WMSred2(1:15))
hold on
scatter(tinred(1:15), wmstest.Fitted(1:15))
hold on
scatter(tinred(16:30), WMSred2(16:30))
hold on
scatter(tinred(16:30), wmstest.Fitted(16:30))
title('WMS saturated')
xlabel('Time (hrs)')
ylabel('DOCw (mg)')
legend('low obs','low model','high obs','high model')
saveas(gcf,['WMSsat'],'png');

%% REDUCED DATA SET - PCF
PCFred = xlsread('obs_analytical_PC.xlsx');
PCFred2 = PCFred(1:30,1);
tinred = PCFred(1:30,2);
y0z = [[y0all(7,2) y0all(5,2)];[y0all(7,1) y0all(5,1)]];
%beta = [rates(1,11) rates(2,11)];
beta = [1 1];
qm = qmax(21);

pcftest=fitnlm(tinred,PCFred2,@fitsatredfixS0,beta);
pcfmodout = modout;

figure()
scatter(tinred(1:15),PCFred2(1:15))
hold on
scatter(tinred(1:15), pcftest.Fitted(1:15))
hold on
scatter(tinred(16:30), PCFred2(16:30))
hold on
scatter(tinred(16:30), pcftest.Fitted(16:30))
title('PCF saturated')
xlabel('Time (hrs)')
ylabel('DOCw (mg)')
legend('low obs','low model','high obs','high model')
saveas(gcf,['PCFsat'],'png');

%% PCS
PCSred = xlsread('obs_analytical_PC.xlsx');
PCSred2 = PCSred(31:60,1);
tinred = PCSred(31:60,2);
y0z = [[y0all(8,2) y0all(6,2)];[y0all(8,1) y0all(6,1)]];
beta = [rates(1,12) rates(2,12)];
qm = qmax(22);

pcstest=fitnlm(tinred,PCSred2,@fitsatredfixS0,beta);
pcsmodout = modout;

figure()
scatter(tinred(1:15),PCSred2(1:15))
hold on
scatter(tinred(1:15), pcstest.Fitted(1:15))
hold on
scatter(tinred(16:30), PCSred2(16:30))
hold on
scatter(tinred(16:30), pcstest.Fitted(16:30))
title('PCS saturated')
xlabel('Time (hrs)')
ylabel('DOCw (mg)')
legend('low obs','low model','high obs','high model')
saveas(gcf,['PCSsat'],'png');

%% REDUCED DATA SET - PIF
PIFred = xlsread('obs_analytical_PI.xlsx');
PIFred2 = PIFred(1:33,1);
tinred = PIFred(1:33,2);
y0z = [[y0all(11,2) y0all(9,2)];[y0all(11,1) y0all(9,1)]];
beta = [rates(1,13) rates(2,13)];
qm = qmax(25);

piftest=fitnlm(tinred,PIFred2,@fitsatredfixS0,beta);
pifmodout = modout;

figure()
scatter(tinred(1:15),PIFred2(1:15))
hold on
scatter(tinred(1:15), piftest.Fitted(1:15))
hold on
scatter(tinred(16:33), PIFred2(16:33))
hold on
scatter(tinred(16:33), piftest.Fitted(16:33))
title('PIF saturated')
xlabel('Time (hrs)')
ylabel('DOCw (mg)')
legend('low obs','low model','high obs','high model')
saveas(gcf,['PIFsat'],'png');

%% PIS
PISred = xlsread('obs_analytical_PI.xlsx');
PISred2 = PISred(34:63,1);
tinred = PISred(34:63,2);
y0z = [[y0all(12,2) y0all(10,2)];[y0all(12,1) y0all(10,1)]];
%beta = [rates(1,14) rates(2,14)];
beta = [0.1 0.1];
qm = qmax(26);

pistest=fitnlm(tinred,PISred2,@fitsatredfixS0,beta);
pismodout = modout;

figure()
scatter(tinred(1:15),PISred2(1:15))
hold on
scatter(tinred(1:15), pistest.Fitted(1:15))
hold on
scatter(tinred(16:30), PISred2(16:30))
hold on
scatter(tinred(16:30), pistest.Fitted(16:30))
title('PIS saturated')
xlabel('Time (hrs)')
ylabel('DOCw (mg)')
legend('low obs','low model','high obs','high model')
saveas(gcf,['PISsat'],'png');

%% REDUCED DATA SET - PMF
PMFred = xlsread('obs_analytical_PM.xlsx');
PMFred2 = PMFred(1:36,1);
tinred = PMFred(1:36,2);
y0z = [[y0all(15,2) y0all(13,2)];[y0all(15,1) y0all(13,1)]];
beta = [rates(1,15) rates(2,15)];
qm = qmax(29);

pmftest=fitnlm(tinred,PMFred2,@fitsatredfixS0,beta);
pmfmodout = modout;

figure()
scatter(tinred(1:15),PMFred2(1:15))
hold on
scatter(tinred(1:15), pmftest.Fitted(1:15))
hold on
scatter(tinred(16:36), PMFred2(16:36))
hold on
scatter(tinred(16:36), pmftest.Fitted(16:36))
title('PMF saturated')
xlabel('Time (hrs)')
ylabel('DOCw (mg)')
legend('low obs','low model','high obs','high model')
saveas(gcf,['PMFsat'],'png');

%% PMS
PMSred = xlsread('obs_analytical_PM.xlsx');
PMSred2 = PMSred(37:66,1);
tinred = PMSred(37:66,2);
y0z = [[y0all(16,2) y0all(14,2)];[y0all(16,1) y0all(14,1)]];
%beta = [rates(1,16) rates(2,16)];
beta = [0.1 0.1];
qm = qmax(30);

pmstest=fitnlm(tinred,PMSred2,@fitsatredfixS0,beta);
pmsmodout = modout;

figure()
scatter(tinred(1:15),PMSred2(1:15))
hold on
scatter(tinred(1:15), pmstest.Fitted(1:15))
hold on
scatter(tinred(16:30), PMSred2(16:30))
hold on
scatter(tinred(16:30), pmstest.Fitted(16:30))
title('PMS saturated')
xlabel('Time (hrs)')
ylabel('DOCw (mg)')
legend('low obs','low model','high obs','high model')
saveas(gcf,['PMSsat'],'png');

%% Data Output

jbfstat = [jbftest.Rsquared.Adjusted(1,1) jbftest.ModelCriterion.AIC(1,1) jbftest.ModelCriterion.AICc(1,1) jbftest.RMSE(1,1)]; 
jbfstat2 = renamevars(array2table(jbfstat'),["Var1"],["Estimate"]);
jbfall = renamevars([jbftest.Coefficients(:,1); jbfstat2],["Estimate"],["JBF"]);
jbfall = [rownames jbfall];

jbsstat = [jbstest.Rsquared.Adjusted(1,1) jbstest.ModelCriterion.AIC(1,1) jbstest.ModelCriterion.AICc(1,1) jbstest.RMSE(1,1)]; 
jbsstat2 = renamevars(array2table(jbsstat'),["Var1"],["Estimate"]);
jbsall = renamevars([jbstest.Coefficients(:,1); jbsstat2],["Estimate"],["JBS"]);
jball = [jbfall jbsall];

tafstat = [taftest.Rsquared.Adjusted(1,1) taftest.ModelCriterion.AIC(1,1) taftest.ModelCriterion.AICc(1,1) taftest.RMSE(1,1)]; 
tafstat2 = renamevars(array2table(tafstat'),["Var1"],["Estimate"]);
tafall = renamevars([taftest.Coefficients(:,1); tafstat2],["Estimate"],["TAF"]);

tasstat = [tastest.Rsquared.Adjusted(1,1) tastest.ModelCriterion.AIC(1,1) tastest.ModelCriterion.AICc(1,1) tastest.RMSE(1,1)]; 
tasstat2 = renamevars(array2table(tasstat'),["Var1"],["Estimate"]);
tasall = renamevars([tastest.Coefficients(:,1); tasstat2],["Estimate"],["TAS"]);
taall = [tafall tasall];

wcfstat = [wcftest.Rsquared.Adjusted(1,1) wcftest.ModelCriterion.AIC(1,1) wcftest.ModelCriterion.AICc(1,1) wcftest.RMSE(1,1)]; 
wcfstat2 = renamevars(array2table(wcfstat'),["Var1"],["Estimate"]);
wcfall = renamevars([wcftest.Coefficients(:,1); wcfstat2],["Estimate"],["WCF"]);

wcsstat = [wcstest.Rsquared.Adjusted(1,1) wcstest.ModelCriterion.AIC(1,1) wcstest.ModelCriterion.AICc(1,1) wcstest.RMSE(1,1)]; 
wcsstat2 = renamevars(array2table(wcsstat'),["Var1"],["Estimate"]);
wcsall = renamevars([wcstest.Coefficients(:,1); wcsstat2],["Estimate"],["WCS"]);
wcall = [wcfall wcsall];

wifstat = [wiftest.Rsquared.Adjusted(1,1) wiftest.ModelCriterion.AIC(1,1) wiftest.ModelCriterion.AICc(1,1) wiftest.RMSE(1,1)]; 
wifstat2 = renamevars(array2table(wifstat'),["Var1"],["Estimate"]);
wifall = renamevars([wiftest.Coefficients(:,1); wifstat2],["Estimate"],["WIF"]);

wisstat = [wistest.Rsquared.Adjusted(1,1) wistest.ModelCriterion.AIC(1,1) wistest.ModelCriterion.AICc(1,1) wistest.RMSE(1,1)]; 
wisstat2 = renamevars(array2table(wisstat'),["Var1"],["Estimate"]);
wisall = renamevars([wistest.Coefficients(:,1); wisstat2],["Estimate"],["WIS"]);
wiall = [wifall wisall];

wmfstat = [wmftest.Rsquared.Adjusted(1,1) wmftest.ModelCriterion.AIC(1,1) wmftest.ModelCriterion.AICc(1,1) wmftest.RMSE(1,1)]; 
wmfstat2 = renamevars(array2table(wmfstat'),["Var1"],["Estimate"]);
wmfall = renamevars([wmftest.Coefficients(:,1); wmfstat2],["Estimate"],["WMF"]);

wmsstat = [wmstest.Rsquared.Adjusted(1,1) wmstest.ModelCriterion.AIC(1,1) wmstest.ModelCriterion.AICc(1,1) wmstest.RMSE(1,1)]; 
wmsstat2 = renamevars(array2table(wmsstat'),["Var1"],["Estimate"]);
wmsall = renamevars([wmstest.Coefficients(:,1); wmsstat2],["Estimate"],["WMS"]);
wmall = [wmfall wmsall];

pcfstat = [pcftest.Rsquared.Adjusted(1,1) pcftest.ModelCriterion.AIC(1,1) pcftest.ModelCriterion.AICc(1,1) pcftest.RMSE(1,1)]; 
pcfstat2 = renamevars(array2table(pcfstat'),["Var1"],["Estimate"]);
pcfall = renamevars([pcftest.Coefficients(:,1); pcfstat2],["Estimate"],["PCF"]);

pcsstat = [pcstest.Rsquared.Adjusted(1,1) pcstest.ModelCriterion.AIC(1,1) pcstest.ModelCriterion.AICc(1,1) pcstest.RMSE(1,1)]; 
pcsstat2 = renamevars(array2table(pcsstat'),["Var1"],["Estimate"]);
pcsall = renamevars([pcstest.Coefficients(:,1); pcsstat2],["Estimate"],["PCS"]);
pcall = [pcfall pcsall];

pifstat = [piftest.Rsquared.Adjusted(1,1) piftest.ModelCriterion.AIC(1,1) piftest.ModelCriterion.AICc(1,1) piftest.RMSE(1,1)]; 
pifstat2 = renamevars(array2table(pifstat'),["Var1"],["Estimate"]);
pifall = renamevars([piftest.Coefficients(:,1); pifstat2],["Estimate"],["PIF"]);

pisstat = [pistest.Rsquared.Adjusted(1,1) pistest.ModelCriterion.AIC(1,1) pistest.ModelCriterion.AICc(1,1) pistest.RMSE(1,1)]; 
pisstat2 = renamevars(array2table(pisstat'),["Var1"],["Estimate"]);
pisall = renamevars([pistest.Coefficients(:,1); pisstat2],["Estimate"],["PIS"]);
piall = [pifall pisall];

pmfstat = [pmftest.Rsquared.Adjusted(1,1) pmftest.ModelCriterion.AIC(1,1) pmftest.ModelCriterion.AICc(1,1) pmftest.RMSE(1,1)]; 
pmfstat2 = renamevars(array2table(pmfstat'),["Var1"],["Estimate"]);
pmfall = renamevars([pmftest.Coefficients(:,1); pmfstat2],["Estimate"],["PMF"]);

pmsstat = [pmstest.Rsquared.Adjusted(1,1) pmstest.ModelCriterion.AIC(1,1) pmstest.ModelCriterion.AICc(1,1) pmstest.RMSE(1,1)]; 
pmsstat2 = renamevars(array2table(pmsstat'),["Var1"],["Estimate"]);
pmsall = renamevars([pmstest.Coefficients(:,1); pmsstat2],["Estimate"],["PMS"]);
pmall = [pmfall pmsall];

satout = [jball taall wcall wiall wmall pcall piall pmall];
writetable(satout,'newSaturatedOutput.xlsx');

%% Model Output

modall = array2table([jbfmodout jbsmodout tafmodout tasmodout wcfmodout wcsmodout wifmodout wismodout wmfmodout wmsmodout pcfmodout pcsmodout pifmodout pismodout pmfmodout pmsmodout]);
modall.Properties.VariableNames = {'JBF' 'JBS' 'TAF' 'TAS' 'WCF' 'WCS' 'WIF' 'WIS' 'WMF' 'WMS' 'PCF' 'PCS' 'PIF' 'PIS' 'PMF' 'PMS'};

writetable(modall,'newSaturatedModel.xlsx');
