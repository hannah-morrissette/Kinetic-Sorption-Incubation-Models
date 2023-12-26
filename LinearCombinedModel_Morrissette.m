%% Code to fit a set of observational DOC (mg) data over time

%%INPUT DATA
global y0all y0F y0S tmod2
% DOC mass on the sediments and in the water at time 0
y0all = xlsread('y0_analytical.xlsx');
% Produced kdes and kads rates from previous run
rates = xlsread('newAnalyticalOutput.xlsx');
kads = rates(1,:);
kdes = rates(2,:);

rownames = {'kads';'kdes';'adjRsq';'AIC';'AICc';'RMSE'};
colnames = {'JBF' 'JBS' 'TAF' 'TAS' 'WCF' 'WCS' 'WIF' 'WIS' 'WMF' 'WMS' 'PCF' 'PCS' 'PIF' 'PIS' 'PMF' 'PMS'};

tmod = 0:0.1:24;
tmod1 = tmod';
tmod2 = [tmod1; tmod1];
%% Equations
% z = exp(-1*(ka + kd)*t);
% S(t) = (1 / (ka + kd))*((S0 * (ka + (kd * z)) + (W0 * ka * (1 - z))));
% W(t) = (1 / (ka + kd))*((S0 * kd * (1 - z)) + (W0 * (ka + (kd * z))));

%% Fitting Scheme

% modeloutput = fitnlm(time vector (x), observational values, @function, [beta1 beta2])

%% JB
% Observational values over time (ordered L's in triplicate then H's)
obs = xlsread('obs_analytical_JB.xlsx');
obsF = obs(1:30,:);
obsS = obs(31:60,:);

%y0 = [LS0 HS0 LW0 HW0]
y0F = [[y0all(3,2) y0all(1,2)];[y0all(3,1) y0all(1,1)]];
beta = [kads(1) kdes(1)];
jbFsol = fitnlm(obsF(:,2), obsF(:,1), @fun_analytical_FfixS0, beta);
jbFest = renamevars(jbFsol.Coefficients(:,1),["Estimate"],["JBF rates fix S0"]);
%model output
%fresh
    ka = beta(1);
    kd = beta(2);
    LS0 = y0F(1,1);  %low S0
    HS0 = y0F(1,2);  %high S0
    LW0 = y0F(2,1);  %low w0
    HW0 = y0F(2,2);  %high w0
    x = tmod2;
    z = exp(-(ka + kd)*x);
    y = zeros(482,1);
    for i = 1:241
        y(i) = (1 / (ka + kd))*((LS0 * kd * (1 - z(i))) + (LW0 * (ka + (kd * z(i)))));
    end
    for i = 242:482
        y(i) = (1 / (ka + kd))*((HS0 * kd * (1 - z(i))) + (HW0 * (ka + (kd * z(i)))));
    end
jbfmod = y;

y0S = [[y0all(4,2) y0all(2,2)];[y0all(4,1) y0all(2,1)]];
beta = [kads(2) kdes(2)];
jbSsol = fitnlm(obsS(:,2), obsS(:,1), @fun_analytical_SfixS0, beta);
jbSest = renamevars(jbSsol.Coefficients(:,1),["Estimate"],["JBS rates fix S0"]);

%model output
%saline
    ka = beta(1);
    kd = beta(2);
    LS0 = y0S(1,1);  %low S0
    HS0 = y0S(1,2);  %high S0
    LW0 = y0S(2,1);  %low w0
    HW0 = y0S(2,2);  %high w0
    x = tmod2;
    z = exp(-(ka + kd)*x);
    y = zeros(482,1);
    for i = 1:241
        y(i) = (1 / (ka + kd))*((LS0 * kd * (1 - z(i))) + (LW0 * (ka + (kd * z(i)))));
    end
    for i = 242:482
        y(i) = (1 / (ka + kd))*((HS0 * kd * (1 - z(i))) + (HW0 * (ka + (kd * z(i)))));
    end
jbsmod = y;
jbmod = [jbfmod jbsmod];

%estimate output
JBest = [jbFest jbSest];
%jbfall = jbFsol.Coefficients(:,1); 
jbfstat = [jbFsol.Rsquared.Adjusted(1,1) jbFsol.ModelCriterion.AIC(1,1) jbFsol.ModelCriterion.AICc(1,1) jbFsol.RMSE(1,1)]; 
jbfstat2 = renamevars(array2table(jbfstat'),["Var1"],["Estimate"]);
jbfall = renamevars([jbFsol.Coefficients(:,1); jbfstat2],["Estimate"],["JBF"]);
jbfall = [rownames jbfall];

jbsstat = [jbSsol.Rsquared.Adjusted(1,1) jbSsol.ModelCriterion.AIC(1,1) jbSsol.ModelCriterion.AICc(1,1) jbSsol.RMSE(1,1)]; 
jbsstat2 = renamevars(array2table(jbsstat'),["Var1"],["Estimate"]);
jbsall = renamevars([jbSsol.Coefficients(:,1); jbsstat2],["Estimate"],["JBS"]);
jball = [jbfall jbsall];

figure()
scatter(obsF(1:15,2),obsF(1:15,1));
hold on
scatter(obsF(1:15,2),jbFsol.Fitted(1:15));
hold on
scatter(obsF(16:30,2),obsF(16:30,1));
hold on
scatter(obsF(16:30,2),jbFsol.Fitted(16:30));
title('JBF Analytical')
xlabel('Time (hrs)')
ylabel('DOCw (mg)')
legend('low obs','low model','high obs','high model')
saveas(gcf,['JBFan'],'png');
figure()
scatter(obsS(1:15,2),obsS(1:15,1));
hold on
scatter(obsS(1:15,2),jbSsol.Fitted(1:15));
hold on
scatter(obsS(16:30,2),obsS(16:30,1));
hold on
scatter(obsS(16:30,2),jbSsol.Fitted(16:30));
title('JBS Analytical')
xlabel('Time (hrs)')
ylabel('DOCw (mg)')
legend('low obs','low model','high obs','high model')
saveas(gcf,['JBSan'],'png');

%% TA
% Observational values over time (ordered L's in triplicate then H's)
obs = xlsread('obs_analytical_TA.xlsx');
obsF = obs(1:30,:);
obsS = obs(31:60,:);

%y0 = [LS0 HS0 LW0 HW0]
y0F = [[y0all(19,2) y0all(17,2)];[y0all(19,1) y0all(17,1)]] ;
beta = [kads(3) kdes(3)];
%beta = [1.23 0.38 6.07 10.87];
taFsol = fitnlm(obsF(:,2), obsF(:,1), @fun_analytical_FfixS0, beta);
taFest = renamevars(taFsol.Coefficients(:,1),["Estimate"],["TAF rates"]);

%model output
%fresh
    ka = beta(1);
    kd = beta(2);
    LS0 = y0F(1,1);  %low S0
    HS0 = y0F(1,2);  %high S0
    LW0 = y0F(2,1);  %low w0
    HW0 = y0F(2,2);  %high w0
    x = tmod2;
    z = exp(-(ka + kd)*x);
    y = zeros(482,1);
    for i = 1:241
        y(i) = (1 / (ka + kd))*((LS0 * kd * (1 - z(i))) + (LW0 * (ka + (kd * z(i)))));
    end
    for i = 242:482
        y(i) = (1 / (ka + kd))*((HS0 * kd * (1 - z(i))) + (HW0 * (ka + (kd * z(i)))));
    end
tafmod = y;

y0S = [[y0all(20,2) y0all(18,2)];[y0all(20,1) y0all(18,1)]];
beta = [kads(4) kdes(4)];
%beta = [1.56 0.87 2.78 4.07];
taSsol = fitnlm(obsS(:,2), obsS(:,1), @fun_analytical_SfixS0, beta);
taSest = renamevars(taSsol.Coefficients(:,1),["Estimate"],["TAS rates"]);

%model output
%saline
    ka = beta(1);
    kd = beta(2);
    LS0 = y0S(1,1);  %low S0
    HS0 = y0S(1,2);  %high S0
    LW0 = y0S(2,1);  %low w0
    HW0 = y0S(2,2);  %high w0
    x = tmod2;
    z = exp(-(ka + kd)*x);
    y = zeros(482,1);
    for i = 1:241
        y(i) = (1 / (ka + kd))*((LS0 * kd * (1 - z(i))) + (LW0 * (ka + (kd * z(i)))));
    end
    for i = 242:482
        y(i) = (1 / (ka + kd))*((HS0 * kd * (1 - z(i))) + (HW0 * (ka + (kd * z(i)))));
    end
tasmod = y;
tamod = [tafmod tasmod];

%estimate output
TAest = [taFest taSest];
tafstat = [taFsol.Rsquared.Adjusted(1,1) taFsol.ModelCriterion.AIC(1,1) taFsol.ModelCriterion.AICc(1,1) taFsol.RMSE(1,1)]; 
tafstat2 = renamevars(array2table(tafstat'),["Var1"],["Estimate"]);
tafall = renamevars([taFsol.Coefficients(:,1); tafstat2],["Estimate"],["TAF"]);

tasstat = [taSsol.Rsquared.Adjusted(1,1) taSsol.ModelCriterion.AIC(1,1) taSsol.ModelCriterion.AICc(1,1) taSsol.RMSE(1,1)]; 
tasstat2 = renamevars(array2table(tasstat'),["Var1"],["Estimate"]);
tasall = renamevars([taSsol.Coefficients(:,1); tasstat2],["Estimate"],["TAS"]);
taall = [tafall tasall];

figure()
scatter(obsF(1:15,2),obsF(1:15,1));
hold on
scatter(obsF(1:15,2),taFsol.Fitted(1:15));
hold on
scatter(obsF(16:30,2),obsF(16:30,1));
hold on
scatter(obsF(16:30,2),taFsol.Fitted(16:30));
title('TAF Analytical')
xlabel('Time (hrs)')
ylabel('DOCw (mg)')
legend('low obs','low model','high obs','high model')
saveas(gcf,['TAFan'],'png');
figure()
scatter(obsS(1:15,2),obsS(1:15,1));
hold on
scatter(obsS(1:15,2),taSsol.Fitted(1:15));
hold on
scatter(obsS(16:30,2),obsS(16:30,1));
hold on
scatter(obsS(16:30,2),taSsol.Fitted(16:30));
title('TAS Analytical')
xlabel('Time (hrs)')
ylabel('DOCw (mg)')
legend('low obs','low model','high obs','high model')
saveas(gcf,['TASan'],'png');
%% WC 
% Observational values over time (ordered L's in triplicate then H's)
obs = xlsread('obs_analytical_WC.xlsx');
obsF = obs(1:30,:);
obsS = obs(31:57,:);

%y0 = [LS0 HS0 LW0 HW0]
y0F = [[y0all(23,2) y0all(21,2)];[y0all(23,1) y0all(21,1)]];
beta = [kads(5) kdes(5)];
wcFsol = fitnlm(obsF(:,2), obsF(:,1), @fun_analytical_FfixS0, beta);
wcFest = renamevars(wcFsol.Coefficients(:,1),["Estimate"],["WCF rates"]);

%model output
%fresh
    ka = beta(1);
    kd = beta(2);
    LS0 = y0F(1,1);  %low S0
    HS0 = y0F(1,2);  %high S0
    LW0 = y0F(2,1);  %low w0
    HW0 = y0F(2,2);  %high w0
    x = tmod2;
    z = exp(-(ka + kd)*x);
    y = zeros(482,1);
    for i = 1:241
        y(i) = (1 / (ka + kd))*((LS0 * kd * (1 - z(i))) + (LW0 * (ka + (kd * z(i)))));
    end
    for i = 242:482
        y(i) = (1 / (ka + kd))*((HS0 * kd * (1 - z(i))) + (HW0 * (ka + (kd * z(i)))));
    end
wcfmod = y;


y0S = [[y0all(24,2) y0all(22,2)];[y0all(24,1) y0all(22,1)]];
beta = [kads(6) kdes(6)];
wcSsol = fitnlm(obsS(:,2), obsS(:,1), @fun_analytical_SfixS0, beta);
wcSest = renamevars(wcSsol.Coefficients(:,1),["Estimate"],["WCS rates"]);

%model output
%saline
    ka = beta(1);
    kd = beta(2);
    LS0 = y0S(1,1);  %low S0
    HS0 = y0S(1,2);  %high S0
    LW0 = y0S(2,1);  %low w0
    HW0 = y0S(2,2);  %high w0
    x = tmod2;
    z = exp(-(ka + kd)*x);
    y = zeros(482,1);
    for i = 1:241
        y(i) = (1 / (ka + kd))*((LS0 * kd * (1 - z(i))) + (LW0 * (ka + (kd * z(i)))));
    end
    for i = 242:482
        y(i) = (1 / (ka + kd))*((HS0 * kd * (1 - z(i))) + (HW0 * (ka + (kd * z(i)))));
    end
wcsmod = y;
wcmod = [wcfmod wcsmod];

%estimate output
WCest = [wcFest wcSest];
wcfstat = [wcFsol.Rsquared.Adjusted(1,1) wcFsol.ModelCriterion.AIC(1,1) wcFsol.ModelCriterion.AICc(1,1) wcFsol.RMSE(1,1)]; 
wcfstat2 = renamevars(array2table(wcfstat'),["Var1"],["Estimate"]);
wcfall = renamevars([wcFsol.Coefficients(:,1); wcfstat2],["Estimate"],["WCF"]);

wcsstat = [wcSsol.Rsquared.Adjusted(1,1) wcSsol.ModelCriterion.AIC(1,1) wcSsol.ModelCriterion.AICc(1,1) wcSsol.RMSE(1,1)]; 
wcsstat2 = renamevars(array2table(wcsstat'),["Var1"],["Estimate"]);
wcsall = renamevars([wcSsol.Coefficients(:,1); wcsstat2],["Estimate"],["WCS"]);
wcall = [wcfall wcsall];

figure()
scatter(obsF(1:15,2),obsF(1:15,1));
hold on
scatter(obsF(1:15,2),wcFsol.Fitted(1:15));
hold on
scatter(obsF(16:30,2),obsF(16:30,1));
hold on
scatter(obsF(16:30,2),wcFsol.Fitted(16:30));
title('WCF Analytical')
xlabel('Time (hrs)')
ylabel('DOCw (mg)')
legend('low obs','low model','high obs','high model')
saveas(gcf,['WCFan'],'png');
figure()
scatter(obsS(1:12,2),obsS(1:12,1));
hold on
scatter(obsS(1:12,2),wcSsol.Fitted(1:12));
hold on
scatter(obsS(13:27,2),obsS(13:27,1));
hold on
scatter(obsS(13:27,2),wcSsol.Fitted(13:27));
title('WCS Analytical')
xlabel('Time (hrs)')
ylabel('DOCw (mg)')
legend('low obs','low model','high obs','high model')
saveas(gcf,['WCSan'],'png');

%% WI
% Observational values over time (ordered L's in triplicate then H's)
obs = xlsread('obs_analytical_WI.xlsx');
obsF = obs(1:33,:);
obsS = obs(34:63,:);

%y0 = [LS0 HS0 LW0 HW0]
y0F = [[y0all(27,2) y0all(25,2)];[y0all(27,1) y0all(25,1)]];
beta = [kads(7) kdes(7)];
wiFsol = fitnlm(obsF(:,2), obsF(:,1), @fun_analytical_FfixS0, beta);
wiFest = renamevars(wiFsol.Coefficients(:,1),["Estimate"],["WIF rates"]);

%model output
%fresh
    ka = beta(1);
    kd = beta(2);
    LS0 = y0F(1,1);  %low S0
    HS0 = y0F(1,2);  %high S0
    LW0 = y0F(2,1);  %low w0
    HW0 = y0F(2,2);  %high w0
    x = tmod2;
    z = exp(-(ka + kd)*x);
    y = zeros(482,1);
    for i = 1:241
        y(i) = (1 / (ka + kd))*((LS0 * kd * (1 - z(i))) + (LW0 * (ka + (kd * z(i)))));
    end
    for i = 242:482
        y(i) = (1 / (ka + kd))*((HS0 * kd * (1 - z(i))) + (HW0 * (ka + (kd * z(i)))));
    end
wifmod = y;

y0S = [[y0all(28,2) y0all(26,2)];[y0all(28,1) y0all(26,1)]];
beta = [kads(8) kdes(8)];
wiSsol = fitnlm(obsS(:,2), obsS(:,1), @fun_analytical_SfixS0, beta);
wiSest = renamevars(wiSsol.Coefficients(:,1),["Estimate"],["WIS rates"]);

%model output
%saline
    ka = beta(1);
    kd = beta(2);
    LS0 = y0S(1,1);  %low S0
    HS0 = y0S(1,2);  %high S0
    LW0 = y0S(2,1);  %low w0
    HW0 = y0S(2,2);  %high w0
    x = tmod2;
    z = exp(-(ka + kd)*x);
    y = zeros(482,1);
    for i = 1:241
        y(i) = (1 / (ka + kd))*((LS0 * kd * (1 - z(i))) + (LW0 * (ka + (kd * z(i)))));
    end
    for i = 242:482
        y(i) = (1 / (ka + kd))*((HS0 * kd * (1 - z(i))) + (HW0 * (ka + (kd * z(i)))));
    end
wismod = y;
wimod = [wifmod wismod];

%estimate output
WIest = [wiFest wiSest];
wifstat = [wiFsol.Rsquared.Adjusted(1,1) wiFsol.ModelCriterion.AIC(1,1) wiFsol.ModelCriterion.AICc(1,1) wiFsol.RMSE(1,1)]; 
wifstat2 = renamevars(array2table(wifstat'),["Var1"],["Estimate"]);
wifall = renamevars([wiFsol.Coefficients(:,1); wifstat2],["Estimate"],["WIF"]);

wisstat = [wiSsol.Rsquared.Adjusted(1,1) wiSsol.ModelCriterion.AIC(1,1) wiSsol.ModelCriterion.AICc(1,1) wiSsol.RMSE(1,1)]; 
wisstat2 = renamevars(array2table(wisstat'),["Var1"],["Estimate"]);
wisall = renamevars([wiSsol.Coefficients(:,1); wisstat2],["Estimate"],["WIS"]);
wiall = [wifall wisall];

figure()
scatter(obsF(1:15,2),obsF(1:15,1));
hold on
scatter(obsF(1:15,2),wiFsol.Fitted(1:15));
hold on
scatter(obsF(16:33,2),obsF(16:33,1));
hold on
scatter(obsF(16:33,2),wiFsol.Fitted(16:33));
title('WIF Analytical')
xlabel('Time (hrs)')
ylabel('DOCw (mg)')
legend('low obs','low model','high obs','high model')
saveas(gcf,['WIFan'],'png');
figure()
scatter(obsS(1:15,2),obsS(1:15,1));
hold on
scatter(obsS(1:15,2),wiSsol.Fitted(1:15));
hold on
scatter(obsS(16:30,2),obsS(16:30,1));
hold on
scatter(obsS(16:30,2),wiSsol.Fitted(16:30));
title('WIS Analytical')
xlabel('Time (hrs)')
ylabel('DOCw (mg)')
legend('low obs','low model','high obs','high model')
saveas(gcf,['WISan'],'png');
%% WM
% Observational values over time (ordered L's in triplicate then H's)
obs = xlsread('obs_analytical_WM.xlsx');
obsF = obs(1:30,:);
obsS = obs(31:60,:);

%y0 = [LS0 HS0 LW0 HW0]
y0F = [[y0all(31,2) y0all(29,2)];[y0all(31,1) y0all(29,1)]];
beta = [kads(9) kdes(9)];
wmFsol = fitnlm(obsF(:,2), obsF(:,1), @fun_analytical_FfixS0, beta);
wmFest = renamevars(wmFsol.Coefficients(:,1),["Estimate"],["WMF rates"]);

%model output
%fresh
    ka = beta(1);
    kd = beta(2);
    LS0 = y0F(1,1);  %low S0
    HS0 = y0F(1,2);  %high S0
    LW0 = y0F(2,1);  %low w0
    HW0 = y0F(2,2);  %high w0
    x = tmod2;
    z = exp(-(ka + kd)*x);
    y = zeros(482,1);
    for i = 1:241
        y(i) = (1 / (ka + kd))*((LS0 * kd * (1 - z(i))) + (LW0 * (ka + (kd * z(i)))));
    end
    for i = 242:482
        y(i) = (1 / (ka + kd))*((HS0 * kd * (1 - z(i))) + (HW0 * (ka + (kd * z(i)))));
    end
wmfmod = y;

y0S = [[y0all(32,2) y0all(30,2)];[y0all(32,1) y0all(30,1)]];
beta = [kads(10) kdes(10)];
wmSsol = fitnlm(obsS(:,2), obsS(:,1), @fun_analytical_SfixS0, beta);
wmSest = renamevars(wmSsol.Coefficients(:,1),["Estimate"],["WMS rates"]);

%model output
%saline
    ka = beta(1);
    kd = beta(2);
    LS0 = y0S(1,1);  %low S0
    HS0 = y0S(1,2);  %high S0
    LW0 = y0S(2,1);  %low w0
    HW0 = y0S(2,2);  %high w0
    x = tmod2;
    z = exp(-(ka + kd)*x);
    y = zeros(482,1);
    for i = 1:241
        y(i) = (1 / (ka + kd))*((LS0 * kd * (1 - z(i))) + (LW0 * (ka + (kd * z(i)))));
    end
    for i = 242:482
        y(i) = (1 / (ka + kd))*((HS0 * kd * (1 - z(i))) + (HW0 * (ka + (kd * z(i)))));
    end
wmsmod = y;
wmmod = [wmfmod wmsmod];

%estimate output
WMest = [wmFest wmSest];
wmfstat = [wmFsol.Rsquared.Adjusted(1,1) wmFsol.ModelCriterion.AIC(1,1) wmFsol.ModelCriterion.AICc(1,1) wmFsol.RMSE(1,1)]; 
wmfstat2 = renamevars(array2table(wmfstat'),["Var1"],["Estimate"]);
wmfall = renamevars([wmFsol.Coefficients(:,1); wmfstat2],["Estimate"],["WMF"]);

wmsstat = [wmSsol.Rsquared.Adjusted(1,1) wmSsol.ModelCriterion.AIC(1,1) wmSsol.ModelCriterion.AICc(1,1) wmSsol.RMSE(1,1)]; 
wmsstat2 = renamevars(array2table(wmsstat'),["Var1"],["Estimate"]);
wmsall = renamevars([wmSsol.Coefficients(:,1); wmsstat2],["Estimate"],["WMS"]);
wmall = [wmfall wmsall];

figure()
scatter(obsF(1:15,2),obsF(1:15,1));
hold on
scatter(obsF(1:15,2),wmFsol.Fitted(1:15));
hold on
scatter(obsF(16:30,2),obsF(16:30,1));
hold on
scatter(obsF(16:30,2),wmFsol.Fitted(16:30));
title('WMF Analytical')
xlabel('Time (hrs)')
ylabel('DOCw (mg)')
legend('low obs','low model','high obs','high model')
saveas(gcf,['WMFan'],'png');
figure()
scatter(obsS(1:15,2),obsS(1:15,1));
hold on
scatter(obsS(1:15,2),wmSsol.Fitted(1:15));
hold on
scatter(obsS(16:30,2),obsS(16:30,1));
hold on
scatter(obsS(16:30,2),wmSsol.Fitted(16:30));
title('WMS Analytical')
xlabel('Time (hrs)')
ylabel('DOCw (mg)')
legend('low obs','low model','high obs','high model')
saveas(gcf,['WMSan'],'png');
%% PC
% Observational values over time (ordered L's in triplicate then H's)
obs = xlsread('obs_analytical_PC.xlsx');
obsF = obs(1:30,:);
obsS = obs(31:60,:);

%y0 = [LS0 HS0 LW0 HW0]
y0F = [[y0all(7,2) y0all(5,2)];[y0all(7,1) y0all(5,1)]];
beta = [kads(11) kdes(11)];
pcFsol = fitnlm(obsF(:,2), obsF(:,1), @fun_analytical_FfixS0, beta);
pcFest = renamevars(pcFsol.Coefficients(:,1),["Estimate"],["PCF rates"]);

%model output
%fresh
    ka = beta(1);
    kd = beta(2);
    LS0 = y0F(1,1);  %low S0
    HS0 = y0F(1,2);  %high S0
    LW0 = y0F(2,1);  %low w0
    HW0 = y0F(2,2);  %high w0
    x = tmod2;
    z = exp(-(ka + kd)*x);
    y = zeros(482,1);
    for i = 1:241
        y(i) = (1 / (ka + kd))*((LS0 * kd * (1 - z(i))) + (LW0 * (ka + (kd * z(i)))));
    end
    for i = 242:482
        y(i) = (1 / (ka + kd))*((HS0 * kd * (1 - z(i))) + (HW0 * (ka + (kd * z(i)))));
    end
pcfmod = y;

y0S = [[y0all(8,2) y0all(6,2)];[y0all(8,1) y0all(6,1)]];
beta = [kads(12) kdes(12)];
pcSsol = fitnlm(obsS(:,2), obsS(:,1), @fun_analytical_SfixS0, beta);
pcSest = renamevars(pcSsol.Coefficients(:,1),["Estimate"],["PCS rates"]);

%model output
%saline
    ka = beta(1);
    kd = beta(2);
    LS0 = y0S(1,1);  %low S0
    HS0 = y0S(1,2);  %high S0
    LW0 = y0S(2,1);  %low w0
    HW0 = y0S(2,2);  %high w0
    x = tmod2;
    z = exp(-(ka + kd)*x);
    y = zeros(482,1);
    for i = 1:241
        y(i) = (1 / (ka + kd))*((LS0 * kd * (1 - z(i))) + (LW0 * (ka + (kd * z(i)))));
    end
    for i = 242:482
        y(i) = (1 / (ka + kd))*((HS0 * kd * (1 - z(i))) + (HW0 * (ka + (kd * z(i)))));
    end
pcsmod = y;
pcmod = [pcfmod pcsmod];

%estimate output
PCest = [pcFest pcSest];
pcfstat = [pcFsol.Rsquared.Adjusted(1,1) pcFsol.ModelCriterion.AIC(1,1) pcFsol.ModelCriterion.AICc(1,1) pcFsol.RMSE(1,1)]; 
pcfstat2 = renamevars(array2table(pcfstat'),["Var1"],["Estimate"]);
pcfall = renamevars([pcFsol.Coefficients(:,1); pcfstat2],["Estimate"],["PCF"]);

pcsstat = [pcSsol.Rsquared.Adjusted(1,1) pcSsol.ModelCriterion.AIC(1,1) pcSsol.ModelCriterion.AICc(1,1) pcSsol.RMSE(1,1)]; 
pcsstat2 = renamevars(array2table(pcsstat'),["Var1"],["Estimate"]);
pcsall = renamevars([pcSsol.Coefficients(:,1); pcsstat2],["Estimate"],["PCS"]);
pcall = [pcfall pcsall];

figure()
scatter(obsF(1:15,2),obsF(1:15,1));
hold on
scatter(obsF(1:15,2),pcFsol.Fitted(1:15));
hold on
scatter(obsF(16:30,2),obsF(16:30,1));
hold on
scatter(obsF(16:30,2),pcFsol.Fitted(16:30));
title('PCF Analytical')
xlabel('Time (hrs)')
ylabel('DOCw (mg)')
legend('low obs','low model','high obs','high model')
saveas(gcf,['PCFan'],'png');
figure()
scatter(obsS(1:15,2),obsS(1:15,1));
hold on
scatter(obsS(1:15,2),pcSsol.Fitted(1:15));
hold on
scatter(obsS(16:30,2),obsS(16:30,1));
hold on
scatter(obsS(16:30,2),pcSsol.Fitted(16:30));
title('PCS Analytical')
xlabel('Time (hrs)')
ylabel('DOCw (mg)')
legend('low obs','low model','high obs','high model')
saveas(gcf,['PCSan'],'png');
%% PI
% Observational values over time (ordered L's in triplicate then H's)
obs = xlsread('obs_analytical_PI.xlsx');
obsF = obs(1:33,:);
obsS = obs(34:63,:);

%y0 = [LS0 HS0 LW0 HW0]
y0F = [[y0all(11,2) y0all(9,2)];[y0all(11,1) y0all(9,1)]];
beta = [kads(13) kdes(13)];
piFsol = fitnlm(obsF(:,2), obsF(:,1), @fun_analytical_FfixS0, beta);
piFest = renamevars(piFsol.Coefficients(:,1),["Estimate"],["PIF rates"]);

%model output
%fresh
    ka = beta(1);
    kd = beta(2);
    LS0 = y0F(1,1);  %low S0
    HS0 = y0F(1,2);  %high S0
    LW0 = y0F(2,1);  %low w0
    HW0 = y0F(2,2);  %high w0
    x = tmod2;
    z = exp(-(ka + kd)*x);
    y = zeros(482,1);
    for i = 1:241
        y(i) = (1 / (ka + kd))*((LS0 * kd * (1 - z(i))) + (LW0 * (ka + (kd * z(i)))));
    end
    for i = 242:482
        y(i) = (1 / (ka + kd))*((HS0 * kd * (1 - z(i))) + (HW0 * (ka + (kd * z(i)))));
    end
pifmod = y;

y0S = [[y0all(12,2) y0all(10,2)];[y0all(12,1) y0all(10,1)]];
beta = [kads(14) kdes(14)];
piSsol = fitnlm(obsS(:,2), obsS(:,1), @fun_analytical_SfixS0, beta);
piSest = renamevars(piSsol.Coefficients(:,1),["Estimate"],["PIS rates"]);

%model output
%saline
    ka = beta(1);
    kd = beta(2);
    LS0 = y0S(1,1);  %low S0
    HS0 = y0S(1,2);  %high S0
    LW0 = y0S(2,1);  %low w0
    HW0 = y0S(2,2);  %high w0
    x = tmod2;
    z = exp(-(ka + kd)*x);
    y = zeros(482,1);
    for i = 1:241
        y(i) = (1 / (ka + kd))*((LS0 * kd * (1 - z(i))) + (LW0 * (ka + (kd * z(i)))));
    end
    for i = 242:482
        y(i) = (1 / (ka + kd))*((HS0 * kd * (1 - z(i))) + (HW0 * (ka + (kd * z(i)))));
    end
pismod = y;
pimod = [pifmod pismod];

%estimate output
PIest = [piFest piSest];
pifstat = [piFsol.Rsquared.Adjusted(1,1) piFsol.ModelCriterion.AIC(1,1) piFsol.ModelCriterion.AICc(1,1) piFsol.RMSE(1,1)]; 
pifstat2 = renamevars(array2table(pifstat'),["Var1"],["Estimate"]);
pifall = renamevars([piFsol.Coefficients(:,1); pifstat2],["Estimate"],["PIF"]);

pisstat = [piSsol.Rsquared.Adjusted(1,1) piSsol.ModelCriterion.AIC(1,1) piSsol.ModelCriterion.AICc(1,1) piSsol.RMSE(1,1)]; 
pisstat2 = renamevars(array2table(pisstat'),["Var1"],["Estimate"]);
pisall = renamevars([piSsol.Coefficients(:,1); pisstat2],["Estimate"],["PIS"]);
piall = [pifall pisall];

figure()
scatter(obsF(1:15,2),obsF(1:15,1));
hold on
scatter(obsF(1:15,2),piFsol.Fitted(1:15));
hold on
scatter(obsF(16:33,2),obsF(16:33,1));
hold on
scatter(obsF(16:33,2),piFsol.Fitted(16:33));
title('PIF Analytical')
xlabel('Time (hrs)')
ylabel('DOCw (mg)')
legend('low obs','low model','high obs','high model')
saveas(gcf,['PIFan'],'png');
figure()
scatter(obsS(1:15,2),obsS(1:15,1));
hold on
scatter(obsS(1:15,2),piSsol.Fitted(1:15));
hold on
scatter(obsS(16:30,2),obsS(16:30,1));
hold on
scatter(obsS(16:30,2),piSsol.Fitted(16:30));
title('PIS Analytical')
xlabel('Time (hrs)')
ylabel('DOCw (mg)')
legend('low obs','low model','high obs','high model')
saveas(gcf,['PISan'],'png');
%% PM
% Observational values over time (ordered L's in triplicate then H's)
obs = xlsread('obs_analytical_PM.xlsx');
obsF = obs(1:36,:);
obsS = obs(37:66,:);

%y0 = [LS0 HS0 LW0 HW0]
y0F = [[y0all(15,2) y0all(13,2)];[y0all(15,1) y0all(13,1)]];
beta = [kads(15) kdes(15)];
pmFsol = fitnlm(obsF(:,2), obsF(:,1), @fun_analytical_FfixS0, beta);
pmFest = renamevars(pmFsol.Coefficients(:,1),["Estimate"],["PMF rates"]);

%model output
%fresh
    ka = beta(1);
    kd = beta(2);
    LS0 = y0F(1,1);  %low S0
    HS0 = y0F(1,2);  %high S0
    LW0 = y0F(2,1);  %low w0
    HW0 = y0F(2,2);  %high w0
    x = tmod2;
    z = exp(-(ka + kd)*x);
    y = zeros(482,1);
    for i = 1:241
        y(i) = (1 / (ka + kd))*((LS0 * kd * (1 - z(i))) + (LW0 * (ka + (kd * z(i)))));
    end
    for i = 242:482
        y(i) = (1 / (ka + kd))*((HS0 * kd * (1 - z(i))) + (HW0 * (ka + (kd * z(i)))));
    end
pmfmod = y;

y0S = [[y0all(16,2) y0all(14,2)];[y0all(16,1) y0all(14,1)]];
beta = [kads(16) kdes(16)];
pmSsol = fitnlm(obsS(:,2), obsS(:,1), @fun_analytical_SfixS0, beta);
pmSest = renamevars(pmSsol.Coefficients(:,1),["Estimate"],["PMS rates"]);

%model output
%saline
    ka = beta(1);
    kd = beta(2);
    LS0 = y0S(1,1);  %low S0
    HS0 = y0S(1,2);  %high S0
    LW0 = y0S(2,1);  %low w0
    HW0 = y0S(2,2);  %high w0
    x = tmod2;
    z = exp(-(ka + kd)*x);
    y = zeros(482,1);
    for i = 1:241
        y(i) = (1 / (ka + kd))*((LS0 * kd * (1 - z(i))) + (LW0 * (ka + (kd * z(i)))));
    end
    for i = 242:482
        y(i) = (1 / (ka + kd))*((HS0 * kd * (1 - z(i))) + (HW0 * (ka + (kd * z(i)))));
    end
pmsmod = y;
pmmod = [pmfmod pmsmod];

%estimate output
PMest = [pmFest pmSest];
pmfstat = [pmFsol.Rsquared.Adjusted(1,1) pmFsol.ModelCriterion.AIC(1,1) pmFsol.ModelCriterion.AICc(1,1) pmFsol.RMSE(1,1)]; 
pmfstat2 = renamevars(array2table(pmfstat'),["Var1"],["Estimate"]);
pmfall = renamevars([pmFsol.Coefficients(:,1); pmfstat2],["Estimate"],["PMF"]);

pmsstat = [pmSsol.Rsquared.Adjusted(1,1) pmSsol.ModelCriterion.AIC(1,1) pmSsol.ModelCriterion.AICc(1,1) pmSsol.RMSE(1,1)]; 
pmsstat2 = renamevars(array2table(pmsstat'),["Var1"],["Estimate"]);
pmsall = renamevars([pmSsol.Coefficients(:,1); pmsstat2],["Estimate"],["PMS"]);
pmall = [pmfall pmsall];

figure()
scatter(obsF(1:15,2),obsF(1:15,1));
hold on
scatter(obsF(1:15,2),pmFsol.Fitted(1:15));
hold on
scatter(obsF(16:36,2),obsF(16:36,1));
hold on
scatter(obsF(16:36,2),pmFsol.Fitted(16:36));
title('PMF Analytical')
xlabel('Time (hrs)')
ylabel('DOCw (mg)')
legend('low obs','low model','high obs','high model')
saveas(gcf,['PMFan'],'png');
figure()
scatter(obsS(1:15,2),obsS(1:15,1));
hold on
scatter(obsS(1:15,2),pmSsol.Fitted(1:15));
hold on
scatter(obsS(16:30,2),obsS(16:30,1));
hold on
scatter(obsS(16:30,2),pmSsol.Fitted(16:30));
title('PMS Analytical')
xlabel('Time (hrs)')
ylabel('DOCw (mg)')
legend('low obs','low model','high obs','high model')
saveas(gcf,['PMSan'],'png');
%% Data Output

AnRates = [JBest TAest WCest WIest WMest PCest PIest PMest];
writetable(AnRates,'newAnalyticalRates.csv')

anout = [jball taall wcall wiall wmall pcall piall pmall];
writetable(anout,'newAnalyticalOutput.xlsx')

anmodall = array2table([jbmod tamod wcmod wimod wmmod pcmod pimod pmmod]);
anmodall.Properties.VariableNames = {'JBF' 'JBS' 'TAF' 'TAS' 'WCF' 'WCS' 'WIF' 'WIS' 'WMF' 'WMS' 'PCF' 'PCS' 'PIF' 'PIS' 'PMF' 'PMS'};
writetable(anmodall,'newAnalyticalModel.xlsx');
