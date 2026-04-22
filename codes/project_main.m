%% PROJECT - ISCONVERSIONAL METHODS
clear all
close all
clc

% ----- PLOT SETUP
set(0, 'defaultTextInterpreter', 'latex')
set(groot, 'defaultAxesTickLabelInterpreter', 'latex')
set(0, 'defaultLegendInterpreter', 'latex')
set(0, 'defaultFigureColormap', turbo(256));
set(0, 'defaultSurfaceEdgeAlpha', 0.3);
set(0, 'defaultLineLineWidth', 1);
set(0, 'defaultFigureColor', [1, 1, 1]);
set(0, 'defaultAxesColor', 'none');
set(groot, 'DefaultAxesFontSize', 17);

% ----- ADDPATHS 
% RawData folder with experimental data 
addpath("RawData\AlWO3_ACTPOWDERS\");
addpath("RawData\AlWO3_LOOSEPOWDER\");
addpath("RawData\MgSiO2_DLP\");
addpath("RawData\MgSiO2_LOOSEPOWDER\");

% ----- LOAD MATRICES
load("OFW_table.mat");                                       % OFW tables: iterative procedure with integration

% AlWO3_LOOSEPOWDER
data.AlWO3_LOOSE_2p5 = readmatrix("ExpDat_20260220_Al_WO3_60_40_loose_Ar_2p5.txt",'CommentStyle', '#');
data.AlWO3_LOOSE_5   = readmatrix("ExpDat_20260220_Al_WO3_60_40_loose_Ar_5p0.txt",'CommentStyle', '#');
data.AlWO3_LOOSE_10  = readmatrix("ExpDat_20260220_Al_WO3_60_40_loose_Ar_10p0.txt",'CommentStyle', '#');
data.AlWO3_LOOSE_20  = readmatrix("ExpDat_20260220_Al_WO3_60_40_loose_Ar_20p0.txt",'CommentStyle', '#');

% AlWO3_ACTPOWDERS @ 20 K/min
data1.AlWO3_ACT_40_60    = readmatrix("ExpDat_20260225_Al_WO3_40_60_act_Ar_20p0.txt",'CommentStyle', '#');
data1.AlWO3_ACT_30_70_R1 = readmatrix("ExpDat_20260225_Al_WO3_30_70_act_Ar_20p0_R1.txt",'CommentStyle', '#');
data1.AlWO3_ACT_30_70_R2 = readmatrix("ExpDat_20260225_Al_WO3_30_70_act_Ar_20p0_R2.txt",'CommentStyle', '#');

% MgSiO2_LOOSEPOWDER
data2.MgSiO2_LOOSE_2p5  = readmatrix("ExpDat_20260220_Mg_SiO2_stoich_loose_Ar_2p5.txt",'CommentStyle', '#');
data2.MgSiO2_LOOSE_5_R1 = readmatrix("ExpDat_20260224_Mg_SiO2_stoich_loose_Ar_5p0_R1.txt",'CommentStyle', '#');
data2.MgSiO2_LOOSE_5_R2 = readmatrix("ExpDat_20260224_Mg_SiO2_stoich_loose_Ar_5p0_R2.txt",'CommentStyle', '#');
data2.MgSiO2_LOOSE_10   = readmatrix("ExpDat_20260224_Mg_SiO2_stoich_loose_Ar_10p0.txt",'CommentStyle', '#');
data2.MgSiO2_LOOSE_20   = readmatrix("ExpDat_20260220_Mg_SiO2_stoich_loose_Ar_20p0.txt",'CommentStyle', '#');

% MgSiO2_DLP
data3.MgSiO2_ACT_2p5 = readmatrix("ExpDat_20260218_Mg_SiO2_stoich_60wt_high_speed_Ar_2p5.txt",'CommentStyle', '#');
data3.MgSiO2_ACT_5   = readmatrix("ExpDat_20260218_Mg_SiO2_stoich_60wt_high_speed_Ar_5p0.txt",'CommentStyle', '#');
data3.MgSiO2_ACT_10  = readmatrix("ExpDat_20260211_Mg_SiO2_stoich_60wt_high_speed_Ar_10p0.txt",'CommentStyle', '#');
data3.MgSiO2_ACT_20  = readmatrix("ExpDat_20260211_Mg_SiO2_stoich_60wt_high_speed_Ar_20p0.txt",'CommentStyle', '#');

% ----- OTHER DATA

% Heating rates
gammaA = 2.5;                                                              % [°C/min], [K/min]
gammaB = 5;                                                                % [°C/min], [K/min]
gammaC = 10;                                                               % [°C/min], [K/min]
gammaD = 20;                                                               % [°C/min], [K/min]
gammaVect = [gammaA gammaB gammaC gammaD];

% Lifetime Isothermal COndition
T_LFT = 650 + 273.15;                                                      % [K]

% Conversion Degree
% alpha = (m_max - m(t))/(m_max - m_min)
NAlpha = 100; % arbitrary
alpha = linspace(0.1, 0.9, NAlpha);                                        % [-] -- usually this span, attention to ortega method 

% ----- SORT DATA

% DATA
names = fieldnames(data);

for i = 1:numel(names)
    name = names{i};
    datas = data.(name);

    inp.(name).T        = datas(:,1);
    inp.(name).t        = datas(:,2);
    inp.(name).DTAcurve = datas(:,3);
    inp.(name).massPerc = datas(:,4);
    inp.(name).TGA      = datas(:,4) / datas(1,4);

end

% ----- PLOTS DTA
names = fieldnames(inp);

figure
hold on
grid on

for i = 1:numel(names)
    name = names{i};
    plot(inp.(name).T, inp.(name).DTAcurve)
end

legend(names, 'Interpreter', 'none')
xlabel('T')
ylabel('DTA')
title('DTA vs Temperature')

% ----- PLOTS MASS%
figure
hold on
grid on

for i = 1:numel(names)
    name = names{i};
    plot(inp.(name).T, inp.(name).massPerc)
end

legend(names, 'Interpreter', 'none')
xlabel('T')
ylabel('Mass \%')
title('Mass \% vs Temperature')

% DATA 1
names = fieldnames(data1);

for i = 1:numel(names)
    name = names{i};
    datas = data1.(name);

    inp1.(name).T        = datas(:,1);
    inp1.(name).t        = datas(:,2);
    inp1.(name).DTAcurve = datas(:,3);
    inp1.(name).massPerc = datas(:,4);
    inp1.(name).TGA      = datas(:,4) / datas(1,4);

end

% ----- PLOTS DTA
names = fieldnames(inp1);

figure
hold on
grid on

for i = 1:numel(names)
    name = names{i};
    plot(inp1.(name).T, inp1.(name).DTAcurve)
end

legend(names, 'Interpreter', 'none')
xlabel('T')
ylabel('DTA')
title('DTA vs Temperature')

% ----- PLOTS MASS%
figure
hold on
grid on

for i = 1:numel(names)
    name = names{i};
    plot(inp1.(name).T, inp1.(name).massPerc)
end

legend(names, 'Interpreter', 'none')
xlabel('T')
ylabel('Mass \%')
title('Mass \% vs Temperature')

% DATA 2
names = fieldnames(data2);

for i = 1:numel(names)
    name = names{i};
    datas = data2.(name);

    inp2.(name).T        = datas(:,1);
    inp2.(name).t        = datas(:,2);
    inp2.(name).DTAcurve = datas(:,3);
    inp2.(name).massPerc = datas(:,4);
    inp2.(name).TGA      = datas(:,4) / datas(1,4);

end

% ----- PLOTS DTA
names = fieldnames(inp2);

figure
hold on
grid on

for i = 1:numel(names)
    name = names{i};
    plot(inp2.(name).T, inp2.(name).DTAcurve)
end

legend(names, 'Interpreter', 'none')
xlabel('T')
ylabel('DTA')
title('DTA vs Temperature')

% ----- PLOTS MASS%
figure
hold on
grid on

for i = 1:numel(names)
    name = names{i};
    plot(inp2.(name).T, inp2.(name).massPerc)
end

legend(names, 'Interpreter', 'none')
xlabel('T')
ylabel('Mass \%')
title('Mass \% vs Temperature')

% DATA 3
names = fieldnames(data3);

for i = 1:numel(names)
    name = names{i};
    datas = data3.(name);

    inp3.(name).T        = datas(:,1);
    inp3.(name).t        = datas(:,2);
    inp3.(name).DTAcurve = datas(:,3);
    inp3.(name).massPerc = datas(:,4);
    inp3.(name).TGA      = datas(:,4) / datas(1,4);

end

% ----- PLOTS DTA
names = fieldnames(inp3);

figure
hold on
grid on

for i = 1:numel(names)
    name = names{i};
    plot(inp3.(name).T, inp3.(name).DTAcurve)
end

legend(names, 'Interpreter', 'none')
xlabel('T')
ylabel('DTA')
title('DTA vs Temperature')

% ----- PLOTS MASS%
figure
hold on
grid on

for i = 1:numel(names)
    name = names{i};
    plot(inp3.(name).T, inp3.(name).massPerc)
end

legend(names, 'Interpreter', 'none')
xlabel('T')
ylabel('Mass \%')
title('Mass \% vs Temperature')

% % ----- PLOTS MASS %
% figure
% hold on
% grid on
% 
% for i = 1:numel(names)
%     name = names{i};
%     plot(inp.(name).T, inp.(name).TGA)
% end
% 
% legend(names, 'Interpreter', 'none')
% xlabel('T')
% ylabel('TGA \%')
% title('TGA vs Temperature')




%%
% Strcuct miscInputs
miscInput.OFW_intTables = OFW_intTables.ofwtable;
miscInput.NinterpPts    = 10^5;
miscInput.refTemp       = 110;
miscInput.T_LFT         = T_LFT;
miscInput.deltaAlpha    = 0.2;
method = 'OFW';
method1 = 'KAS';
method2 = 'ORT';

%% RUN OFW MODEL

% actEnergyEval: takes the integral model you want, the inputs of the solid
% and the environment fluid, the inputs  

% OFW MODEL: AIR
[out.OFW.S2p5W1CB.AIR] = actEnergyEval(method, inp.S2p5W1CB.AIR, miscInput, alpha);
[out.OFW.SEBSma.AIR]   = actEnergyEval(method, inp.SEBSma.AIR, miscInput, alpha);
[out.OFW.SASOL.AIR]    = actEnergyEval(method, inp.SASOL.AIR, miscInput, alpha);

% OFW MODEL : Ar
[out.OFW.S2p5W1CB.Ar]  = actEnergyEval(method, inp.S2p5W1CB.Ar, miscInput, alpha);
[out.OFW.SEBSma.Ar]    = actEnergyEval(method, inp.SEBSma.Ar, miscInput, alpha);
[out.OFW.SASOL.Ar]     = actEnergyEval(method, inp.SASOL.Ar, miscInput, alpha);

%% RUN KAS MODEL

% KAS MODEL: AIR
[out.KAS.S2p5W1CB.AIR] = actEnergyEval(method1, inp.S2p5W1CB.AIR, miscInput, alpha);
[out.KAS.SEBSma.AIR]   = actEnergyEval(method1, inp.SEBSma.AIR, miscInput, alpha);
[out.KAS.SASOL.AIR]    = actEnergyEval(method1, inp.SASOL.AIR, miscInput, alpha);

% KAS MODEL : Ar
[out.KAS.S2p5W1CB.Ar]  = actEnergyEval(method1, inp.S2p5W1CB.Ar, miscInput, alpha);
[out.KAS.SEBSma.Ar]    = actEnergyEval(method1, inp.SEBSma.Ar, miscInput, alpha);
[out.KAS.SASOL.Ar]     = actEnergyEval(method1, inp.SASOL.Ar, miscInput, alpha);

%% RUN ORT MODEL

% ORT MODEL: AIR
[out.ORT.S2p5W1CB.AIR] = actEnergyEval(method1, inp.S2p5W1CB.AIR, miscInput, alpha);
[out.ORT.SEBSma.AIR]   = actEnergyEval(method1, inp.SEBSma.AIR, miscInput, alpha);
[out.ORT.SASOL.AIR]    = actEnergyEval(method1, inp.SASOL.AIR, miscInput, alpha);

% ORT MODEL : Ar
[out.ORT.S2p5W1CB.Ar]  = actEnergyEval(method1, inp.S2p5W1CB.Ar, miscInput, alpha);
[out.ORT.SEBSma.Ar]    = actEnergyEval(method1, inp.SEBSma.Ar, miscInput, alpha);
[out.ORT.SASOL.Ar]     = actEnergyEval(method1, inp.SASOL.Ar, miscInput, alpha);
