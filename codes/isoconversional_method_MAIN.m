%% ISCONVERSIONAL METHODS
clear all
close all
clc

%% ADDPATHS 

%% PLOT SETTINGS
set(0, 'defaultTextInterpreter', 'latex')
set(groot, 'defaultAxesTickLabelInterpreter', 'latex')
set(0, 'defaultLegendInterpreter', 'latex')
set(0, 'defaultFigureColormap', turbo(256));
set(0, 'defaultSurfaceEdgeAlpha', 0.3);
set(0, 'defaultLineWidth', 1);
set(0, 'defaultFigureColor', [1, 1, 1]);
set(0, 'defaultAxesColor', 'none');
set(groot, 'DefaultAxesFontSize', 17);

% LOAD DTA and dsc MATRIX
% time, temperature, TG, trace, and DTA trace columns 
% A = load("..\..\") % DTA.mat
% B = load() % DSC.mat
OFW_intTables = load("..\..\secondo semeste\ctp\codes\ofw_table.mat") % OFW tables: iterative procedure with integration

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

% Gather Inputs and create matrices: S2p5W1CB
% Air
inp.S2p5W1CB.AIR.gammaVect = gammaVect;
inp.S2p5W1CB.AIR.TGAcurve  = [A.S2p5W1CB_2p5Kmin(:, 3) A.S2p5W1CB_5Kmin(:, 3) A.S2p5W1CB_10Kmin(:, 3) A.S2p5W1CB_20Kmin(:, 3)]; % TGA curve at different heating rate - so you get a matrix, each column corrispond to a gamma
inp.S2p5W1CB.AIR.TGAtemp   = [A.S2p5W1CB_2p5Kmin(:, 2) A.S2p5W1CB_5Kmin(:, 2) A.S2p5W1CB_10Kmin(:, 2) A.S2p5W1CB_20Kmin(:, 2)]; % TGA associated temperature since at differnt gamma, you have different temperature
inp.S2p5W1CB.AIR.DTAcurve  = [A.S2p5W1CB_2p5Kmin(:, 6) A.S2p5W1CB_5Kmin(:, 6) A.S2p5W1CB_10Kmin(:, 6) A.S2p5W1CB_20Kmin(:, 6)];
% Argon
inp.S2p5W1CB.Ar.gammaVect  = gammaVect;
inp.S2p5W1CB.Ar.TGAcurve   = [B.S2p5W1CB_2p5Kmin_(1:end-1, 3) B.S2p5W1CB_5Kmin_(:, 3) B.S2p5W1CB_10Kmin_(1:end-3, 3) B.S2p5W1CB_20Kmin_(:, 3)];
inp.S2p5W1CB.Ar.TGAtemp    = [B.S2p5W1CB_2p5Kmin_(1:end-1, 1) B.S2p5W1CB_5Kmin_(:, 1) B.S2p5W1CB_10Kmin_(1:end-3, 1) B.S2p5W1CB_20Kmin_(:, 1)];
inp.S2p5W1CB.Ar.DTAcurve   = [B.S2p5W1CB_2p5Kmin_(1:end-1, 4) B.S2p5W1CB_5Kmin_(:, 4) B.S2p5W1CB_10Kmin_(1:end-3, 4) B.S2p5W1CB_20Kmin_(:, 4)];

% Gather Inputs and create matrices: SEBSma
% Air
inp.SEBSma.AIR.gammaVect = gammaVect;
inp.SEBSma.AIR.TGAcurve  = [A.SEBSma_2p5Kmin(:, 3) A.SEBSma_5Kmin(:, 3) A.SEBSma_10Kmin(:, 3) A.SEBSma_20Kmin(:, 3)]; % TGA Cur
inp.SEBSma.AIR.TGAtemp   = [A.SEBSma_2p5Kmin(:, 2) A.SEBSma_5Kmin(:, 2) A.SEBSma_10Kmin(:, 2) A.SEBSma_20Kmin(:, 2)]; % TGA Ter
inp.SEBSma.AIR.DTAcurve  = [A.SEBSma_2p5Kmin(:, 6) A.SEBSma_5Kmin(:, 6) A.SEBSma_10Kmin(:, 6) A.SEBSma_20Kmin(:, 6)]; % TGA Te
% Argon
inp.SEBSma.Ar.gammaVect  = gammaVect;
inp.SEBSma.Ar.TGAcurve   = [B.SEBSma_2p5Kmin_(:, 3) B.SEBSma_5Kmin_(1:end-1, 3) B.SEBSma_10Kmin_(1:end-1, 3) B.SEBSma_20Kmin_(1:end-1, 3)];
inp.SEBSma.Ar.TGAtemp    = [B.SEBSma_2p5Kmin_(:, 1) B.SEBSma_5Kmin_(1:end-1, 1) B.SEBSma_10Kmin_(1:end-1, 1) B.SEBSma_20Kmin_(1:end-1, 1)];
inp.SEBSma.Ar.DTAcurve   = [B.SEBSma_2p5Kmin_(:, 4) B.SEBSma_5Kmin_(1:end-1, 4) B.SEBSma_10Kmin_(1:end-1, 4) B.SEBSma_20Kmin_(1:end-1, 4)];

% Gather Inputs and create matrices: W1
% Air
inp.SASOL.AIR.gammaVect = gammaVect;
inp.SASOL.AIR.TGAcurve  = [A.SASOL_2p5Kmin(1:948, 3) A.SASOL_5Kmin(1:948, 3) A.SASOL_10Kmin(1:948, 3) A.SASOL_20Kmin(:, 3)]; %
inp.SASOL.AIR.TGAtemp   = [A.SASOL_2p5Kmin(1:948, 2) A.SASOL_5Kmin(1:948, 2) A.SASOL_10Kmin(1:948, 2) A.SASOL_20Kmin(:, 2)]; %
inp.SASOL.AIR.DTAcurve  = [A.SASOL_2p5Kmin(1:948, 6) A.SASOL_5Kmin(1:948, 6) A.SASOL_10Kmin(1:948, 6) A.SASOL_20Kmin(:, 6)]; %
% Argon
inp.SASOL.Ar.gammaVect  = gammaVect;
inp.SASOL.Ar.TGAcurve   = [B.SASOL_2p5Kmin_(1:end-1, 3) B.SASOL_5Kmin_(:, 3) B.SASOL_10Kmin_(1:end-2, 3) B.SASOL_20Kmin_(1:end-1, 3)];
inp.SASOL.Ar.TGAtemp    = [B.SASOL_2p5Kmin_(1:end-1, 1) B.SASOL_5Kmin_(:, 1) B.SASOL_10Kmin_(1:end-2, 1) B.SASOL_20Kmin_(1:end-1, 1)];
inp.SASOL.Ar.DTAcurve   = [B.SASOL_2p5Kmin_(1:end-1, 4) B.SASOL_5Kmin_(:, 4) B.SASOL_10Kmin_(1:end-2, 4) B.SASOL_20Kmin_(1:end-1, 4)];

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
