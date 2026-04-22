function out = actEnergyEval(method, inp, miscInput, alphaTarget)

% EXTRACT QUANTITIES: inp
gammaVect    = inp.gammaVect; % [row vector 1x4]
TGAcurve     = inp.TGAcurve;
TGAtemp      = inp.TGAtemp;

% EXTRACT QUANTITIES: miscInput
OFW_intTables = miscInput.OFW_intTables;
NinterpPts    = miscInput.NinterpPts;
refTemp       = miscInput.refTemp;
T_LFT         = miscInput.T_LFT;
deltaAlpha    = miscInput.deltaAlpha;

% Curves Extraction and Interpolation: Temperature
% to increase the resulution of the trace by increasing the number of point in the array
TA = TGAtemp(:, 1);                                   % 1st Heating Program Temperature (2.5 K/min) Vector, [K]
TB = TGAtemp(:, 2);                                   % 2nd Heating Program Temperature (5 K/min) Vector, [K]
TC = TGAtemp(:, 3);                                   % 3rd Heating Program Temperature (10 K/min) Vector, [K]
TD = TGAtemp(:, 4);                                   % 4th Heating Program Temperature (20 K/min) Vector, [K]
TA_spl = linspace(TA(1), TA(end), NinterpPts);      % 1st Temperature Vector Definition, [K]
TB_spl = linspace(TB(1), TB(end), NinterpPts);      % 2nd Temperature Vector Definition, [K]
TC_spl = linspace(TC(1), TC(end), NinterpPts);      % 3rd Temperature Vector Definition, [K]
TD_spl = linspace(TD(1), TD(end), NinterpPts);      % 4th Temperature Vector Definition, [K]
T_spl  = [TA_spl; TB_spl; TC_spl; TD_spl];         % Create Temperature Matrix with higher resolution

%% START EVALUATION PROCEDURE

% Curves Extraction and Interpolation: TGA
% estrai le curve tg dell'heating program prendi gli estremi e fai la spline con i punti scelti per estrarre
TGAA = TGAcurve(:, 1);                             % 1st Heating Program TG Curve, [%]
TGAB = TGAcurve(:, 2);                             % 2nd Heating Program TG Curve, [%]
TGAC = TGAcurve(:, 3);                             % 3rd Heating Program TG Curve, [%]
TGAD = TGAcurve(:, 4);                             % 4th Heating Program TG Curve, [%]
TGAA_spl = spline(TA, TGAA, TA_spl);             % 1st Heating Program TG Curve Interpolation [%]
TGAB_spl = spline(TB, TGAB, TA_spl);             % 1st Heating Program TG Curve Interpolation [%]
TGAC_spl = spline(TC, TGAC, TA_spl);             % 1st Heating Program TG Curve Interpolation [%]
TGAD_spl = spline(TD, TGAD, TA_spl);             % 1st Heating Program TG Curve Interpolation [%]

% Extract Interest Values for each TGAx curve

% ------------------------------------
% FOR EACH HEATING RATE CASE: 
%     1) find a max value in the TGA, at the refTemp (100 C) since you may have
%        moisture removal to do before that reference temperature so do not take the first value of the trace
%     2) take the minimum
%     3) get the delta between max and min
%     4) do the conversion: product between delta and alphaTarget (to
%        get the increment vertically) for each conversion degree
%     5) compute mass loss percentage on the TG Curve: TGAA_max - convA 
%     6) do a backward interpolation to have the temperature associated to each
%        point of the TG curve
% NOTE: as the gamma (heating rate) increase, the TGA curve moves to higher temperatures
% ------------------------------------

%  Extract Interest Values for TGAA Curves
TGAA_max   = interp1(TA_spl, TGAA_spl, refTemp); % Mass Loss Percentage @ Onset tempRef, [%] 
TGAA_min   = min(TGAA_spl);                      % Mass Loss Percentage @ end of Experiment, [%]
deltaTGAA  = TGAA_max - TGAA_min;                % Maximum Delta TGA A, [%]
convA      = deltaTGAA .* alphaTarget;           % Conversion Degree Associated to the TGA @ Fixed Conversion Degree
TGAA_conv  = TGAA_max - convA;                   % Mass Loss Percentage on the TG curve @ Fixed Conversion Degree
tempA_conv = interp1(TGAA_spl, TA_spl, TGAA_conv) + 273.15; % Mass Loss Conversion Temperature, [K]

% Extract Interest Values for TGAB Curves
TGAB_max   = interp1(TB_spl, TGAB_spl, refTemp); % Mass Loss Percentage @ Onset tempRef, [%]
TGAB_min   = min(TGAB_spl);                      % Mass Loss Percentage @ end of Experiment, [%]
deltaTGAB  = TGAB_max - TGAB_min;                % Maximum Delta TGA B, [%]
convB      = deltaTGAB .* alphaTarget;           % Conversion Degree Associated to the TGA @ Fixed Conversion Degree
TGAB_conv  = TGAB_max - convB;                   % Mass Loss Percentage on the TG curve @ Fixed Conversion Degree
tempB_conv = interp1(TGAB_spl, TB_spl, TGAB_conv) + 273.15; % Mass Loss Conversion Temperature, [K]

% Extract Interest Values for TGAC Curves
TGAC_max   = interp1(TC_spl, TGAC_spl, refTemp); % Mass Loss Percentage @ Onset tempRef, [%]
TGAC_min   = min(TGAC_spl);                      % Mass Loss Percentage @ end of Experiment, [%]
deltaTGAC  = TGAC_max - TGAC_min;                % Maximum Delta TGA C, [%]
convC      = deltaTGAC .* alphaTarget;           % Conversion Degree Associated to the TGA @ Fixed Conversion Degree
TGAC_conv  = TGAC_max - convC;                   % Mass Loss Percentage on the TG curve @ Fixed Conversion Degree
tempC_conv = interp1(TGAC_spl, TC_spl, TGAC_conv) + 273.15; % Mass Loss Conversion Temperature, [K]

% Extract Interest Values for TGAD Curves
TGAD_max   = interp1(TD_spl, TGAD_spl, refTemp); % Mass Loss Percentage @ Onset tempRef, [%]
TGAD_min   = min(TGAD_spl);                      % Mass Loss Percentage @ end of Experiment, [%]
deltaTGAD  = TGAD_max - TGAD_min;                % Maximum Delta TGA C, [%]
convD      = deltaTGAD .* alphaTarget;           % Conversion Degree Associated to the TGA @ Fixed Conversion Degree
TGAD_conv  = TGAD_max - convD;                   % Mass Loss Percentage on the TG curve @ Fixed Conversion Degree
tempD_conv = interp1(TGAD_spl, TD_spl, TGAD_conv) + 273.15; % Mass Loss Conversion Temperature, [K]

% Gather Max Value for TGAs Curve
%   (each column is related to a specific heating rate)
maxTGAEns = [TGAA_max TGAB_max TGAC_max TGAD_max];
deltaTGAEns = [deltaTGAA deltaTGAB deltaTGAC deltaTGAD];

% Gather Isoconversional Temperature in a single matrix 
%   (each column is related to a specific heating rate)
tempConv_Vect = [tempA_conv' tempB_conv' tempC_conv' tempD_conv'];

%% CORE: ISOCONVERSIONAL METHOD EXECUTION

switch method
    
    case 'OFW'
        
        [out_OFW] = OFW_isoMethod(gammaVect, tempConv_Vect, OFW_intTables, alphaTarget);
        out = out_OFW;
        
    case 'KAS'
        
        [out_KAS] = KAS_isoMethod(gammaVect, tempConv_Vect, alphaTarget, deltaTGAC, TGAC_max, TGAC_spl, TC_spl);
        out = out_KAS;
        
    case 'ORT'
        
        % Nota: la riga seguente è parzialmente tagliata nello screen -- CHECK PARAMETRI DA REC
        [out_ORT] = ORT_isoMethod(gammaVect, TGA_spl, T_spl, deltaTGAEns, tempConv_Vect, alphaTarget, maxTGAEns);
        out = out_ORT;
        
end
