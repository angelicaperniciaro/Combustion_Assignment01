function out_ORT = ORT_isoMethod(gammaVect, TGA_spl, T_spl, deltaTGAEns, tempConv_Vect, alphaTarget, maxTGAEns)

% FROM PRESENTAION: Ortega Model, slides 37/38

% Universal Gas Constant                                    
R = 8.314;                                                                 % [J/(mol*K)]

% Set Integration Step
deltaAlpha = (alphaTarget(end) - alphaTarget(1)) / length(alphaTarget);

% check the incrememnt width: hp <= 0.02
if deltaAlpha > 0.02
    alphaTarget = [alphaTarget(1) : 0.02 : alphaTarget(end)];
end

% Loop over conversion degree array
for idx = 2 : length(alphaTarget)
    
    % alpha
    alpha            = alphaTarget(idx);
    % alpha - delta_alpha )alpha at previous index)
    alpha_deltaAlpha = alphaTarget(idx-1);
    
    deltaTGA                  = deltaTGAEns;                      
    conv_alpha                = deltaTGA .* alpha;                
    conv_alphaDeltaAlpha      = deltaTGA .* alpha_deltaAlpha;     
    TGA_conv_alpha            = maxTGAEns - conv_alpha;           
    TGA_conv_alphaDeltaAlpha  = maxTGAEns - conv_alphaDeltaAlpha; 
    
    % Link Conversion Degree with Temperature
    for jdx = 1 : size(TGA_spl, 1)
        deltaTempConv(jdx, :) = interp1(TGA_spl(jdx, :), T_spl(jdx, :), TGA_conv_alphaDeltaAlpha(jdx)) + 273.15;
        deltaTemp_tmp(jdx, :) = tempConv_Vect(:, jdx) - deltaTempConv(jdx, :); % T
    end
    
    deltaTemp(idx, :) = deltaTemp_tmp(:, idx)'; % Temperature Difference used at the denominator fo the LHS
    
end

% Linear Fitting for ORT Activation Energy
for ii = 1 : length(deltaTemp)
    % LHS for Linear Interpolation
    LHS_interp = log(gammaVect ./ deltaTemp(ii, :));
    
    % RHS for Linear Interpolation
    RHS_interp = 1 ./ tempConv_Vect(ii, :);
    
    % Slope coefficient
    p(ii, :) = polyfit(RHS_interp, LHS_interp, 1); % first order approx
end

% Activation energy
E_ORT = -p(:,1) .* R;                                                      % [J/mol]

% Assess the Reaction Models - g(alpha) for F1, A2, A3, A4
F1 = -log(1 - alphaTarget);                         % [-]
A2 = sqrt(-log(1 - alphaTarget));                   % [-]
A3 = (-log(1 - alphaTarget)) .^ (1 / 3);            % [-]
A4 = (-log(1 - alphaTarget)) .^ (1 / 4);            % [-]

% Linear Fitting for each reaction model
n1 = 1.89466100;
n2 = 3.63504095;
n3 = 1.001451033;

% Assess Mass Variation
TGAmaxRef = maxTGAEns(3);
massDecr = alphaTarget .* deltaTGARef;              % Mass Loss Percentage w.r.t Maximum
varFromMax = TGAmaxRef - massDecr;                  % Mass Loss Percentage on the TG curve 

% Extract Temperatures
TGACurveRef = TGA_spl(3, :);
tempRefVect = T_spl(3, :);
TAlpha = interp1(TGACurveRef, tempRefVect, varFromMax) + 273.15;           % Mass Loss Conversion Temperature

% Select Reference Heating Rate
gammaRef = gammVect(3);

% F1
pF1 = polyfit(1./TAlpha, log(F1 ./ TAlpha .^ n1)); % polyfitti una cosa lineare in 1/T, per cui la slope = p1(1) intercept = p1(2)
E_F1 = -pF1(1) * R / n3;
A_F1 = exp(pF1(2) - log(E_F1 / (gammaRef * R)) - n2 + n1*log(E_F1));

% A2 
pA2 = polyfit(1./TAlpha, log(A2 ./ TAlpha .^ n1)); 
E_A2 = -pA2(1) * R / n3;
A_A2 = exp(pA2(2) - log(E_A2 / (gammaRef * R)) - n2 + n1*log(E_A2));

% A3
pA3 = polyfit(1./TAlpha, log(A3 ./ TAlpha .^ n1)); 
E_A3 = -pA3(1) * R / n3;
A_A3 = exp(pA3(2) - log(E_A3 / (gammaRef * R)) - n2 + n1*log(E_A3));

% A4
pA4 = polyfit(1./TAlpha, log(A4 ./ TAlpha .^ n1)); 
E_A4 = -pA2(1) * R / n3;
A_A4 = exp(pA4(2) - log(E_A4 / (gammaRef * R)) - n2 + n1*log(E_A4));

% Build interpolative vectors
Evec = [E_F1 E_A2 E_A3 E_A4]; % dummy
Avec = [A_F1 A_A2 A_A3 A_A4]; % dummy

% Tang's Model Fitting
pTang = polyfit(Evec, log(Avec), 1);
A_ORT = exp(pTang(1) * E_ORT + pTang(2));

% Evaluate Lifetime
% gammaTarget = gammaVect(3);
initTemp    = 25 + 273.15;
tempVect    = tempConv_Vect(:, 3);
T_LFT       = 650 + 273.15;

tAlphaVect  = zeros(length(alphaTarget), 1);
gFunVect    = zeros(length(alphaTarget), 1);

funToInt = @(x) exp(- Ea_kas(1) ./ (8.314 .* x));
dTime(1) = integral(funToInt, initTemp, tempVect(1)) / (gammaVect(3) * exp(- Ea_KAS(1) ./ (8.314 .* T_LFT)));
dg(1) = A_KAS(1) / (gammaVect(3)) * integral(funToInt, initTemp, tempVect(1));
idx   = 1;

while idx <= length(alphaTarget)
    
    for jdx = 2 : idx
        
        funToInt   = @(x) exp(- E_KAS(jdx) ./ (8.314 .* x));
        int        = integral(funToInt, tempVect(jdx-1), tempVect(jdx));
        dTime(jdx) = int / (gammaVect(3) * exp(- E_KAS(jdx) ./ (8.314 .* T_LFT)));
        dg(idx)    = A_KAS(jdx) / (gammaVect(3)) * int;
        
    end
    
    tAlphaVect(idx) = sum(dTime) * 60;
    gFunVect(idx)   = sum(dg);
    idx             = idx + 1;
    
end

g_alpha  = (A_KAS .* R) ./ (E_KAS .* exp(pol(:, 2)));
fFunVect = 1 ./ (diff(gFunVect) ./ diff(alphaTarget'));
% MANCA UN PEZZO CHECK

% output
out_ORT.Ea = Ea_ORT;
out_ORT.A  = A_ORT;
