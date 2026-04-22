function [out_KAS] = KAS_isoMethod(gammaVect, tempConv_Vect, alphaTarget, deltaTGARef, TGAmaxRef, TGACurveRef, tempRefVect)

% FROM PRESENTAION: Kissinger - Akahira - Sunose Model, slides 35/36
                               
% Universal Gas Constant                                    
R = 8.314;                                                                 % [J/(mol*K)]
                                    
% Evaluate Mid-Point Heating Rate
% At least you need 3 values of gamma
n_gamma = length(gammaVect);                                               % number of gamma used
gammaAVG = 1 / n_gamma * sum(gammaVect);                                   % Effective Mid-Point Heating Rate, [K/min]

% Evaluate Temperature at the Mid-Point Heating Rate
for ii = 1 : length(tempConv_Vect)
    
    tempConvMid(ii, :) = interp1(gammaVect, tempConv_Vect(ii, :), gammaAVG); % [K]
        % interpolates to find the temperature that corresponds to the average
        % (= midpoint) heate rate gamma since the reduce activation energy has
        % to be evaluated at that midpoint (iterative procedure needed)
end

% Evaluate LHS for Linear Interpolation
%       LHS is gamma/T^2
LHS_interp = log(gammaVect ./ tempConv_Vect .^ 2);

% Evaluate RHS for Linear Interpolation
RHS_interp = 1 ./ tempConv_Vect;

% Linear Fitting for Ea
for ii = 1 : length(tempConv_Vect)
    
    pol(ii, :) = polyfit(RHS_interp(ii, :), LHS_interp(ii, :), 1);         % linear interpolation LHS and RHS
    m(ii)      = pol(ii, 1);                                               % you get the slope of the funcion = (Ea/R)
    
end

% Slope Coefficient
m = m'; 

% Activation Energy
E_KAS = -m .* R; % [J/mol]

% COMPENSATION EFFECT METHOD
% Assess the Reaction Models - g(alpha) for F1, A2, A3, A4
F1 = -log(1 - alphaTarget);                         % [-]
A2 = sqrt(-log(1 - alphaTarget));                   % [-]
A3 = (-log(1 - alphaTarget)) .^ (1 / 3);            % [-]
A4 = (-log(1 - alphaTarget)) .^ (1 / 4);            % [-]

% Select Reference Heating Rate in the vector
%   the most widely employed
gammaRef = gammaVect(3);

% Assess Mass Variation
massDecr = alphaTarget .* deltaTGARef;              % Mass Loss Percentage w.r.t Maximum
varFromMax = TGAmaxRef - massDecr;                  % Mass Loss Percentage on the TG curve 

% Extract Temperatures
TAlpha = interp1(TGACurveRef, tempRefVect, varFromMax) + 273.15;           % Mass Loss Conversion Temperature

% Linear Fitting for each reaction model
%       you get the dummy Ea and A to find a and b to finally get the real
%       A knowing the real Ea

% Tang's Model Coefficients
n1 = 1.89466100;
n2 = 3.63504095;
n3 = 1.001451033;

% F1
pF1 = polyfit(1./TAlpha, log(F1 ./ TAlpha .^ n1)); % polyfitti una cosa lineare in 1/T, per cui la slope è p1(1) e la q è p1(2)
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

% Build interpolative vectors with dummy Ea and dummy A
Evec = [E_F1 E_A2 E_A3 E_A4]; % dummy
Avec = [A_F1 A_A2 A_A3 A_A4]; % dummy

% Tang's Model Fitting
pTang = polyfit(Evec, log(Avec), 1); % 1st order linear fit
A_KAS = exp(pTang(1) * E_KAS + pTang(2)); % pTang(1) = slope = a;  pTang(2) = intercept = b

% Evaluate Lifetime
gammaTarget = gammaVect(3);
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
out_KAS.Ea = Ea_KAS;
out_KAS.A  = A_KAS;