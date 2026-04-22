function [out_OFW] = OFW_isoMethod(gammaVect, tempConv_Vect, OFW_intTables, alphaTarget)

% FROM PRESENTAION: Osawa - FLynn - Wall Model, slides 33/34
%                   OFW is an iterative procedure that requires tables

% Extract OFW Integration Tables Vectors
redEffEnVect = OFW_intTables(:, 1);                                        % Array of Reduced Effective Activation Energy, [-]
aVect        = OFW_intTables(:, 2);                                        % Array of a Integration Parameter, [-]
bVect        = OFW_intTables(:, 3);                                        % Array of b Integration Parameter, [-]
b            = 0.457;                                                      % 1st Guess on b, from ASTM, [-]
                                    
% Universal Gas Constant                                    
R = 8.314;                                                                 % [J/(mol*K)]
                                    
% Evaluate Mid-Point (average) Heating Rate   
% At least you need 3 values of gamma
n_gamma = length(gammaVect);                                               % number of gamma used
gammaAVG = 1 / n_gamma * sum(gammaVect);                                   % Effective Mid-Point Heating Rate, [K/min]

% Evaluate Temperature at the Mid-Point Heating Rate
for ii = 1 : length(tempConv_Vect)
    
    tempConvMid(ii, :) = interp1(gammaVect, tempConv_Vect(ii, :), gammaAVG); % [K]
        % interpolates to find the temperature that corresponds to the average
        % (=midpoint) heate rate gamma since the reduce activation energy has
        % to be evaluated at that midpoint (iterative procedure needed)
    
end

% Compute the deltaLogGamma / deltaInvTemp 
% Evaluate LHS (left hand side = log(gamma_i)) for Linear Interpolation
LHS_interp = log10(gammaVect);

% Evaluate RHS for Linear Interpolation
        % you use the slope wrt 1/T
RHS_interp = 1 ./ tempConv_Vect;

% Linear Fitting for Ea
%   loop over the rows to perform the linear fitting between LHS and RHS;
%   extract then the slope to get Ea
for ii = 1 : lenght(tempConv_Vect)

    pol(ii, :) = polyfit(RHS_interp(ii, :), LHS_interp, 1);   % linear interpolation LHS and RHS
    m(ii) = - pol(ii,1); % you get the slope of the funcion = b*(Ea/R)

end

% Slope from the Linear Regression Considering Multiple Heating Rates
% Pt 4.3 in OFW-ASTM
m = m';

% Start loop over Ea
%     In each while loop you create a column vector of Ea (jj-th col) and stop
%     when tolerance is respected; then take the last column as the Ea that you
%     need. In the loop you look for a fittin b each time
jj = 1;
difference = 1;
toll = 1e-3;

while max(difference) > toll

    % Compute Activation Energy 
    Ea(:, jj) = R .* m ./ b';                                              %[J/mol] ---> m = b * Ea/R
    
    % Compute reduced activation energy
    %     you use the new w to get a new estimation of b from the tables of OFW
    %     as said, w is calculated at the midpoint heating rate temperauture
    redEffEnNew(:,jj) = Ea(:,jj) ./ (R .* tempConvMid);                    %[-] ---> w = Ea/(R*Tmid)

    % Find Closest Point to New Value for the Reduced Activation Energy and extract the b
    for ii = 1 : length(redEffEnNew(:,jj))

        [~, idxMin(ii)] = min(abs(redEffEnVect - redEffEnVect(ii,jj)));    %[-]
        b(ii) = bVect(idxMin(ii));
    end

    % Compute updated value for Ea with the new b
    Ea(:, jj) = R .* m ./ b';                                              %[J/mol]

    % State tollerance
    difference = abs(Ea(:,jj) - Ea(:, jj+1)) ./abs(Ea(:,jj));              %[-]
    % update counter
    jj = jj + 1;

end

% Extract Activation Energy 
%    (last column of the matrix obtained in the loop = last iteration value)
Ea_OFW = Ea(:, end);                                                       %[J/mol]

% Extract Multiplicative Factor for pre-exponential coefficient
b = b';                                                                    %[-]

for ii = 1 : length(b)
    a(ii) = aVect(idxMin(ii));
end

% Evaluate OFW pre-exponential coefficient
% see formula slide 34, but code giambelli has a minus bevore gammaAVG  WHY?? CHECK
A_OFW = - gammaAVG .* R ./ Ea_OFW' .* log(1 - alphaTarget) .* 10 .^ a;

% Evaluate Lifetime
gammaTarget = gammaVect(3);
initTemp    = 25 + 273.15;
tempVect    = tempConv_Vect(:, 3);
T_LFT       = 650 + 273.15;

tAlphaVect  = zeros(length(alphaTarget), 1);
gFunVect    = zeros(length(alphaTarget), 1);

funToInt = @(x) exp(- Ea_OFW(1) ./ (8.314 .* x));
dTime(1) = integral(funToInt, initTemp, tempVect(1)) / (gammaVect(3) * exp(- Ea_OFW(1) ./ (8.314 .* T_LFT)));
dg(1) = A_KAS(1) / (gammaVect(3)) * integral(funToInt, initTemp, tempVect(1));
idx   = 1;

while idx <= length(alphaTarget)
    
    for jdx = 2 : idx
        
        funToInt   = @(x) exp(- E_OFW(jdx) ./ (8.314 .* x));
        int        = integral(funToInt, tempVect(jdx-1), tempVect(jdx));
        dTime(jdx) = int / (gammaVect(3) * exp(- E_OFW(jdx) ./ (8.314 .* T_LFT)));
        dg(idx)    = A_OFW(jdx) / (gammaVect(3)) * int;
        
    end
    
    tAlphaVect(idx) = sum(dTime) * 60;
    gFunVect(idx)   = sum(dg);
    idx             = idx + 1;
    
end

g_alpha  = (A_OFW .* R) ./ (E_OFW .* exp(pol(:, 2)));
fFunVect = 1 ./ (diff(gFunVect) ./ diff(alphaTarget'));
% MANCA QUALCOSA CREDO -- CHECK

% output
out_OFW.Ea = Ea_OFW;
out_OFW.A  = A_OFW;
