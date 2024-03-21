% This fits two regimes to all MSDs for relevant time lags for each
% molecule

% Load time vector `t` and `MSDs`
% Each column of `MSDs` should correspond to the MSD for a given molecule

N1 = 6; % Regime 1 from t = 1 to N1
N2 = 80; % Regime 2 from t = N2 to end

% Set up fittype and options.
ft = fittype( 'a*x^b', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.01 1];

coeffvals1 = cell(size(MSDs,2),1);
coeffvals2 = cell(size(MSDs,2),1);

for ii = 1:size(MSDs,2)

    % Fit model to data

    % Regime 1
    [xData, yData] = prepareCurveData( t(1:N1), MSDs(1:N1,ii) );
    [fitresult1, ~] = fit( xData, yData, ft, opts );   
    coeffvals1{ii} = coeffvalues(fitresult1);

    % Regime 2
    [xData, yData] = prepareCurveData( t(N2:end), MSDs(N2:end,ii) );
    [fitresult2, ~] = fit( xData, yData, ft, opts );   
    coeffvals2{ii} = coeffvalues(fitresult2);

end

Reg1 = cell2mat(coeffvals1); % Use column 2
Reg2 = cell2mat(coeffvals2); % Use column 2
