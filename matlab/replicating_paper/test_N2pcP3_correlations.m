%% FUNCTION test_N2pcP3_correlations
% ==> Test whether the larger correlation coefficient between the N2pc and
%   P3 in correct trials than in intrusion trials is statistically
%   significant
% ==> Inputs: 
%     r_cor -> corr coeff of correct trials
%     r_int -> corr coeff of intrusion trials
%     n_cor -> number of samples in correct trials
%     n_int -> number of samples in intrusion trials
% ==> Sources:
%     [for code start]                - https://de.mathworks.com/matlabcentral/fileexchange/44658-compare_correlation_coefficients 
%     [for underlying logic]          - https://core.ecu.edu/wuenschk/docs30/CompareCorrCoeff.pdf 
%     [for how 2 test spearman corrs] - https://luckytoilet.wordpress.com/2019/04/02/hypothesis-testing-for-difference-in-pearson-spearman-correlations/ 

function p = test_N2pcP3_correlations(correlationtype, r_cor,r_int,n_cor,n_int)
t_r1 = 0.5*log((1+r_cor)/(1-r_cor));
t_r2 = 0.5*log((1+r_int)/(1-r_int));
switch correlationtype
    case 'Pearson'
        z = (t_r1-t_r2)/sqrt(1/(n_cor-3)+1/(n_int-3));
    case 'Spearman'
        z = (t_r1-t_r2)/sqrt(1.06/(n_cor-3)+1.06/(n_int-3)); % 1.06 derived empirically by Fieller et al 1957 (https://www.jstor.org/stable/2332878)
    otherwise
        msgbox('correlationtype can only be Pearson or Spearman - cancelling!')
        return
end
p = (1-normcdf(abs(z),0,1))*2;
end