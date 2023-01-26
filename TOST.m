function [TestResult, pValues, pooled_SD] = TOST(x1, x2, min_effectsize)
%% 
% calculated on basis of "two one sided tests" ("equivalence tests" doi: 10.1177/1948550617697177)
%
% INPUTS:
% x1, x2:           Dataarrays
% alpha:            Significancelevel
% min_effectsize:   Smallest effect size of interest (effect size: Cohens
% d: d = (mu_1 - mu_2)/(pooled std) ), qualitatively: how many STDs in
% difference between the means should be resolved


%% execute correct test
x1_m = nanmean(x1);
x2_m = nanmean(x2);
x1_sd= nanstd(x1);
x2_sd= nanstd(x2);
x1_l = length(x1(~isnan(x1)));
x2_l = length(x2(~isnan(x2)));
pooled_SD = (((x1_l-1)*x1_sd^2 + (x2_l - 1)*x2_sd^2)/(x1_l+x2_l - 2))^(0.5);
d_l = -min_effectsize*pooled_SD;
d_u = min_effectsize*pooled_SD;


%unequal variances
t_l = (x1_m - x2_m - d_l)/(x1_sd^2/x1_l + x2_sd^2/x2_l)^(0.5); %lower t
t_u = (d_u - (x1_m - x2_m))/(x1_sd^2/x1_l + x2_sd^2/x2_l)^(0.5); %upper t
df = (x1_sd^2/x1_l + x2_sd^2/x2_l)^2/((x1_sd^2/x1_l)^2/(x1_l-1) + (x2_sd^2/x2_l)^2/(x2_l - 1)); % satherthwaite correction: An approximate distribution of estimates of variance components.SATTERTHWAITE FEBiometrics. 1946 Dec; 2(6):110-4.
alpha = 0.05;
if t_u <= -tinv( 1-alpha ,df) && t_l >= tinv( 1-alpha ,df) %reject 0-hypothesis, that there is a true effect of "d_l" and "d_u" -> statistically equal
    TestResult = "Statistically Equal";
    pValues(1) = 1 - tcdf(t_u,df);
    pValues(2) = 1 - tcdf(t_l,df);
else
    TestResult = "Statistically Unequal";
    pValues(1) = 1 - tcdf(t_u,df);
    pValues(2) = 1 - tcdf(t_l,df);
end
%check CIs

% 
% 
% %%Example with x1~N(1,1) and x2~N(1.1, 1)
% NoRange = [10^3, 10^4, 10^5, 10^6];
% ColorInd = 0;
% Colors{1} = [1 0 0];Colors{2} = [0 1 0];Colors{3} = [0 0 1] ;Colors{4} = [1 1 0] ;
% figure
% for No =NoRange
%     ColorInd = ColorInd +1;
%     mu1 = 1;
%     mu2 = 1.1;
%     sigma1 = 1;
%     sigma2 = 1; %for simplicity
%     x1 = normrnd(mu1,sigma1, 1,No);
%     x2 = normrnd(mu2,sigma2, 1,No);
%     alpha = 0.03;
%     range = 0.01:0.001:0.2;
%     AllP = zeros(1,length(range));
%     ind = 0;
%     for min_effectsize = range
%         ind = ind +1;
%         [TestResult, p] = TOST(x1, x2, min_effectsize, alpha, 'unequal');
%         AllP(ind) =  max(p);
%     end
%     hold on
%     plot(range, AllP, '*', 'Color', Colors{ColorInd})
% end
% hold on
% line([abs(mu1-mu2)/((sigma1 + sigma2)/2) abs(mu1-mu2)/((sigma1 + sigma2)/2)], [0 1])
% ylabel('p-Value')
% xlabel('Minimal Effect Size')
% legend('Sample: 10^3', 'Sample: 10^4', 'Sample: 10^5', 'Sample: 10^6', '(\mu_1-\mu_2)/\sigma')
% ax = gca;
% ax.FontSize = 16;