function [future_mean, future_variance] = predicted_mean_variance(theta, x, eta)

Ts=1; 
% [GradFa] = gradientF(theta, x); 
% % this is the F_i(k+1|k) 
% 
% future_thetaA = (1-eta.*GradFa).*theta;
% % theta(k+1|k) 
% 
% future_mean = mean(future_thetaA);
% % \bar{theta}(k+1|k) 
% 
% future_variance = var(future_thetaA);


%% improved version

[GradFa] = gradientF(theta, x);
% this is the F_i(k+1|k)


future_thetaA1 = theta(:,1) - eta.*GradFa.*(-x.^2)*Ts./(1+x.^4*eta);
% future_thetaA2 = theta(:,2) - eta.*GradFa.*(x)*Ts./(1+x.^2*eta);


future_thetaA = [future_thetaA1 ];

future_mean = [mean(future_thetaA1) ];

future_variance = var(future_thetaA);


end