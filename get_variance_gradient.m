
function [variance_gradient] = get_variance_gradient(theta, x, eta)

Ts =1;
delta_x = 0.01;
% delta_x denotes the size of step 

[GradFa] = gradientF(theta, x);

% future_thetaA = (1-eta.*GradFa).*theta;
% future_variance = var(future_thetaA) ;
%
%
% [GradFa_delta] = gradientF(theta, x+delta_x);
%
% future_thetaA_delta = (1-eta.*GradFa_delta).*theta;
% future_variance_delta = var(future_thetaA_delta) ;
%
%
% variance_gradient = (future_variance_delta-future_variance)./delta_x;

%%

future_thetaA1 = theta(:,1) - eta.*GradFa.*(-x.^2)*Ts./(1+x.^4*eta);
% future_thetaA2 = theta(:,2) - eta.*GradFa.*(x)*Ts./(1+x.^2*eta);


future_thetaA = [future_thetaA1];
future_mean = [mean(future_thetaA1)];
future_variance = var(future_thetaA) ;


[GradFa_delta] = gradientF(theta, x+delta_x);

future_thetaA1_delta = theta(:,1) - eta.*GradFa_delta.*(-x.^2)*Ts./(1+x.^4*eta);
% future_thetaA2_delta = theta(:,2) - eta.*GradFa_delta.*(x)*Ts./(1+x.^2*eta);


future_thetaA_delta = [future_thetaA1_delta ];
future_mean = [mean(future_thetaA1_delta) ];

future_variance_delta = var(future_thetaA_delta);

variance_gradient1 = (future_variance_delta-future_variance)./delta_x;
variance_gradient = sum(variance_gradient1);


end