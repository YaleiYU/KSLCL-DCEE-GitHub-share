
function [GradFa] = gradientF(theta, x)

theta_size = size(theta);

if  theta_size(1) == 1
    currentMean = mean(theta,1);
else currentMean = mean(theta);
end

predictedJ = reward(currentMean, x); % predicted mean of future reward based on current belief

individualJ = reward(theta, x);

GradFa = (individualJ - predictedJ);

% GradFa = (individualJ - predictedJ).*(-x.^2);
% according to the paper should be this GradFa =
% (individualJ-predictedJ)*phi(y)

end