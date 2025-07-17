function [ J ] = reward(theta,y)


% x1 = x(1); 
% x2 = x(2);
theta1 = theta(:,1);

% theta2 = theta(:,2);
theta2 = 2;

% J = 2*x1 + 2*x2 - theta1.*x1^2 - theta2.*x2^2;

% J = 2*y - theta1*y^2; 
J = -theta1.*y.^2 + theta2.*y; 


end