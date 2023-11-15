function [ydot] = systemEqns(t,Y)
alpha = 0.1;
beta = 0.8; 
gamma = 0.58;
delta = 0.05;
epsil = 1;
eta = 1.5;
theta = 0.001;
E_bar = 0.21;
% P = [1-alpha; 2-beta; 3-gamma; 4-delta; 5-epsilon; 6-eta; 7-theta; 8-E_bar]
P = [alpha, beta, gamma, delta, epsil, eta, theta, E_bar];
f = @(x) P(5)*((1-x)^((P(5)-P(6)*(1+P(5)))/P(6)))*(x^(1-P(2)))/P(2); 

diff_f = @(y) -(P(5)*((1-y)^((P(5)-P(6)*(1+P(5)))/P(6)))...
    *(((P(5)+P(5))*P(6)-P(5))*y+(1-P(2))*P(6)))/(P(2)*P(6)*(y-1)*(y^P(2)));

Kdot = (Y(1)^P(1))*(Y(3)^(P(2)-1))*(Y(2)^P(3))*((P(2)+P(5))*Y(3)-P(2))/P(5);
Edot = Y(2).*(P(8)-Y(2))-P(4)*(Y(1).^P(1)).*(Y(3).^P(2)).*(Y(2).^P(3));
Ldot = (f(L)/diff_f(L)).*(P(1)*(Y(1).^(P(1)-1)).*(Y(3).^(P(2)-1)).*(Y(2).^P(3))...
    .*((P(2)+P(5))*Y(3)-P(2))/P(5) + P(3)*(P(8)-Y(2)-P(4)*(Y(1).^P(1)).*(Y(3).^P(2)).*(Y(2).^(P(3)-1)))...
    +(P(7)-P(1)*(Y(1).^(P(1)-1)).*(Y(3).^P(2)).*(Y(2).^P(3)))/P(6));

ydot = [Kdot;Edot;Ldot];
end 