%>>> Code finds the fixed points of the main dynamical system 
% based on chosen parameters
clear

% params
alpha = 0.1;
beta = 0.8; 
gamma = 0.58;
delta = 0.05;
epsil = 1;
eta = 1.5;
theta = 0.001;
E_bar = 0.17;

syms x y z
eqn1 = z-(beta/(beta+epsil)) == 0;
eqn2 = x-(alpha/(delta*theta))*y*(E_bar-y) == 0;
eqn3 = y+delta*(beta/(beta+epsil))^(beta/(1-alpha))*...
    (alpha/theta)^(alpha/(1-alpha))*y^((alpha+gamma-1)/(1-alpha))-E_bar == 0;

% f = @(x) epsil*((1-x)^((epsil-eta*(1+epsil))/eta))*(x^(1-beta))/beta; 
% diff_f = @(x) -(epsil*((1-x)^((epsil-eta*(1+epsil))/eta))...
%     *(((beta+epsil)*eta-epsil)*x+(1-beta)*eta))/(beta*eta*(x-1)*(x^beta));

% eqn1 = (x^alpha)*(z^(beta-1))*(y^gamma)*((beta+epsil)*z-beta)/epsil == 0;
% eqn2 = y*(E_bar-y)-delta*(x^alpha)*(z^beta)*(y^gamma) == 0;
% eqn3 = (f(z)/diff_f(z))*(alpha*(x^(alpha-1))*(z^(beta-1))*(y^gamma)...
%     *((beta+epsil)*z-beta)/epsil + gamma*(E_bar-y-delta*(x^alpha)*(z^beta)*(y^(gamma-1)))...
%     +(theta-alpha*(x^(alpha-1))*(z^beta)*(y^gamma))/eta) == 0;


sol = solve([eqn1, eqn2, eqn3], [x, y, z],'ReturnConditions', true,'Real',true);
xSol = vpa(sol.x) % K
ySol = vpa(sol.y) % E
zSol = vpa(sol.z) % L
%K1 = xSol(1)