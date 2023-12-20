% Code finds the fixed points of the main dynamical system 
% based on chosen parameters
clear

%% params
alpha = 0.1;
beta = 0.8; 
gamma = 0.58;
delta = 0.05;
epsil = 1;
eta = 1.5;
theta = 0.001;
E_bar = 0.17;

%% equations
syms x y z % respectively K E L
eqn1 = z-(beta/(beta+epsil)) == 0;
eqn2 = x-(alpha/(delta*theta))*y*(E_bar-y) == 0;
eqn3 = y+delta*(beta/(beta+epsil))^(beta/(1-alpha))*...
    (alpha/theta)^(alpha/(1-alpha))*y^((alpha+gamma-1)/(1-alpha))-E_bar == 0;

%% results
sol = solve([eqn1, eqn2, eqn3], [x, y, z],'ReturnConditions', true,'Real',true);
xSol = vpa(sol.x); % K
ySol = vpa(sol.y); % E
zSol = vpa(sol.z); % L
if size(xSol,1) == 2        % 2 fixed points 
    fprintf('--> P1*=(K,E,L)=(%4.8f, %4.8f, %4.8f)\n', xSol(1),ySol(1),zSol(1));
    fprintf('--> P2*=(K,E,L)=(%4.8f, %4.8f, %4.8f)\n', xSol(2),ySol(2),zSol(2));
elseif size(xSol,1) == 1    % 1 fixed point
    fprintf('--> P1*=(K,E,L)=(%4.8f, %4.8f, %4.8f)\n', xSol(1),ySol(1),zSol(1));
else % case the solver gives me also spurious sol or no real sol
    disp(xSol)
    disp(xSol)
    disp(xSol)
end