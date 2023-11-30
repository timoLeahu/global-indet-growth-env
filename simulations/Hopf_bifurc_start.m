close all; clear

% code that find the zeros of a function with different initial values
% finds zeros arround the declared initial values

% params
alpha = 0.1;
beta = 0.8; 
gamma = 0.58;
delta = 0.05;
epsil = 1;
eta = 1.5;
theta = 0.001;
E_bar = 0.21;

star_L = beta/(beta+epsil);

lb = 0.01;             % Set a lower bound for the function.
ub = 0.1;          % Set an upper bound for the function.
x = NaN*ones(10,1);             % Initializes x.
f = @(x) x+delta*(beta/(beta+epsil))^(beta/(1-alpha))*...
    (alpha/theta)^(alpha/(1-alpha))*x^((alpha+gamma-1)/(1-alpha))-E_bar;
starting_points=linspace(lb,ub,10);
for i=1:10
        % Look for the zeros in the function's current window.
        x(i)=fzero(f, starting_points(i));
end

star_K=(alpha/(delta*theta))*x(1)*(E_bar-x(1));
disp(star_K)

