% computation of the Jacobian matrix of fixed points, J*, determinant of
% J*, trace of J* and other formula that characterize the stability
% coming from the Routh–Hurwitz Criterion 

clear
%% parameters
alpha = 0.1;
beta = 0.8; 
gamma = 0.58;
delta = 0.05;
epsil = 1;
eta = 1.5;
theta = 0.001;
% E_bar = 0.167159; % value which g'(E_bar)
E_bar = 0.17;
star_L = beta/(beta+epsil);

%% >> def fixed point P* with E_bar=0.21
% K=4.5625;
% E=0.0115;
% L=star_L;

%% >> def first fixed point P* with E_bar=0.17
% K=8.802396759;
% E=0.031860545540;
% L=star_L;
% >> def second fixed point P* with E_bar=0.17
K=13.11009169;
E=0.0591165274446;
L=star_L;

%% equations
f = @(x) epsil*((1-x)^((epsil-eta*(1+epsil))/eta))*(x^(1-beta))/beta; 
diff_f = @(x) -(epsil*((1-x)^((epsil-eta*(1+epsil))/eta))...
    *(((beta+epsil)*eta-epsil)*x+(1-beta)*eta))/(beta*eta*(x-1)*(x^beta));
dK_dL =(beta+epsil)/(delta*epsil)*E*(E_bar-E);
dE_dK =-delta*theta;
dE_dE =E_bar*(1-gamma)-E*(2-gamma);
dE_dL =-(beta+epsil)*E*(E_bar-E);
dL_dK =(f(L)/diff_f(L))*delta*theta*(theta*(1-alpha)/(alpha*eta*(E_bar-E))-gamma)/E;
dL_dE =(f(L)/diff_f(L))*gamma*((1-gamma)*(E_bar-E)-E-(theta/eta))/E;
dL_dL =(f(L)/diff_f(L))*(beta+epsil)*((theta*(beta+epsil)/epsil)-(theta/eta)-gamma*(E_bar-E));
E_A = (2-2*alpha-gamma)*((delta^(1-alpha))/((1-alpha)^(1-alpha)*(1-alpha-gamma)^(1-alpha-gamma))...
    *L^(beta)*(alpha/theta)^alpha)^(1/(2-2*alpha-gamma));
fprintf("\nE_bar = %f vs. E_A = %f\n", E_bar, E_A)

% >------------------------------------------------------
% compute the trace analytically - check of correctness
a1 = (eta*((1-gamma)*(beta+epsil)-beta*gamma*epsil)-beta*epsil*(1-gamma))/...
    (eta*(beta+epsil)-beta*epsil);
b1 = (beta*theta*(eta*(beta+epsil)-epsil))/(eta*(beta+epsil)-beta*epsil);
tr_J_anal =a1*(E_bar-E)-E+b1;
% >------------------------------------------------------

%% results:
J = [0 0 dK_dL; dE_dK dE_dE dE_dL; dL_dK dL_dE dL_dL]
det_J=det(J)
tr_J=trace(J) %OK
sigma_J=dE_dE*dL_dL-dE_dL*dL_dE-dK_dL*dL_dK
rho_J=-sigma_J*tr_J+det_J

% eigenvalues
eigenvalues = eig(J)
