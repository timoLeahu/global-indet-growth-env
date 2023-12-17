%%%% Simulation of the 3 dynamical system eqns dot{K}, dot{E}, dot{L}
% realization of two trajectories in the space starting form an arbitrary
% point and the fixed point P*. Code that generates the Fig 3 in report

clear; close all;
t0=0;
T=2500;
h=0.1;
%% parameters of the model
alpha = 0.1;
beta = 0.8; 
gamma = 0.58;
delta = 0.05;
epsil = 1;
eta = 1.5;
theta = 0.001;
E_bar = 0.17;
star_L = beta/(beta+epsil);
% repr. ratio f(L*)/f'(L*)
%f_ratio=(beta*epsil*eta)/((beta+epsil)*(eta*(beta+epsil)-(beta*epsil)));

%% equations
f = @(x) epsil*((1-x)^((epsil-eta*(1+epsil))/eta))*(x^(1-beta))/beta; 
diff_f = @(x) -(epsil*((1-x)^((epsil-eta*(1+epsil))/eta))...
    *(((beta+epsil)*eta-epsil)*x+(1-beta)*eta))/(beta*eta*(x-1)*(x^beta));

K_f = @(K, E, L) (K^alpha)*(L^(beta-1))*(E^gamma)*((beta+epsil)*L-beta)/epsil;
E_f = @(K, E, L) E*(E_bar-E)-delta*(K^alpha)*(L^beta)*(E^gamma);
L_f = @(K, E, L) (f(L)/diff_f(L))*(alpha*(K^(alpha-1))*(L^(beta-1))*(E^gamma)...
    *((beta+epsil)*L-beta)/epsil + gamma*(E_bar-E-delta*(K^alpha)*(L^beta)*(E^(gamma-1)))...
    +(theta-alpha*(K^(alpha-1))*(L^beta)*(E^gamma))/eta);

%% Runge-Kutta
%>>> initial point - set 1
Klv(1)=9.002396759;
Elv(1)=0.035860545540;
Llv(1)=0.4059105885;   
% Klv(1)=8.002396759;
% Elv(1)=0.025860545540;
% Llv(1)=0.3059105885; 

N=(T-t0)/h;% timesteps
tt= t0:h:T;
% Runge-Kutta coefficients calcualtion of steps
for n = 1:N
    K1k = K_f(Klv(n), Elv(n), Llv(n));
    K1e = E_f(Klv(n), Elv(n), Llv(n));
    K1l = L_f(Klv(n), Elv(n), Llv(n));
    K2k = K_f(Klv(n) + h*K1k/2, Elv(n) + h*K1e/2, Llv(n) + h*K1l/2);
    K2e = E_f(Klv(n) + h*K1k/2, Elv(n) + h*K1e/2, Llv(n) + h*K1l/2);
    K2l = L_f(Klv(n) + h*K1k/2, Elv(n) + h*K1e/2, Llv(n) + h*K1l/2);
    K3k = K_f(Klv(n) + h*K2k/2, Elv(n) + h*K2e/2, Llv(n) + h*K2l/2);
    K3e = E_f(Klv(n) + h*K2k/2, Elv(n) + h*K2e/2, Llv(n) + h*K2l/2);
    K3l = L_f(Klv(n) + h*K2k/2, Elv(n) + h*K2e/2, Llv(n) + h*K2l/2);
    K4k = K_f(Klv(n) + h*K3k, Elv(n) + h*K3e, Llv(n) + h*K3l);
    K4e = E_f(Klv(n) + h*K3k, Elv(n) + h*K3e, Llv(n) + h*K3l);
    K4l = L_f(Klv(n) + h*K3k, Elv(n) + h*K3e, Llv(n) + h*K3l);
    Klv(n+1) = Klv(n) + h*(K1k + 2*K2k + 2*K3k + K4k)/6;
    Elv(n+1) = Elv(n) + h*(K1e + 2*K2e + 2*K3e + K4e)/6;
    Llv(n+1) = Llv(n) + h*(K1l + 2*K2l + 2*K3l + K4l)/6;
end
%% graph plotting
figure(1)
hold on
plot(tt,Klv,'red')
plot(tt,Elv,'green')
plot(tt,Llv,'blue')
legend('K','E','L')
hold off

figure(2)
plot3(Klv,Elv,Llv)
hold on 
plot3(8.802396759,0.031860545540,0.4444,'Marker','*','Color','r') %P1*
plot3(13.11009169,0.0591165274446,0.4444,'o') %P2*
xlabel('K')
ylabel('E')
zlabel('L')

for j=1:2
    %% Runge-Kutta 2nd trajectory
    %>>> initial point - set 1
    % Klv(1)=3.1+j;
    % Elv(1)=0.005;
    % Llv(1)=0.22;   
    Klv(1)=8.802396759;
    Elv(1)=0.031860545540;
    Llv(1)=0.3745910588+j*0.01; 
    
    N=(T-t0)/h;
    tt= t0:h:T;
    % Runge-Kutta coefficients calcualtion of steps
    for n = 1:N
        K1k = K_f(Klv(n), Elv(n), Llv(n));
        K1e = E_f(Klv(n), Elv(n), Llv(n));
        K1l = L_f(Klv(n), Elv(n), Llv(n));
        K2k = K_f(Klv(n) + h*K1k/2, Elv(n) + h*K1e/2, Llv(n) + h*K1l/2);
        K2e = E_f(Klv(n) + h*K1k/2, Elv(n) + h*K1e/2, Llv(n) + h*K1l/2);
        K2l = L_f(Klv(n) + h*K1k/2, Elv(n) + h*K1e/2, Llv(n) + h*K1l/2);
        K3k = K_f(Klv(n) + h*K2k/2, Elv(n) + h*K2e/2, Llv(n) + h*K2l/2);
        K3e = E_f(Klv(n) + h*K2k/2, Elv(n) + h*K2e/2, Llv(n) + h*K2l/2);
        K3l = L_f(Klv(n) + h*K2k/2, Elv(n) + h*K2e/2, Llv(n) + h*K2l/2);
        K4k = K_f(Klv(n) + h*K3k, Elv(n) + h*K3e, Llv(n) + h*K3l);
        K4e = E_f(Klv(n) + h*K3k, Elv(n) + h*K3e, Llv(n) + h*K3l);
        K4l = L_f(Klv(n) + h*K3k, Elv(n) + h*K3e, Llv(n) + h*K3l);
        Klv(n+1) = Klv(n) + h*(K1k + 2*K2k + 2*K3k + K4k)/6;
        Elv(n+1) = Elv(n) + h*(K1e + 2*K2e + 2*K3e + K4e)/6;
        Llv(n+1) = Llv(n) + h*(K1l + 2*K2l + 2*K3l + K4l)/6;
    end
    
    plot3(Klv,Elv,Llv)
    
end
hold off