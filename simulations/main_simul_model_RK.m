clear; close all;
t0=0;
T=9000;
h=0.01;

%% parameters
alpha = 0.1;
beta = 0.8; 
gamma = 0.58;
delta = 0.05;
epsil = 1;
eta = 1.5;
theta = 0.001;
E_bar = 0.21;
star_L = beta/(beta+epsil);
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
Klv(1)=1.1;
Elv(1)=0.01;
Llv(1)=0.254;
N=(T-t0)/h;
tt= t0:h:T;
%Runge-Kutta coefficients
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
plot3(4.5625,0.0115,0.4444,'Marker','*','Color','k') %P1*
% plot3(21.24629,0.12504,0.4444,'*') %P2*
xlabel('K')
ylabel('E')
zlabel('L')

%% another trajectory starting from the fixed point
% Runge-Kutta
Klv(1)=4.5625;
Elv(1)=0.0115;
Llv(1)=0.4444;
% N=(T-t0)/h;
% tt= t0:h:T;
%Runge-Kutta coefficients
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

hold on
plot3(Klv,Elv,Llv)
hold off