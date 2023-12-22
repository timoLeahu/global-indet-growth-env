%%%% Simulation of the 3 eqns dynamical system dot{K}, dot{E}, dot{L}
% realization of a trajectory in the space starting form a specific point.
% In this case I show an example of trajectory when the conditions of 
% Lemma 7 in the theory does not hold

clear; close all;
t0=0;
T=3000;
h=0.1;
%% parameters of the model
alpha = 0.19;
beta = 0.8; 
gamma = 0.8;
delta = 12.8620934285269;
epsil = 0.2;
eta = 4/9;
theta = 0.19;
E_bar = 10.75;
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
%>>> initial point
Elv(1)=0.25781852;%E2*
Klv(1)=0.21031403;%K2*
Llv(1)=0.815; %L_0   

N=(T-t0)/h;
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
plot3(Elv,Klv,Llv)
set(gca,'YDir','reverse')
hold on 
plot3(0.05444843,0.04527692,star_L,'o','Color','k','MarkerFaceColor','k') %P1*
hold on
plot3(0.25781852,0.21031403,star_L,'o','Color','r','MarkerFaceColor','r') %P2*
xlabel('E')
ylabel('K')
zlabel('L')
legend('Trajectory','Fixed point P_1^*','Fixed point P_2^*')
box on
hold off
