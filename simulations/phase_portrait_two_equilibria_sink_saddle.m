%%%% Simulation of the 3 dynamical system eqns dot{K}, dot{E}, dot{L}
% realization of a series of trajectories in the space starting form 
% different initial values points arbitrary chosen
clear

t0=0;
T=2000; % time horizont
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

%% equations of dynamics
f = @(x) epsil*((1-x)^((epsil-eta*(1+epsil))/eta))*(x^(1-beta))/beta; 
diff_f = @(x) -(epsil*((1-x)^((epsil-eta*(1+epsil))/eta))...
    *(((beta+epsil)*eta-epsil)*x+(1-beta)*eta))/(beta*eta*(x-1)*(x^beta));

K_f = @(K, E, L) (K^alpha)*(L^(beta-1))*(E^gamma)*((beta+epsil)*L-beta)/epsil;
E_f = @(K, E, L) E*(E_bar-E)-delta*(K^alpha)*(L^beta)*(E^gamma);
L_f = @(K, E, L) (f(L)/diff_f(L))*(alpha*(K^(alpha-1))*(L^(beta-1))*(E^gamma)...
    *((beta+epsil)*L-beta)/epsil + gamma*(E_bar-E-delta*(K^alpha)*(L^beta)*(E^(gamma-1)))...
    +(theta-alpha*(K^(alpha-1))*(L^beta)*(E^gamma))/eta);
%% ploting
plot3(8.802396759,0.031860545540,0.4444,'o', 'MarkerFaceColor','black') %P1*
hold on
plot3(13.11009169,0.0591165274446,0.4444,'o', 'MarkerFaceColor','red') %P2*
xlabel('K')
ylabel('E')
zlabel('L')

%% Runge-Kutta iterations with different init points
N=(T-t0)/h;% timesteps
tt= t0:h:T;
for a=8.602396759:0.2:9.602396759
    for b=[0.028860545540 0.031860545540 0.033860545540]
        for c=[0.3559105885 0.37459105885 0.3959105885]
            % initial dyn values
            Klv(1)=a;
            Elv(1)=b;
            Llv(1)=c;
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
            hold on
        end %for c
    end %for b
end %for a

%% plot from another region of the space (close to P2*)
for a=16.602396759:0.2:17.602396759
    for b=[0.028860545540 0.036860545540 0.041860545540]
        for c=[0.3559105885 0.37459105885 0.3959105885]
            % initial dyn values
            Klv(1)=a;
            Elv(1)=b;
            Llv(1)=c;
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
        end %for c
    end %for b
end %for a
legend('Fixed point P_1^* (sink)','Fixed point P_2^* (saddle)')
hold off
