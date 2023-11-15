clear; clc; close all;

%% Specify parameter set
alpha = 0.1;
beta = 0.8; 
gamma = 0.58;
delta = 0.05;
epsil = 1;
eta = 1.5;
theta = 0.001;
E_bar = 0.21;

N=200;
star_L = beta/(beta+epsil);

%% specify the region of the plot for vector plot
% x==K, y==E, z==L corresp in our model
[K, E] = meshgrid(1:1:20, 0.01:0.01:0.1);
% f = @(x) epsil*((1-x).^((epsil-eta*(1+epsil))/eta)).*(x.^(1-beta))/beta; 
% diff_f = @(x) -(epsil*((1-x).^((epsil-eta*(1+epsil))/eta))...
%     .*(((beta+epsil)*eta-epsil).*x+(1-beta)*eta))./(beta*eta*(x-1).*(x.^beta));
L=star_L;
Kdot = (K.^alpha).*(L.^(beta-1)).*(E.^gamma).*((beta+epsil)*L-beta)/epsil; %Note the use of .* and .^
Edot = E.*(E_bar-E)-delta*(K.^alpha).*(L.^beta).*(E.^gamma);
% Ldot = (f(L)./diff_f(L)).*(alpha*(K.^(alpha-1)).*(L.^(beta-1)).*(E.^gamma)...
%     .*((beta+epsil)*L-beta)/epsil + gamma*(E_bar-E-delta*(K.^alpha).*(L.^beta).*(E.^(gamma-1)))...
%     +(theta-alpha*(K.^(alpha-1)).*(L.^beta).*(E.^gamma))/eta);
%first plot the vector plot with quiver
figure
quiver(K,E,Kdot,Edot)



a=1; b=2; c=a*(b^2-1/(b^2));
f = @(t,x)[x(2)*x(3)-a*x(1)
    x(1)*(x(3)-c)-a*x(2)
    1-x(1)*x(2)];
t=linspace(0,50,1000);
[t,x]=ode45(f,t,[1 0 1]);

cm = jet(numel(t));

figure
hold all
for k = 1:size(x,1)-1
    hl = plot3([x(k,1),x(k+1,1)], [x(k,2),x(k+1,2)], [x(k,3),x(k+1,3)]);
    set(hl, 'LineStyle','-.', 'Color',cm(k,:));
end
hold off
grid on
view(50, 30)