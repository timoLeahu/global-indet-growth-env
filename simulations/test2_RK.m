clear; close all;
for j = 1:150
t0=0;
T=j;
h=0.01;
%parameters
alpha=0.45;
beta=0.1;
kappa=0.003;
gamma=0.3;
sigma=0.8;
mu=0.004;
delta=0.266;
p1=0.0004;
p2=0.0004;
%equations
P = @(p, w, x) (alpha-beta-kappa*x)*p;
W = @(p, w, x) gamma*p-sigma*w-mu*w;
X = @(p, w, x) mu*w-delta*x+p1*p*x+p2*w*x;
Tspec = @(pop) (mu*gamma*pop)/(delta*(sigma+mu)-(p1*(sigma+mu)+p2*gamma)*pop);
%Runge-Kutta
plv(1)=1300;
wlv(1)=0;
xlv(1)=0;
N=(T-t0)/h;
tt= t0:h:T;
%Runge-Kutta coefficients
for n = 1:N
K1p = P(plv(n), wlv(n), xlv(n));
K1w = W(plv(n), wlv(n), xlv(n));
K1x = X(plv(n), wlv(n), xlv(n));
K2p = P(plv(n) + h*K1p/2, wlv(n) + h*K1w/2, xlv(n) + h*K1x/2);
K2w = W(plv(n) + h*K1p/2, wlv(n) + h*K1w/2, xlv(n) + h*K1x/2);
K2x = X(plv(n) + h*K1p/2, wlv(n) + h*K1w/2, xlv(n) + h*K1x/2);
K3p = P(plv(n) + h*K2p/2, wlv(n) + h*K2w/2, xlv(n) + h*K2x/2);
K3w = W(plv(n) + h*K2p/2, wlv(n) + h*K2w/2, xlv(n) + h*K2x/2);
K3x = X(plv(n) + h*K2p/2, wlv(n) + h*K2w/2, xlv(n) + h*K2x/2);
K4p = P(plv(n) + h*K3p, wlv(n) + h*K3w, xlv(n) + h*K3x);
K4w = W(plv(n) + h*K3p, wlv(n) + h*K3w, xlv(n) + h*K3x);
K4x = X(plv(n) + h*K3p, wlv(n) + h*K3w, xlv(n) + h*K3x);
plv(n+1) = plv(n) + h*(K1p + 2*K2p + 2*K3p + K4p)/6;
wlv(n+1) = wlv(n) + h*(K1w + 2*K2w + 2*K3w + K4w)/6;
xlv(n+1) = xlv(n) + h*(K1x + 2*K2x + 2*K3x + K4x)/6;
end
%graph plotting
hold on
plot(tt,plv,'red')
plot(tt,wlv,'green')
plot(tt,xlv,'blue')
legend('Population','Waste','Toxicity')
%writing the results in csv file
res1=plv(N);    
res2=wlv(N);
res3=xlv(N);
RES=[res1 res2 res3];
%Tabe=[IDmesh,x(i),y(i),u_ex(i),x_chol(i),x_jac(i)];
% filename ="D:\mathematical engineering MSc\sem1\dynamic systems\test1_t_"+j+".csv";
% dlmwrite(filename,RES,’-append’);
hold off
end
figure(2)
plot3(plv,wlv,xlv); % works!
xlabel('Population')
ylabel('Waste')
zlabel('Toxicity')
grid on