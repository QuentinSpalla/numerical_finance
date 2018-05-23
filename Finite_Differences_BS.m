%% Prepare Matlab

clc
clf
clear all
close all

%% Parameters

% Black-Scholes parameters
S = 100
K = 100
sigma = 0.2;
r = 0.1;
T = 0.5;
t = 0;

% Notation
b = r-sigma^2/2;

% Underlying
X = log(S);
l = abs(b)*T + 3*sigma*sqrt(T); % localization
XL = X-l; % lower bound
XU = X+l; % upper bound
XN = 101

% Time
tL = 0; % lower bound
tU = T; % upper bound
tN = 50

%% Grid

% Underlying
Xg = linspace(XL,XU,XN);
Sg = exp(Xg);

% Time
tg = linspace(tL,tU,tN);

% Grid sizes
k = Xg(2)-Xg(1);
h = tg(2)-tg(1);

%% Scheme

% theta = 0; % explicit
theta = 0.5; % Crank-Nicolson
% theta = 1; % fully implicit

u = zeros(XN,tN);

% Terminal condition
for j=1:XN;
    u(j,tN) = max(Sg(j)-K,0);
end

% Coefficients
alpha = sigma^2/(2*k^2)-b/(2*k);
beta = -sigma^2/k^2-r;
gamma = sigma^2/(2*k^2)+b/(2*k);
Ak = diag(alpha*ones(XN-1,1),-1)+...
      diag(beta*ones(XN,1),0)+...
      diag(gamma*ones(XN-1,1),1);
A = eye(XN)-h*theta*Ak;
B = eye(XN)+h*(1-theta)*Ak;
v = zeros(XN,1);

ws = waitbar(0,'Running scheme: 0% done...');
for i = tN:(-1):2
    
    % Vector v
    v(1) = alpha*0;
    v(end) = gamma*(theta*(Sg(XN)*exp(k)-K*exp(-r*(T-tg(i-1)))) ...
             + (1-theta)*(Sg(XN)*exp(k)-K*exp(-r*(T-tg(i)))));
    
    % Systems
    u(:,i-1) = A\(B*u(:,i)+h*v);
    
    waitbar((tN-i)/(tN-1),ws,sprintf('Running scheme: %d%% done...',round((tN-i)/(tN-1)*100)));
    
end
close(ws)

% Plots
figure(1)
    mesh(T-tg,Sg,u)
    title('Finite Differences for Black-Scholes')
    xlabel('Time to maturity')
    ylabel('S_0')
    zlabel('Price')
figure(2)
    plot(Sg,u(:,1),Sg,BSCall(Sg,K,sigma,T,r))
    legend(['theta = ' num2str(theta)],'Theoretical')

%% Prices

j = find(Sg>=S,1); % corresponding j
i = find(tg>=t,1); % corresponding i

Sc = Sg(j); % closest to S
tc = tg(i); % closest to t

EstimatedC = u(j,i) % estimated price
TheoreticalC = BSCall(Sc,K,sigma,T-tc,r) % theoretical price
Error = TheoreticalC-EstimatedC % error
