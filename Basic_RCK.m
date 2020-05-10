%% Basic Ramsey-Cass-Koopmans Model 
% (Written by Kian Abbas Nejad)
% This code plots the stationary loci, saddle path to equilibrium and the
% vector field of a simple RCK model with CRRA utility
% $U(c)=\frac{c^{1-\sigma}}{1-\sigma}$ and Cobb-Douglas production function
% $F(k)=A k^a$. Solution of this model is given by the ODE: 
%
% $\frac{\dot{c}}{c}=\frac{1}{\sigma}(F'(k)-\delta-rho)$
%
% $\dot{k}=F(k)-\delta k - c$
%% Housekeeping
clearvars
close all
clc

%Default Parameters: A: Technology, a: Share of capital in production 
% sigma:risk-aversion parameter delta:depreciation rho:discount rate dx:
% maximum distance from equilibrium for the saddle path algorithm N: number
% of iterations from equilibrium
params={'A','a','sigma','delta','rho','dx','N'};
paramvals=[1,0.3,2,0.05,0.05,0.01,100];
printparam(params,paramvals);

% Prompt for Changing Parameters
paramsel = input('Override default parameters (y/n)? ','s');
if paramsel == 'y'
     for i=1:length(params)
         paramvals(i)=input(sprintf('Enter the value for %s:',params{i}));
     end   
     printparam(params,paramvals);
end
%% Stationary Loci
% Asking for min and max values of k and c
boundary = {'kmin', 'kmax'; 'cmin', 'cmax'};
boundaryval = [0 40; 0 2];
printparam(boundary,boundaryval);
boundsel = input('Override default boundary values (y/n)?','s');
if boundsel == 'y'
    for i=1:numel(boundary)
         boundaryval(i)=input(sprintf('Enter the value for %s:',boundary{i}));
    end 
     printparam(boundary,boundaryval)
end
k=boundaryval(1,1):0.1:boundaryval(1,2);

% Plotting 
fig=figure(1);
fig.Name='RCK_Phase_Diagram';
hold on
grid on
xlabel('$k$','FontSize',16,'interpreter','latex');
ylabel('$c$','FontSize',16,'interpreter','latex');
axis(reshape(boundaryval',1,4));
title('RBC Phase Diagram','interpreter','latex');
xline(kss(paramvals), 'LineWidth',2');
yline(css(kss(paramvals),paramvals),'--','LineWidth',1.2);
plot(k,css(k,paramvals),'LineWidth',2);

%% Saddle Path 

% Since the generic RCK equilibrium is known to have a unique stable
% manifold (also obvious from the vector field), I use the backward ODE to
% plot the stable manifold instead of other methods such as a shooting algorithm.
N=paramvals(end);
dx=paramvals(end-1);
% Defining the ODE as an anonymous function
dxdt = @(t,x) [F(x(1),paramvals)-paramvals(4).*x(1)-x(2) ; ...
    (x(2)./paramvals(3)).*(paramvals(1).*paramvals(2).*x(1).^(paramvals(2)-1)-paramvals(5)-paramvals(4))];

% creating a polygon/circle object with radius dx
% to be iterated through the manifold
cir = circle(kss(paramvals),css(kss(paramvals),paramvals),dx); 
plot(cir(1,:),cir(2,:),'r--')

for i=1:length(cir)
    [~,x]=ode45(dxdt, [0, -N],cir(:,i)); 
    plot(x(:,1),x(:,2),'r--','LineWidth',1.5) ;
end

%% Vector Field
[K,C]=meshgrid(linspace(0.2,boundaryval(1,2),20),linspace(0,boundaryval(2,2),20));
veck = zeros(size(K));  
vecc = zeros(size(C));

%Creating the Vector field and quasi-normalising the vectors 
for i = 1:numel(K)
    [k_dot,c_dot] = deriv(K(i), C(i),paramvals);
    nor= norm([k_dot; c_dot]).^0.6;
    veck(i) = k_dot./nor;
    vecc(i) = c_dot./nor;
end
%plotting the vector field
quiver(K,C,veck,vecc,0.5,'Color',[0.7 0.2 0],'MaxHeadSize',0.008,'AutoScale','on');
%% Functions
function [f, derf] = F(k, paramvals)
% Production function $F(k)=A k^a$ and its derivative
% Format: F(k,parameter value vector) 
A = paramvals(1);
a= paramvals(2);
f = A.*k.^a;
derf = A.*a.*k.^(a-1);        
end

function k_ss = kss(paramvals)
% Returns the Steady State value of k
% Format: kss(parameter value vector)
delta=paramvals(4);
rho=paramvals(5);
A=paramvals(1);
a = paramvals(2);
k_ss=((a.*A)./(rho+delta)) .^ (1/(1-a));
end

function c_ss = css(k,paramvals)
% Returns the Steady State Values of c 
% Format: css(k, parameter value vector) 
delta=paramvals(4);
c_ss= F(k,paramvals)-delta.*k;
end

function [k_dot,c_dot]=deriv(k,c,paramvals)
% The ODE system that returns kdot and cdot for each value of c and k
% Format: deriv(k,c, parmeter value vector)
delta=paramvals(4);
sigma=paramvals(3);
rho=paramvals(5);
k_dot=F(k,paramvals)-delta.*k-c;
[~,derf]=F(k,paramvals);
c_dot=(c./sigma).*(derf-rho-delta);
end
function printparam(params,paramvals)
% prints parameters and their values for the user
% Format: printparam(cell of parameter strings, parameter value vector)
n = numel(params);
disp('The parameters are:')
for i=1:n
    fprintf('%s: %0.2f \n',params{i},paramvals(i))
end
end

function cir = circle(x,y,r)
% CIRCLE gives the points of a polygon at (x,y) with radiur r
% used for iterating through the stable manifold
th = 0:pi/25:2*pi;
cir = [r * cos(th) + x;r * sin(th) + y];
end

