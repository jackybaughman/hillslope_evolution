%% Hillslope evolution finite difference diffusion solution
% Computer Modeling, Feb 2016, Week 4, JSB

clear all
figure(1)
clf
%% Initialize 

% Constants and variables

Rho_r = 2700; % rock density kg/m3
Rho_s = 1600; % soil density kg/m3
kappa=0.003; %m^2/yr 
k=kappa*Rho_s; % efficientcy 
edot = 4 * (10^-5); % constant channel erosion rate m/year
wdotnaught = 2 * (10^-5); % m/year initial weathering rate
Hstar = 0.4; % m

% Arrays

% time array setup
dt = 100; % time step (years)
tmax = 10000 * 1000; % max time (ka)
t = 0:dt:tmax; % time array

% distance array setup
dx = 1; % distance step (m)
L = 500; % max dist (m) % min dist (m)
x = -L:dx:L; % dist array

% Initial mobile regolith
Hnaught = 4; % H at time zero 
H = Hnaught * ones(size(x)); %initial H

% Initial topography 
zbmax = 100; % max bedrock
S = zbmax/L; % Slope of initial bedrock
zb = zbmax-(S*abs(x)); % initial bedrock topo
z = zb+H; % Initial hill slope topo

% plot animation
n=100; %number of plots
tplot = tmax/n;
time = 0;


%% RUN

imax = length(t);
% In the time loop
for i = 1:imax;
time = time+dt; % accumulated time of model
wdot = wdotnaught*exp(-(H/Hstar)); % weathering flux exponential function
dzdx = diff(z)/dx; % slope
Q = dzdx*-k; % flux inputs
Q = [Q(1) Q Q(end)]; % filling in the Q array
dQdx = diff(Q)/dx; % rate of change of flux
dHdt = ((Rho_r/Rho_s)*wdot)-((1/Rho_s)*dQdx); % rate of change of mobile regolith
H = H + (dHdt*dt); % mobile regolith height
H = max(H,0); % cannot have negative mobile regolith
zb(2:end-1) = zb(2:end-1)-(wdot(2:end-1)*dt); % change of bedrock due to weathering
zb(1)=zb(1)-(edot*dt); % fixed boundary condition channel erosion
zb(end)=zb(end)-(edot*dt); % fixed boundary condition channel erosion
z = zb+H; % hillslope topography

if (rem(t(i),tplot)==0) % decrease plotting frequency - speed up animation
figure(3)
plot(x,z,'-.','linewidth',3) % plot hill slope over time
hold on
plot(x,zb,'r','linewidth',3) % plot bedrock over time 
hold off

   xlabel('Distance (m)','fontname','arial','fontsize',21) % x label
   ylabel('Height (m)','fontname','arial','fontsize',21) % y label
   set(gca,'fontsize',18,'fontname','arial') % axes number labels
   title(['Hillslope evolution after ',num2str((time/1000)-.1),'ka']) % title - accumulates model time
   axis([-L L 0 110]) % hold axes constant
   pause(.05)

end
end




 