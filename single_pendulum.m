aaa = 1;
function singlePendulumIntegrator(omega, x0, v0, x2, v2, tfinal)
close all
% DEFINE PARAMETERS

% set up start and end times of the integration
tstart = 0;
tend = tfinal;

% timestep 
dt = 0.001;

% initial condition 
x_0 = [x0 v0]; 
x_0second = [x2 v2];
% INTEGRATE AND PLOT
t = linspace(tstart,tend,floor((tend-tstart)/dt));


[t_out,x] = ode45(@(t,x) odefcn(t,x,omega),t,x_0);
[t_out,xsecond] = ode45(@(t,xsecond) odefcn(t,xsecond,omega),t,x_0second);

nmax = tfinal/dt

dft = fft(x(:,1));
power_spectra = abs(dft).^2;
n = 0:1:nmax-1;
size(n)
size(x(:,1))
freq = n/nmax;

figure
plot(freq(1:nmax/2),log(power_spectra)(1:nmax/2))
xlim([0 0.02])
xlabel('Frequency')
ylabel('x_k^2')

figure
plot(t_out,x(:,1));
xlabel('Time')
ylabel('x')
set(gca, 'linewidth', 2, 'fontsize', 20);


figure
plot(x(:,1),x(:,2));
xlabel('theta')
ylabel('omega')
set(gca, 'linewidth', 2, 'fontsize', 20);
axis equal

figure
plot(t_out,abs(x(:,1)-xsecond(:,1)))
xlabel('time')
ylabel('|x_1 - x_2|')

figure
scatter(t_out, log(abs(x(:,1) - xsecond(:,1))))
title('Lyaponov')
set(gca, 'linewidth', 2, 'fontsize', 20);
end

% DEFINE ODE FUNCTION
function dxdt = odefcn(t,x,omega)
    dxdt = zeros(2,1);
    dxdt(1) = x(2);
    dxdt(2) = -omega^2 * sin(x(1));
end