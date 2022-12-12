aaa = 1;
function a = map_to_2pi(angle)
    a = mod(angle, 2*pi);
    if (a>pi)
        a = a-2*pi;
    end
end

function pendulumIntegrator(omega, theta1, theta2, omega1, omega2, tfinal, theta12, theta22, omega12, omega22)
close all
% DEFINE PARAMETERS

% set up start and end times of the integration
tstart = 0;
tend = tfinal;

% timestep 
dt = 0.01;

% initial condition 
x_0 = [theta1 theta2 omega1 omega2]; 
energy00 = -omega^2*3
energy0pi = -omega^2
energypi0 = omega^2
energypipi = 3*omega^2
CurrentEenergy = 0.5 * (2 * omega1^2 + omega2^2 + 2 * omega1 * omega2 * cos(theta2 - theta1)) - omega^2 * (2*cos(theta1) + cos(theta2))

% INTEGRATE AND PLOT
t = linspace(tstart,tend,floor((tend-tstart)/dt));

tic;

[t_out,x] = ode45(@(t,x) odefcn(t,x,omega),t,x_0);

%Second Pendulum
xSecond_0 = [theta12 theta22 omega12 omega22]; 
[t_out,xSecond] = ode45(@(t,xSecond) odefcn(t,xSecond,omega),t,xSecond_0);

nmax = tfinal/dt

%Poincare section
iter=0;
thetasection = zeros(nmax,1);
thetasectionmapped = zeros(nmax,1);
velocitysection = zeros(nmax,1);

for ii = 1:nmax
    if (abs(asin(sin(x(ii,1)))) < 0.01)
            thetasectionmapped(ii) = map_to_2pi(x(ii,2));
            thetasection(ii) = x(ii,2);
            velocitysection(ii) = x(ii,4);
            iter+=1;
    else
        thetasection(ii) = NaN;
        thetasectionmapped(ii) = NaN;
        velocitysection(ii) = NaN;
    end
end
iter
toc;

figure
scatter(thetasectionmapped/(2*pi), velocitysection)
xlabel('theta_2')
ylabel('omega_2')
set(gca, 'linewidth', 2, 'fontsize', 20);
xlim([-1 1])

figure
scatter(thetasection/(2*pi), velocitysection)
xlabel('theta_2')
ylabel('omega_2')
set(gca, 'linewidth', 2, 'fontsize', 20);
axis tight
axis auto

figure
scatter(t_out, abs(x(:,1) - xSecond(:,1)))
title('Lyaponov')
set(gca, 'linewidth', 2, 'fontsize', 20);

figure
scatter(t_out, abs(x(:,2) - xSecond(:,2)))
title('Lyaponov')
set(gca, 'linewidth', 2, 'fontsize', 20);

figure
scatter(t_out, log(abs(x(:,2) - xSecond(:,2))))
title('Lyaponov')
set(gca, 'linewidth', 2, 'fontsize', 20);

% dft1 = fft(x(:,1));
% power_spectra1 = abs(dft1).^2;
% n = 0:1:nmax-1;
% size(n)
% size(x(:,1))
% set(gca, 'linewidth', 2, 'fontsize', 20);
% freq = n/nmax;

% figure
% plot(freq(1:nmax/2),log(power_spectra1)(1:nmax/2))
% xlim([0 0.02])
% xlabel('Frequency')
% set(gca, 'linewidth', 2, 'fontsize', 20);
% ylabel('x_k')

% dft2 = fft(x(:,2));
% power_spectra2 = abs(dft2).^2;
% n = 0:1:nmax-1;
% size(n)
% size(x(:,1))
% freq = n/nmax;

% figure
% plot(freq(1:nmax/2),log(power_spectra2)(1:nmax/2))
% xlim([0 0.02])
% xlabel('Frequency')
% ylabel('x_k')
% set(gca, 'linewidth', 2, 'fontsize', 20);

% figure
% plot(t_out,x(:,1));
% xlabel('Time')
% ylabel('theta_1')
% set(gca, 'linewidth', 2, 'fontsize', 20);
% set(gca, 'linewidth', 2, 'fontsize', 20);

% figure
% plot(t_out,x(:,2));
% xlabel('Time')
% ylabel('theta_2')
% set(gca, 'linewidth', 2, 'fontsize', 20);
% set(gca, 'linewidth', 2, 'fontsize', 20);

% figure
% plot(x(:,1),x(:,3));
% xlabel('theta_1')
% ylabel('omega_1')
% set(gca, 'linewidth', 2, 'fontsize', 20);
% set(gca, 'linewidth', 2, 'fontsize', 20);

% figure
% plot(x(:,2),x(:,4));
% xlabel('theta_2')
% ylabel('omega_2')
% set(gca, 'linewidth', 2, 'fontsize', 20);
% set(gca, 'linewidth', 2, 'fontsize', 20);

figure
plot(x(:,1),x(:,2));
xlabel('theta_1')
ylabel('theta_2')
hold on;
plot(xSecond(:,1),xSecond(:,2))
set(gca, 'linewidth', 2, 'fontsize', 20);

end

% DEFINE ODE FUNCTION
function dxdt = odefcn(t,x,omega)
    dxdt = zeros(2,1);
    dxdt(1) = x(3);
    dxdt(2) = x(4);
    dxdt(3) = -omega^2 * (2*sin(x(1)) - sin(x(2))*cos(x(1)-x(2)))/(2 - cos(x(1) - x(2))^2) - sin(x(1) - x(2)) * (x(4)^2 + x(3)^2 * cos(x(1) - x(2)))/(2 - cos(x(1) - x(2))^2);
    dxdt(4) = sin(x(1) - x(2)) * (2 * omega^2 * cos(x(1)) + 2 * x(3)^2 + 2 * x(4)^2 * cos(x(1) - x(2)))/(2 - cos(x(1) - x(2))^2);
end