aaa = 1;
function a = map_to_2pi(angle)
    a = mod(angle, 2*pi);
    if (a>pi)
        a = a-2*pi;
    end
end

function DKpendulumIntegrator(g, nu, a, gamma, theta1, theta2, omega1, omega2, tfinal)
close all
% DEFINE PARAMETERS

% set up start and end times of the integration
tstart = 0;
tend = tfinal;

% timestep 
dt = 0.005;

% initial condition 
x_0 = [theta1 theta2 omega1 omega2]; 
energy00 = -3*g
energy0pi = -g
energypi0 = g
energypipi = 3*g
CurrentEenergy = 0.5 * (2 * omega1^2 + omega2^2 + 2 * omega1 * omega2 * cos(theta2 - theta1)) - g* (2*cos(theta1) + cos(theta2))

% INTEGRATE AND PLOT
t = linspace(tstart,tend,floor((tend-tstart)/dt));

tic;

[t_out,x] = ode15s(@(t,x) odefcn(t,x,g, nu, a, gamma),t,x_0);

nmax = tfinal/dt

%Poincare section
iter=0;
thetasection = zeros(nmax,1);
thetasectionmapped = zeros(nmax,1);
velocitysection = zeros(nmax,1);

for ii = 1:nmax

    if (abs(map_to_2pi(x(ii,1))) < 0.01)
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


figure(1)
scatter(thetasectionmapped/(2*pi), velocitysection)
xlabel('$\theta_2/(2 \pi)$')
ylabel('$\omega_2$')
title('PoincarĂ© Section at $\theta_1 = 0$')
xlim([-1 1])
set(gca, 'linewidth', 2, 'fontsize', 22);
axis tight
axis auto

figure(2)
scatter(thetasection/(2*pi), velocitysection)
xlabel('$\theta_2/(2 \pi)$')
ylabel('$\omega_2$')
title('PoincarĂ© Section at $\theta_1 = 0$')
set(gca, 'linewidth', 2, 'fontsize', 22);
axis tight
axis auto

dft1 = fft(x(:,1));
power_spectra1 = abs(dft1).^2;
n = 0:1:nmax-1;
freq = n/nmax;

figure(3)
plot(freq(1:nmax/2),log(power_spectra1)(1:nmax/2))
xlabel('Frequency')
ylabel('log(Power Spectra)')
title('Power Spectra of $\theta_1$')
set(gca, 'linewidth', 2, 'fontsize', 20);

dft2 = fft(x(:,2));
power_spectra2 = abs(dft2).^2;
n = 0:1:nmax-1;
freq = n/nmax;

figure(4)
plot(freq(1:nmax/2),log(power_spectra2)(1:nmax/2))
xlabel('Frequency')
ylabel('log(Power Spectra)')
title('Power Spectra of $\theta_2$')
set(gca, 'linewidth', 2, 'fontsize', 20);

figure(5)
plot(t_out,x(:,1));
xlabel('Time')
ylabel('$\theta_1$')
title('Time Series of $\theta_1$')
set(gca, 'linewidth', 2, 'fontsize', 22);

figure(6)
plot(t_out,x(:,2));
xlabel('Time')
ylabel('$\theta_2$')
title('Time Series of $\theta_2$')
set(gca, 'linewidth', 2, 'fontsize', 22);

figure(7)
plot(x(:,1),x(:,3));
xlabel('$\theta_1$')
ylabel('$\omega_1$')
axis tight
set(gca, 'linewidth', 2, 'fontsize', 22);

figure(8)
plot(x(:,2),x(:,4));
xlabel('$\theta_2$')
ylabel('$\omega_2$')
axis tight
set(gca, 'linewidth', 2, 'fontsize', 22);

figure(9)
plot(x(:,1),x(:,2));
xlabel('$\theta_1$')
ylabel('$\theta_2$')
set(gca, 'linewidth', 2, 'fontsize', 22);

figure(10)
plot(t_out,-cos(x(:,1) + a*cos(nu*t_out)));
xlabel('Time')
ylabel('Height')
title('Height of bob 1')
set(gca, 'linewidth', 2, 'fontsize', 22);
ylim([-1-a, 1+a])


figure(11)
plot(t_out,-cos(x(:,2))-cos(x(:,1) + a*cos(nu*t_out)));
xlabel('Time')
ylabel('Height')
title('Height of bob 2')
set(gca, 'linewidth', 2, 'fontsize', 22);
ylim([-2-a, 2+a])


print(figure(1),'-dpdflatexstandalone','DoubleKapitzaPoincareMapped')
print(figure(2),'-dpdflatexstandalone','DoubleKapitzaPoincare')
print(figure(3),'-dpdflatexstandalone','DoubleKapitzaFourierTheta1')
print(figure(4),'-dpdflatexstandalone','DoubleKapitzaFourierTheta2')
print(figure(5),'-dpdflatexstandalone','DoubleKapitzaTimeSeriesTheta1')
print(figure(6),'-dpdflatexstandalone','DoubleKapitzaTimeSeriesTheta2')
print(figure(7),'-dpdflatexstandalone','DoubleKapitzaPhasePortrait1')
print(figure(8),'-dpdflatexstandalone','DoubleKapitzaPhasePortrait2')
print(figure(9),'-dpdflatexstandalone','DoubleKapitzaPhasePortraitTheta1vsTheta2')
print(figure(10),'-dpdflatexstandalone','DoubleKapitzaHeight1')
print(figure(11),'-dpdflatexstandalone','DoubleKapitzaHeight2')


system('pdflatex DoubleKapitzaPoincareMapped')
system('pdflatex DoubleKapitzaPoincare')
system('pdflatex DoubleKapitzaFourierTheta1')
system('pdflatex DoubleKapitzaFourierTheta2')
system('pdflatex DoubleKapitzaTimeSeriesTheta1')
system('pdflatex DoubleKapitzaTimeSeriesTheta2')
system('pdflatex DoubleKapitzaPhasePortrait1')
system('pdflatex DoubleKapitzaPhasePortrait2')
system('pdflatex DoubleKapitzaPhasePortraitTheta1vsTheta2')
system('pdflatex DoubleKapitzaHeight1')
system('pdflatex DoubleKapitzaHeight2')



system('rm *.log *.aux')
system('mv *.pdf Double_Pendulum_Kapitza/')
system('mv *.tex Double_Pendulum_Kapitza/')

end

% DEFINE ODE FUNCTION
function dxdt = odefcn(t, x, g, nu, a, gamma)
    dxdt = zeros(2,1);
    gPrime = g - a * nu^2 * cos(nu * t);
    dxdt(1) = x(3);
    dxdt(2) = x(4);
    dxdt(3) = -gPrime * (2*sin(x(1)) - sin(x(2))*cos(x(1)-x(2)))/(2 - cos(x(1) - x(2))^2) - sin(x(1) - x(2)) * (x(4)^2 + x(3)^2 * cos(x(1) - x(2)))/(2 - cos(x(1) - x(2))^2) - gamma*x(3);
    dxdt(4) = sin(x(1) - x(2)) * (2 * gPrime * cos(x(1)) + 2 * x(3)^2 + 2 * x(4)^2 * cos(x(1) - x(2)))/(2 - cos(x(1) - x(2))^2) - gamma*x(4);
end