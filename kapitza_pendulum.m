aaa = 1;
function KapitzaPendulumIntegrator(g,nu,a,gamma, x0, v0, tfinal)
close all
tic
% DEFINE PARAMETERS

% set up start and end times of the integration
tstart = 0;
tend = tfinal;

% timestep 
dt = 0.001;

% initial condition 
x_0 = [x0 v0]; 
% INTEGRATE AND PLOT
t = linspace(tstart,tend,floor((tend-tstart)/dt));


[t_out,x] = ode15s(@(t,x) odefcn(t,x,g,nu,a,gamma),t,x_0);

nmax = tfinal/dt

dft = fft(x(:,1));
power_spectra = abs(dft).^2;
n = 0:1:nmax-1;
freq = n/nmax;
toc

figure(1)
plot(freq(1:nmax/2),log(power_spectra)(1:nmax/2))
xlim([0 0.02])
xlabel('Frequency')
ylabel('$|x_k|^2$')
set(gca, 'linewidth', 2, 'fontsize', 20);


figure(2)
plot(t_out,x(:,1)/pi);
xlabel('Time')
ylabel('$\theta/\pi$')
set(gca, 'linewidth', 2, 'fontsize', 22);

figure(3)
plot(t_out,cos(x(:,1))+a*cos(nu* t_out));
xlabel('Time')
ylabel('height')
set(gca, 'linewidth', 2, 'fontsize', 22);
ylim([-1-a 1+a])

figure(4)
plot(x(:,1),x(:,2));
xlabel('$\theta$')
ylabel('$\omega$')
set(gca, 'linewidth', 2, 'fontsize', 22);
axis equal


print(figure(1),'-dpdflatexstandalone','KapitzaFourier_Spectrum')
print(figure(2),'-dpdflatexstandalone','KapitzathetaTimeSeries')
print(figure(3),'-dpdflatexstandalone','KapitzaheightTimeSeries')
print(figure(4),'-dpdflatexstandalone','KapitzascatterPlot')

system('pdflatex KapitzaFourier_Spectrum')

system('pdflatex KapitzathetaTimeSeries')
system('pdflatex KapitzaheightTimeSeries')
system('pdflatex KapitzascatterPlot')

system('rm *.log *.aux')
system('mv *.pdf Kapitza/')
system('mv *.tex Kapitza/')

end

% DEFINE ODE FUNCTION
function dxdt = odefcn(t,x,g,nu,a,gamma)
    dxdt = zeros(2,1);
    dxdt(1) = x(2);
    dxdt(2) = g * sin(x(1)) - a*nu^2 * sin(x(1)) * cos(nu*t) - gamma * x(2);

end