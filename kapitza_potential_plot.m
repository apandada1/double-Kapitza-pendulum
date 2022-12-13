n=5;

theta = linspace(0, 2*pi,200);
y = -cos(theta) + n^2 * sin(theta).^2;

figure(1)
plot(theta/pi, y, 'color', 'red', 'linewidth', 2)
xlabel('$\theta/\pi$')
title('Effective Potential')
yticks([])
set(gca, 'linewidth', 2, 'fontsize', 40);
axis tight;

print(figure(1),'-dpdflatexstandalone','EffectivePotential')
system('pdflatex EffectivePotential')
system('rm *.log *.aux EffectivePotential.tex EffectivePotential-inc.pdf')