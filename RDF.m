clear
clc 
% Defaults for this blog post
width = 1.5;     % Width in inches
height = 3;    % Height in inches
alw = 0.75;    % AxesLineWidth
fsz = 11;      % Fontsize
lw = 2.3;      % LineWidth
msz = 8;       % MarkerSize


A = dlmread('SimulationData/RDF.dat');
B = dlmread('SimulationData/densityandenergy.dat');
number_density = mean(B((length(B)/2):end,2)) ; %%<- number density of argon at 85 K and 0.0001 GPa from the simulation

figure(1)
subplot(2,1,1)
pos = get(gcf, 'Position');
set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties
plot(A(1:end,1),A(1:end,2),'b-.','LineWidth',lw,'MarkerSize',msz);
axis([1 10.5 0 3])
legend('Atomic density = 0.0208 Å^{-3}, T = 85 K and P = 0.0001 GPa')
title('Pair correlation function of argon','FontSize',12,'FontWeight','bold','Color','k');
set(gca,'FontSize',9);
ylabel('g(r)','FontSize',9,'FontWeight','bold','Color','k');
xlabel('r (Å)','FontSize',9,'FontWeight','bold','Color','k');

subplot(2,1,2)
[S] = dftt(A(1:end,2), A(1:end,1), number_density);
pos = get(gcf, 'Position');
set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties
plot(A(1:end,1),S,'b-.','LineWidth',lw,'MarkerSize',msz);
axis([1 10.5 0 2.5])
title('Structure factor calculated from g(r) using Fourier transformation','FontSize',12,'FontWeight','bold','Color','k');
set(gca,'FontSize',9);
legend('Atomic density = 0.0208 Å^{-3}, T = 85 K and P = 0.0001 GPa')
ylabel('S(q)','FontSize',9,'FontWeight','bold','Color','k');
xlabel('q (Å^{-1})','FontSize',9,'FontWeight','bold','Color','k');
print('structure', '-dpng', '-r300');
print -depsc2 structure.eps
close;







