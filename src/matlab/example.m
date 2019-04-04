addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'build'));
rng('default') 
number_of_points = 50;
variance = 0.05;
theta = rand(number_of_points,1)*pi/4;
x = sin(theta) + 10 + randn(number_of_points,1)*variance;
y = cos(theta) - 10 + randn(number_of_points,1)*variance;
figure(1); clf;
scatter(x,y); %, 'LineStyle', 'none', 'Marker', '.')
hold on
plot([10,11],[-10,-10], 'Color', 'black', 'Marker', '*');
viscircles([10,-10], 1, 'Color', 'black');
rxy = circle_fit(x,y);
viscircles(rxy(2:3)', rxy(1), 'Color', 'red');
plot([rxy(2),rxy(2)+rxy(1)],[rxy(3),rxy(3)], 'Marker', 'diamond', 'Color', 'red');
rxy = circle_fit_taubin(x,y);
viscircles(rxy(2:3)', rxy(1), 'Color', 'green');
plot([rxy(2),rxy(2)+rxy(1)],[rxy(3),rxy(3)], 'Marker', 'square', 'Color', 'green');
xlim([8.5, 11.5]);
ylim([-11.5, -8.5]);
axis square
legend('Random samples','Original','LM','GWAF Taubin');
title("Estimate circle with AWGN \sigma^2="+variance);
print -dpng example.png
