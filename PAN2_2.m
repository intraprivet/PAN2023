clc
clear variable

H = [500,1000,2000,3000,4000,5000,6000,7000,8000];
M =[0.35759	0.40814	0.47185 0.54851	0.62756	0.69344	0.78566];
alpha=[0,5,7.5,10,12.5,15,15,];

figure(2);
subplot(3,2,1);
plot(H,func_ro(H),'r')
title('График плотности воздушной среды от высоты')
legend ({'Исходные','ro(H)'},'Location','northeast');
grid on;
xlabel('Высота, м')
ylabel('Плотность воздушной среды, кг/м^3')

figure(2);
subplot(3,2,2);
plot(H,func_a(H),'r')
title('График зависимости скорости звука от высоты')
legend ({'Исходные','a(H)'},'Location','northeast');
grid on;
xlabel('Высота, м')
ylabel('Cкорость звука, м/c')

figure(2);
subplot(3,2,3);
plot(alpha,func_cy(alpha),'r')
title('График зависимости Сy от угла атаки')
legend ({'Исходные','Cy(alfa)'},'Location','northeast');
grid on;
xlabel('Угол атаки, град')
ylabel('Коэффициент Cy')


figure(2);
subplot(3,2,4);
plot(alpha,func_cx(alpha),'r')
title('График зависимости Сx от угла атаки')
legend ({'Исходные','Cx(alfa)'},'Location','northeast');
grid on;
xlabel('Угол атаки, град')
ylabel('Коэффициент Cx')

figure(2);
subplot(3,2,5);
plot(M,func_P(M),'r')
title('График зависимости тяги от скорости и высоты')
legend ({'Исходные','P(V,H)'},'Location','northeast');
grid on;
xlabel('Мах')
ylabel('Тяга, Н')

% Полученные ранее полиномы
function res = func_ro(H)
res = -0.0001*H + 1.2247;
end

function res = func_a(H)
res = -0.0038*H + 340.2948;
end

function res = func_cy(alpha)
res = 0.0723*alpha+0.1143;
end

function res = func_cx(alpha)
res =  0.0817*alpha - 0.0943;
end

function res = func_P(M)
res = 2.6146*10^3*M.^3 - 5.2416*10^3*M.^2 + 2.1611*10^3*M + 4.9096*10^3;
end 

