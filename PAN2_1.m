clc
clear variable

H = [500,1000,2000,3000,4000,5000,6000,7000,8000];
ro=[1.16727,1.11166,1.00655,0.909254,0.819347,0.736429,0.660111,0.590018,0.526783];
a=[338.370,336.435,332.532,328.584,324.589,320.545,316.452,312.306,308.105];
P=[5136.66081,5092.762,5024.47725,4958.63039,4868.39579,4734.26328,4646.46745];
Cx=[0.02,0.3,0.4,0.6,0.9,1.3];
Cy=[0.2,0.4,0.6,0.8,1.05,1.25];
M =[0.35759	0.40814	0.47185 0.54851	0.62756	0.69344	0.78566];
alpha=[0,5,7.5,10,12.5,15];


%Плотность воздушной среды
[koefroH,znro] = func_y(H,ro,1);
koefroH
figure(1);
subplot(3,2,1); 
plot(H,ro,'m',H,znro,'b--s')
title('График плотности воздушной среды от высоты')
legend ({'Апраксимированные','Исходные'},'Location','northeast');
grid on;
xlabel('Высота, м')
ylabel('Плотность воздушной среды, кг/м^3')


%Скорость звука
[koefaH,zna] = func_y(H,a,1);
koefaH
figure(1);
subplot(3,2,2); 
plot(H,a,'m',H,zna,'b--s')
title('График зависимости скорости звука от высоты')
legend ({'Апраксимированные','Исходные'},'Location','northeast');
grid on;
xlabel('Высота, м')
ylabel('Cкорость звука, м/c')

n=2;
% %Коэффициент подъемной силы
[koefCy,znCy] = func_y(alpha, Cy,1);
koefCy
figure(1);
subplot(3,2,3); 
plot(alpha,znCy,'m',alpha,Cy,'b--s')
title('График зависимости Сy от угла атаки')
legend ({'Апраксимированные','Исходные'},'Location','northeast');
grid on;
xlabel('Угол атаки, град')
ylabel('Коэффициент Cy')


% %Коэффициент силы лобового сопротивления
[koefCx,znCx] = func_y(alpha, Cx,1);
koefCx
figure(1);
subplot(3,2,4); 
plot(alpha,znCx,'m',alpha,Cx,'b--s')
title('График зависимости Сx от угла атаки')
legend ({'Апраксимированные','Исходные'},'Location','northeast');
grid on;
xlabel('Угол атаки, град')
ylabel('Коэффициент Cx')


%Тяга от скорости и высоты
[koefPM,znPM] = func_y(M, P,3);
koefPM
figure(1);
subplot(3,2,5); 
plot(M,znPM,'m',M,P,'b--s')
title('График зависимости тяги от скорости и высоты')
legend ({'Апраксимированные','Исходные'},'Location','northeast');
grid on;
xlabel('Мах')
ylabel('Тяга, Н')
%Создание функции

function [p1,y1] = func_y(y,x,z)

p1 = polyfit (y,x,z);

y1 = polyval(p1,y);
end

