
H = linspace(700, 8000, n);
V = linspace(350, 880, n+1);
t1 = -50;%показывать на графике все времена больше t1
t2 = 200;%показывать на графике все времена меньшк t2
data = Tp;

[X, Y] = meshgrid(V, H);

figure;
scatter3(X, Y, data, 'filled');
xlabel('Скорость (V)');
ylabel('Высота (H)');
zlabel('время');
title('Врямя подъема в зависимости от скорости и высоты');
axis([-inf inf -inf inf t1 t2]); % Ограничение оси z от -10 до 100





H = linspace(700, 8000, n+1);
V = linspace(350, 880, n);

data = Tr;

[X, Y] = meshgrid(V, H);

figure;
scatter3(X, Y, data, 'filled');
xlabel('Скорость (V)');
ylabel('Высота (H)');
zlabel('время');
title('Время разгона в зависимости от скорости и высоты');
axis([-inf inf -inf inf t1 t2]); % Ограничение оси z от -10 до 100


H = linspace(700, 8000, n);
V = linspace(350, 880, n);

data = Tpr;

[X, Y] = meshgrid(V, H);

figure;
scatter3(X, Y, data, 'filled');
xlabel('Скорость (V)');
ylabel('Высота (H)');
zlabel('время');
title('Время подъем-разгона в зависимости от скорости и высоты');
axis([-inf inf -inf inf t1 t2]); % Ограничение оси z от -10 до 100

