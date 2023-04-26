clear ,clc
close all force


%__________________________________________________________________________
% условие задачи
n = 5; % количество элементарных разбиений
Vn = 350 * 1000 / 3600; % начальная скорость в м/с
Vk = 880 * 1000 / 3600; % конечная скорость в м/с
Hn = 700; % начальная высота
Hk = 8000; % конечная высота 
massa = 47000; % масса ЛА кг
S = 127; % площадь крыльев ЛА м^2
Vtek = Vn; % текущая скорость для элементарной области 
Htek = Hn; % текущая выоста для элементарной области 
dH = (Hk - Hn) / n; %изменение высоты в элементарной области
dV = (Vk - Vn) / n; %изменение скорости в элементарной области


%__________________________________________________________________________
%Разгон
for i = 1 : (n + 1) %"строка"
     for j = 1 : n  %"столбец"
       Vtek = Vn + j * dV;
       traz = Razgon(Htek, Vtek - dV, Vtek, massa, S);
       Tr(i, j) = traz;
     end 
     Htek = Hn + i * dH;
end

%Подъем
for i = 1 : n %"строка"
    Htek = Hn + i * dH;
    for j = 1 : (n + 1)  %"столбец"
        tnab = Podyem(Htek - dH, Htek, Vtek, massa, S);
        Vtek = Vn + j * dV;
        Tp(i, j) = tnab;
    end 
end

%Разгон+подъем
for i = 1 : n %"строка"
    Htek = Hn + i * dH;
    for j = 1 : n %"столбец"
        Vtek = Vn + j * dV;
        tpodraz = PodyemRazgon(Htek - dH, Htek, Vtek - dV, Vtek, massa, S);
        Tpr(i, j) = tpodraz;
    end 
end
%__________________________________________________________________________
% создание сетки минимального времени
s_time(n + 1, n + 1) = 0;
% 1 - разгон 2 - подем 3 - разгон-подъем 
s_direction(n + 1, n + 1) = 0;

% заполнение при Нmax
for j = n : - 1 : 1
    s_time(n + 1, j) = s_time(n + 1, j + 1) + Tr(n, j);
    s_direction(n + 1, j) = 1;
end

% заполнение при Vmax
for i = n : - 1 : 1
    s_time(i, n + 1) = s_time(i + 1, n + 1) + Tp(i, n + 1);
    s_direction(i, n + 1) = 2;
end

% заполнеие оставшихся значений
for i = n : - 1 : 1
    for j = n : - 1 : 1
        s_time(i, j) = s_time(i + 1, j) + Tp(i, j);
        s_direction(i, j) = 2;
        if (s_time(i, j + 1) + Tr(i, j) < s_time(i, j)) || (s_time(i + 1, j + 1) + Tpr(i, j) < s_time(i, j))
            if s_time(i + 1, j + 1) + Tpr(i, j) < s_time(i, j + 1) + Tr(i, j)
                s_time(i, j) = s_time(i + 1, j + 1) + Tpr(i, j);
                s_direction(i, j) = 3;
            else
                s_time(i, j) = s_time(i, j + 1) + Tr(i, j);
                s_direction(i, j) = 1;
            end
        end
    end
end
disp(s_time)

% построение графика оптимальных маневров
Vn = 350;
Vk = 880;
dV = (Vk - Vn) / n;
V_optimal(1) = Vn;
H_optimal(1) = Hn;
i = 0;
j = 0;
k = 2;
while s_direction(i + 1, j + 1) ~= 0
    switch s_direction(i + 1, j + 1)
        case 1 % разгон
            j = j + 1;
            V_optimal(k) = Vn + dV * j;
            H_optimal(k) = Hn + dH * i;
        case 2 % подъем
            i = i + 1;
            H_optimal(k) = Hn + dH * i;
            V_optimal(k) = Vn + dV * j;
        case 3 % разгон-подъем
            j = j + 1;
            i = i + 1;
            V_optimal(k) = Vn + dV * j;
            H_optimal(k) = Hn + dH * i;
    end
    k = k + 1;
end

figure(3);
hold all
ax = gca;
xlim([350 880]);
ylim([700 8000]);
ax.XTick = 350 : (880-350)/n : 880;
ax.YTick = 700 : (8000-700)/n: 8000;
xlabel('V');
ylabel('H');
plot(V_optimal, H_optimal, 'b-s');
grid on
%__________________________________________________________________________

disp('Разгон')
disp(Tr)
disp('Подъем')
disp(Tp)
disp('Разгон-Подъем')
disp(Tpr)
%__________________________________________________________________________

%Построение оптимальной траектории

%__________________________________________________________________________
function tRazgon = Razgon(H, V1, V2, massa, S) 
% определение времени в области разгон

    H = H; % Средняя высота
    a = func_a(H); % Скорость звука
    ro = func_ro(H); % Плотность
    V = (V1 + V2) / 2; % Средняя скорость
    M = V / a; % число маха
    P = func_P(M); % тяга двигателя в Н
    g = 9.8066; % Ускорение свободного падения м/c/c
    % Расчет угла атаки
    Cy0 = 0.1143;   % Cy = 0.0723*alpha+0.1143;
    Cya = 0.0723;
    alpha = (massa*g-Cy0*(ro*V^2)/2*S)/(P/57.3+Cya*(ro*V^2)/2*S);

    % коэффициент продольной силы
    Cx = func_cx(alpha);

    % время разгона
    tRazgon = abs((massa*(V2-V1))./(P.*cos(alpha)-Cx.*((V^2*ro)/2)*S)); 

end

% определение времени в области подъем
function tPodyem = Podyem(H1,H2,V,massa,S)

    H = (H1 + H2) / 2; % Средняя высота
    a = func_a(H); % Скорость звука
    ro = func_ro(H); % Плотность
    M = V / a; % число маха
    P = func_P(M); % тяга двигателя в Н
    g = 9.8066; % Ускорение свободного падения м/c/c

    %Расчет угла атаки
    Cy0 = 0.1143;   % Cy = 0.0723*alpha+0.1143;
    Cya = 0.0723;
    alpha = (massa*g-Cy0*(ro*V.^2)/2*S)./(P/57.3+Cya*(ro*V.^2)/2*S);
    % коэффициент продольной силы
    Cx = func_cx(alpha);

    teta = ((P - Cx * ro * V .^ 2 / 2 * S) * 57.3) / (massa * g);  

    % время подъема
    tPodyem = abs((57.3 * (H2 - H1) / (V * teta))); 



end

% определение времени в  области подъем-разгон
function tPodyemRazgon = PodyemRazgon(H1,H2,V1,V2,massa,S)

    H = (H1 + H2) / 2; % Средняя высота
    a = func_a(H); % Скорость звука
    ro = func_ro(H); % Плотность
    V = (V1 + V2) / 2; % Средняя скорость
    M = V / a; % число маха
    P = func_P(M); % тяга двигателя в Н
    g = 9.8066; % Ускорение свободного падения м/c/c

    % Расчет угла атаки
    Cy0 = 0.1143;   % Cy = 0.0723*alpha+0.1143;
    Cya = 0.0723;
    alpha = (massa*g-Cy0*(ro*V^2)/2*S)/(P/57.3+Cya*(ro*V^2)/2*S);

    % коэффициент продольной силы
    Cx = func_cx(alpha);
    k = (V1 - V2)/(H2 - H1);
    sin_teta = (P * cos(alpha) - Cx * ro * V .^2 / 2 * S) / (massa * (k * V + g));
    % время подъем-разгон
    tPodyemRazgon = abs(( 1 / (k * sin_teta) * log(V2 / V1))) ; 


end



%__________________________________________________________________________
% Полученные ранее полиномы
function res = func_ro(H)
res = -0.0001*H + 1.2247;
end

function res = func_a(H)
res = -0.0038*H + 340.2948;
end

function res = func_cy(alpha)
res =  -0.0002 *alpha.^3+0.0079*alpha.^2+0.0070*alpha+0.2001;
end

function res = func_cx(alpha)
res =  0.0004 *alpha.^3 - 0.0055*alpha.^2 + 0.0694*alpha + 0.0220;
end

function res = func_P(M)
res = 2.6146*10^3*M.^3 - 5.2416*10^3*M.^2 + 2.1611*10^3*M + 4.9096*10^3;
end 

