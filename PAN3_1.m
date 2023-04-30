clear ,clc
close all force

%__________________________________________________________________________
% условие задачи

filter = true;%фильтр 
n = 50; % количество элементарных разбиений
Vn = 350 * 1000 / 3600; % начальная скорость в м/с
Vk = 880 * 1000 / 3600; % конечная скорость в м/с
Hn = 700; % начальная высота
Hk = 8000; % конечная высота 
massa = 47000; % масса ЛА кг
S = 127; % площадь крыла ЛА м^2
Vtek = Vn; % текущая скорость 
Htek = Hn; % текущая выоста  
dH = (Hk - Hn) / n; 
dV = (Vk - Vn) / n; 

%__________________________________________________________________________
%Разгон
Vtek = Vn; 
Htek = Hn;  
for i = 1 : (n + 1) %"строка"
     Vtek = Vn; 
     for j = 1 : n  %"столбец"
       Vtek = Vn + j * dV;
       traz = Razgon(Htek, Vtek - dV, Vtek, massa, S,filter);
       Tr(i, j) = traz;
     end 
     Htek = Hn + i * dH;
end

%Подъем
Vtek = Vn;
Htek = Hn;
for i = 1 : n %"строка"
    Vtek = Vn;      
    Htek = Hn + i * dH;
    for j = 1 : (n + 1)  %"столбец"
        tnab = Podyem(Htek - dH, Htek, Vtek, massa, S,filter);
        Vtek = Vn + j * dV;
        Tp(i, j) = tnab;
    end 
end

%Разгон+подъем 
Htek = Hn; 
for i = 1 : n %"строка"
    Vtek = Vn; 
    Htek = Hn + i * dH;
    for j = 1 : n %"столбец"
        Vtek = Vn + j * dV;
        tpodraz = PodyemRazgon(Htek - dH, Htek, Vtek - dV, Vtek, massa, S,filter);
        Tpr(i, j) = tpodraz;
    end 
end
%__________________________________________________________________________
 %disp('Разгон')
 % disp(Tr)
 % disp('Подъем')
 %disp(Tp)
 %disp('Разгон-Подъем')
 % disp(Tpr)
% %__________________________________________________________________________

%__________________________________________________________________________

% создание сетки минимального времени
s_time(n + 1, n + 1) = 0;
% 1 - разгон 2 - подем 3 - разгон-подъем 
s_direction(n + 1, n + 1) = 0;

% заполнение крайней верхней строки
for j = n : - 1 : 1
    s_time(n + 1, j) = s_time(n + 1, j + 1) + Tr(n, j);
    s_direction(n + 1, j) = 1;
end

% заполнение крайнего правого столбца
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
%disp('сеть времени')
%disp(s_time)

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

%Построение оптимальной траектории

%__________________________________________________________________________
function tRazgon = Razgon(H, V1, V2, massa, S,filter) 
% определение времени в области разгон

    H = H; % Средняя высота
    a = func_a(H); % Скорость звука
    ro = func_ro(H); % Плотность
    V = (V1 + V2) / 2; % Средняя скорость
    M = V / a; % число маха
    P = func_P(M); % тяга двигателя в Н
    g = calculate_gravity(H); % Ускорение свободного падения м/c/c
    % Расчет угла атаки
    Cy0 = 0.1143;   % Cy = 0.0723*alpha+0.1143;
    Cya = 0.0723;
    alpha = (massa*g-Cy0*(ro*V^2)/2*S)/(P/57.3+Cya*(ro*V^2)/2*S);

    % коэффициент продольной силы
    Cx = func_cx(alpha);

    % время разгона
    tRazgon =(massa*(V2-V1))./(P.*cos(alpha/57.3)-Cx.*((V^2*ro)/2)*S);
    if ((tRazgon< 0) && filter)
     tRazgon = 100000000; % запретная область
    end

end

% определение времени в области подъем
function tPodyem = Podyem(H1,H2,V,massa,S,filter)

    H = (H1 + H2) / 2; % Средняя высота
    a = func_a(H); % Скорость звука
    ro = func_ro(H); % Плотность
    M = V / a; % число маха
    P = func_P(M); % тяга двигателя в Н
    g = calculate_gravity(H); % Ускорение свободного падения м/c/c

    %Расчет угла атаки
    Cy0 = 0.1143;   % Cy = 0.0723*alpha+0.1143;
    Cya = 0.0723;
    alpha = (massa*g-Cy0*(ro*V.^2)/2*S)./(P/57.3+Cya*(ro*V.^2)/2*S);
    % коэффициент продольной силы
    Cx = func_cx(alpha);
    X = Cx * ro * V .^ 2 / 2 * S;

    teta = ((P - X) * 57.3) / (massa * g);
    % время подъема89
    tPodyem = ((57.3 * (H2 - H1) / (V * teta))); 
    if ((tPodyem<0) && filter)
    tPodyem = 100000000; % запретная область
    end

end

% определение времени в  области подъем-разгон
function tPodyemRazgon = PodyemRazgon(H1,H2,V1,V2,massa,S,filter)

    H = (H1 + H2) / 2; % Средняя высота
    a = func_a(H); % Скорость звука
    ro = func_ro(H); % Плотность
    V = (V1 + V2) / 2; % Средняя скорость
    M = V / a; % число маха
    P = func_P(M); % тяга двигателя в Н
    g = calculate_gravity(H); % Ускорение свободного падения м/c/c

    % Расчет угла атаки
    Cy0 = 0.1143;   % Cy = 0.0723*alpha+0.1143;
    Cya = 0.0723;
    alpha = (massa*g-Cy0*(ro*V^2)/2*S)/(P/57.3+Cya*(ro*V^2)/2*S);
    % коэффициент продольной силы
    Cx = func_cx(alpha);
    k = (V2 - V1)/(H2 - H1);
    X = Cx * ro * V .^2 / 2 * S; 
    teta = (P * cos(alpha/57.3) - X) / (massa * (k * V + g));
    % время подъем-разгон
    tPodyemRazgon = (( 1 / (k * teta) * log(V2 / V1))) ;
    if ((tPodyemRazgon < 0) && filter)
    tPodyemRazgon = 100000000; % запретная область
    end

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
res =   0.0004 *alpha.^2 + 0.0015*alpha + 0.0312;
end

function res = func_P(M)
res = 2*9.8066*(2.6146*10^3*M.^3 - 5.2416*10^3*M.^2 + 2.1611*10^3*M + 4.9096*10^3);
end 

function g = calculate_gravity(H)
    a3 = -3.4485e-14;
    a2 = 6.4749e-10;
    a1 = -4.4214e-06;
    a0 = 9.8051;
    
    g = a3 * H^3 + a2 * H^2 + a1 * H + a0;
end

