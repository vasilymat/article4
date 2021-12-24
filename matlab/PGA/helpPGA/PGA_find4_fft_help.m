function [ f,s ] = PGA_find4_fft_help( T,sq,a_old,Rex,Global_constants,...
    TypeWindow,radius, Z, is_in_circle, Deg )

% В этом файле делаем процедуру разбиение на локальные радиусы для eig
% radius - это текущий радиус
% R - минимальное расстояние между центрами квадратов
% is_in_circle - передает текущую полную апертуру
% prev_phase - фаза, найденная на предыдущем этапе
% Rex - радиус окружности для исключения
% Deg - cтепень полинома zerfit


dx=Global_constants(4); 
L=Global_constants(5);
Rf = Global_constants(10); %полный радиус
nsq = 200;%оценочно - максимальное количество точек
l = size(T,2);

%the end of block of global constants

%формулы для собств. вектора
l2 = 64;
dx2 = dx*l/l2;
[x2,y2] = meshgrid(-L/2:dx2:L/2-dx2);

x2 = x2 + dx2/2; % Это нужно, чтобы отцентрировать is_in_circle
y2 = y2 + dx2/2;

Ps2 = winmultp(2); %Задаем множитель, который будт определять ширину окна
    function Psm = winmultp(sig)
        is_in_circle_t = circ(sqrt((x2).^2+(y2).^2)/(L/2*sq/l2));
        
        Pst = exp( -1*( (x2).^2+(y2).^2 )/(2*sig^2) )/sqrt(2*pi)/sig ;
        
        Psm = ifft2(fft2(fftshift(Pst)).*fft2(is_in_circle_t));
        m_psm = max(max(Psm));
        Psm = Psm/m_psm;
    end
        

lT = length(T);
tempz = false(lT);

 
 sqtype('full')
    function sqtype(typesq)
        if strcmp(typesq,'piece')
            td = 24;            
            tempz(td:128-td,td:128-td) = true;
        elseif strcmp(typesq,'full')
            td = 30;
            tempz(513-184+td:513+184-td,513-448+td:513+448-td) = true;
        end
    end

%local constants
Npie = 1;   %Это на такое количество отрезков делится размерность
rg = 32/1;  %r_grid - задаем радиус субобласти для eig
rg_o = 4;   %rg_oder, tau_o - это для третьего случая (type_eigpiece = 3) 
tau_o = 4;  %- на сколько элементов делится апертура по tau и r.

type_eigpiece = 3;

A=zeros(200,2);%Массив для настоящих максимумов, откуда будет вырезан квадратик

%the end of block of local constants
%------------------------------

da=fix(sq/2); %dELTa

% в этом куске формируются переменные для eig по частям

is_in_circle2 = circ(sqrt((x2).^2+(y2).^2)/Rf);
is_in_circle4 = circ(sqrt((x2).^2+(y2).^2)/(radius-0*dx2)*1);

[C, en, masl, Ntemp, exclude_is_in, Z2, prev_phase]...
    = eigpiece(type_eigpiece);
    function [C, en, masl, Ntemp, exclude_is_in, Z2, prev_phase]...
            = eigpiece(typep)
        exclude_is_in = circ(sqrt((x2).^2+(y2).^2)/Rex);
        is_in_circle3 = ~exclude_is_in & is_in_circle2; %это незадействованные ранее точки
        C = zeros( length(x2(is_in_circle3)).^2,1); %это строка, cодержащая всевозможные строки изображения
        en = zeros(Npie*Npie+1,1,'uint32');
        [Z2,~,~] = zerfit(x2,x2,y2,Rf,Deg); %тут мы сгенирировали фазу, которую
        prev_phase(is_in_circle2) = Z2*a_old; %нашли на предыдущем этапе
        
        switch typep 
          case 1
            masl = false(l2,l2,Npie*Npie);
            masch = zeros(l2/Npie,Npie,'int8');
            for ii2 = 1:Npie
                masch(:,ii2) = 64*(ii2-1)/Npie+1:64*ii2/Npie;
            end        
    
            for ii2 = 1:Npie    
                for jj = 1:Npie
                    Ntemp = (ii2-1)*Npie+jj;
                    masl(masch(:,ii2),masch(:,jj),Ntemp) = masl(...
                        masch(:,ii2),masch(:,jj),Ntemp) |...
                        is_in_circle3(masch(:,ii2),masch(:,jj));
        
                    %end - конец соответствующего участка
                    en(Ntemp+1) = size(x2(masl(:,:,Ntemp)),1);
         
                end
            end
            en = en.^2;
            en = cumsum(en);
            %конец eig по частям
          case 2 %!альтернативный eig по частям

            dots = is_in_circle2*0;
            sum_dots = dots;


            rs = fix(rg*sqrt(2)*0.99); % расстояние между центрами субобластей в sqrt2 больше радиуса
                           % сделано, чтобы не было пустых областей по
                           % диагонали

            deg = fix(l2/2/rg-10^-10); % это сколько вправо и влево точек сетки 0
                           % -10^-10  -для того, чтобы "вычесть границу"
                           % 32/32 должно быть меньше 1

            masl = false(l2,l2); %это битовый трехмерный массив с масками каждой из субобластей
            Ntemp = 0;
            for ii2 = -deg:deg
                for jj = -deg:deg
                    if (sqrt( (ii2*rs).^2 + (jj*rs).^2) <= l2/4) 
                        Ntemp = Ntemp + 1;
                        dots(33-ii2*rs,33-jj*rs) = true; %записываем 1 в центры областей
                        temp_circ = circ(sqrt((x2/dx2-ii2*rs).^2 ...
                             +(y2/dx2-jj*rs).^2)/rg);
                        temp_circ = temp_circ & is_in_circle2; 
                        masl(:,:,Ntemp) = temp_circ;
                        length(x2(temp_circ));
                        sum_dots = sum_dots + temp_circ; % это массив, на которую 
                                                         % потом отнормируем сум eig
                        en(Ntemp+1) = size(x2(temp_circ),1);
                    end
                end
            end
            sum_dots = sum_dots - 0.5; % тут избавляемся от деления на ноль
            sum_dots = abs(sum_dots);
            sum_dots = sum_dots + 0.5;

            en = en.^2;
            en = cumsum(en);
            %!конец альтернативный eig по частям
            
          case 3
            masl = false(l2,l2); %это битовый трехмерный массив с масками
                                 %каждой из субобластей
            [tau,rn] = cart2pol(x2,y2);
            
            rt = radius/rg_o;     %радиус центральной части
            dtau = 2*pi/tau_o; %шаг по tau
            
            temp_circ = rn <= rt;
            masl(:,:,1) = temp_circ;
            Ntemp = 1;
            en(Ntemp+1) = size(x2(temp_circ),1);
            
            for ii2 = 1:rg_o-1
                temp_circ2 = (rn <= radius*(1/rg_o+ii2/rg_o)) & ~temp_circ;
                for jj = 1:tau_o
                   Ntemp = Ntemp + 1;
                   temp_circ_tau = (tau >= -pi+(jj-1)*dtau) &...
                       (tau < -pi+jj*dtau);
                   masl(:,:,Ntemp) = temp_circ_tau & temp_circ2;
                   %figure(1); imagesc(masl(:,:,Ntemp)); colorbar ; colormap(gray);
                   
                   en(Ntemp+1) = size(x2(masl(:,:,Ntemp)),1);
                   
                end
                temp_circ = temp_circ2 | temp_circ;
            end
            
            sum_dots = temp_circ*0+1;
            en = en.^2;
            en = cumsum(en);
        case 4
            masl = false(l2,l2);
            Ntemp = 0;
            temp_circ = false(l2,l2);
            [tau,rn] = cart2pol(x2,y2);
%             rt = radius/rg_o; %Главный субрадиус
            for ii2 = 0:rg_o-1
                temp_circ2 = (rn <= radius*(1+ii2)/rg_o) & ~temp_circ;
                tau_o = 1 + ii2*2;
                dtau = 2*pi/tau_o; %шаг по tau
                for jj = 1:tau_o
                    Ntemp = Ntemp + 1;
                    temp_circ_tau = (tau >= -pi+(jj-1)*dtau) &...
                       (tau < -pi+jj*dtau);
                    masl(:,:,Ntemp) = temp_circ_tau & temp_circ2;
                    en(Ntemp+1) = size(x2(masl(:,:,Ntemp)),1);
                end
                temp_circ = temp_circ2 | temp_circ;
            end
            en = en.^2;
            en = cumsum(en);            
        end
        
    end

dat = da;
da = fix(da/1);
[x3,y3] = meshgrid(-1:1/(2*da):1-1/(2*da));
x3 = x3 + 1/(4*da); y3 = y3 + 1/(4*da);
exclude_new = circ(x3.^2+y3.^2);
da = dat;
      
mas_max = zeros(nsq,1);
mas_max(1) = 1;
N = 1;        
sT2 = abs(T);
P2 = logical(tempz);
stop_while = true;
G = zeros(l2,l2,200);
sh = 0;
da2 = l2/2;

while stop_while && (199 > N)

    [mas_max(N),P2,cx,cy] = cut_PGA3(P2,sT2,exclude_new,fix(da/1));
    if mas_max(1)/mas_max(N) > 2.5
        stop_while = false;
    end
    A(N,1) = cx(1);
    A(N,2) = cy(1);
    %P1 - массив комплексных чисел, представляет собой вырезанный квадратик
    %P2 - также массив комплексных чисел, содержит всевозможные вырезы

    %Квадратный вырез для окна l2
    %Умножаем на простой Гаусс по ширине окна
%     S3 = zeros(l2,l2);
    S3 = T(cx+sh-da2:cx+sh+da2-1,cy+sh-da2:cy+sh+da2-1).*Ps2;
%     S3 = S3.*Ps2; %Домножение на гауссов множитель
    G(:,:,N) = fftshift(fft2(fftshift(S3)));     
   
N = N + 1;
end
disp(sprintf('Количество точек = %g',N-1));

sum_um = zeros(size(G(:,:,1)));
sum_vm = zeros(size(G(:,:,1)));
for kk = 1:N-1
    Gt = G(:,:,kk);
    for ii = 1:Ntemp        
        xv = Gt(masl(:,:,ii));
        Ct = xv*xv';
        C(en(ii)+1:en(ii+1)) = C(en(ii)+1:en(ii+1)) + Ct(:);
    end
    sum_um = sum_um + Gt.*conj( circshift(Gt,[1 0]) );
    sum_vm = sum_vm + Gt.*conj( circshift(Gt,[0 1]) );
end

aP_x = angle(sum_um).*is_in_circle2;
aP_y = angle(sum_vm).*is_in_circle2;

f3 = r_hand3(aP_x(1:end-1,:),aP_y(:,1:end-1));

P = zeros(l2,l2);

en2 = diff(en);
en2 = int32(sqrt(single(en2)));
for ii = 1:Ntemp %!!! ВОт тут происходит обратный сбор
    dsize = en2(ii);
    Ct = zeros(dsize,dsize);
    Ct(true(dsize,dsize)) = C(en(ii)+1:en(ii+1));
    [V,~] = eig(Ct);
    P(masl(:,:,ii)) = V(:,end);
end

f = T*0;

aP = angle(P);
aP(exclude_is_in) = angle(exp(1i*prev_phase(exclude_is_in)));

[ masl_grad_x, masl_grad_y] = align_pieces( masl, is_in_circle2 );

a_end = grad_svd3( aP, l2, Z2, is_in_circle2, is_in_circle4,...
    masl_grad_x, masl_grad_y );
a_end(1:3) = 0;

f(is_in_circle) = Z*a_end;
% a_end2 = f3(is_in_circle2)\Z2;
% f(is_in_circle) = Z*a_end2';
s = a_end;
end
