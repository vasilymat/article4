function [ f,s ] = PGA_find4_fft( T,sq,tempz,Rex,Global_constants,...
    TypeWindow,radius, Z, is_in_circle, Deg, E2 )

% В этом файле делаем процедуру разбиение на локальные радиусы для eig
% radius - это текущий радиус
% R - минимальное расстояние между центрами квадратов
% is_in_circle - передает текущую полную апертуру
% prev_phase - фаза, найденная на предыдущем этапе
% Rex - радиус окружности для исключения

k=Global_constants(2);
dx=Global_constants(4); 
L=Global_constants(5);
Rf = Global_constants(10); %полный радиус
nsq = 250;%оценочно - максимальное количество точек
l = size(T,2);
f = zeros(size(T));
%the end of block of global constants

%формулы для собств. вектора
l2 = 64;
dx2 = dx*l/l2;
[x2,y2] = meshgrid(-L/2:dx2:L/2-dx2);

x2 = x2 + dx2/2; % Это нужно, чтобы отцентрировать is_in_circle
y2 = y2 + dx2/2;

td = 50;
tempz2 = false(1024,1024);
tempz2(513-184+td:513+184-td,513-448+td:513+448-td) = true;  
% tempz2 = false(1024,1024);
% tempz2(301:807,241:770) = true;

Ps2 = winmultp(2); %Задаем множитель, который будт определять ширину окна
    function Psm = winmultp(sig)
        sig = sig*sq/40;
        is_in_circle_t = circ(sqrt((x2-dx2/2).^2+(y2-dx2/2).^2)/(L/2*sq/l2));
        
        Pst = exp( -1*( (x2-dx2/2).^2+(y2-dx2/2).^2 )/(2*sig^2) )...
            /sqrt(2*pi)/sig ;
        
        Psm = ifft2(fft2(fftshift(Pst)).*fft2(is_in_circle_t));
        m_psm = max(max(Psm));
        Psm = Psm/m_psm;
        
        Pst = point_source(1,2*max(max(x2))/length(x2)/1.5,k,...
            x2,y2,dx2/2,dx2/2);  
        %Psm = Psm - Pst;
    end
 
%local constants

is_eig = true;
Npie = 1;   %Это на такое количество отрезков делится размерность
rg = 32/1;  %r_grid - задаем радиус субобласти для eig
rg_o = 5;   %rg_oder, tau_o - это для третьего случая (type_eigpiece = 3) 
tau_o = 4;  %- на сколько элементов делится апертура по tau и r.

type_eigpiece = 3;

A=zeros(150,2);%Массив для настоящих максимумов, откуда будет вырезан квадратик
%the end of block of local constants
%------------------------------

da=fix(sq/2); %dELTa

% в этом куске формируются переменные для eig по частям

is_in_circle2 = circ(sqrt((x2).^2+(y2).^2)/Rf);
is_in_circle4 = circ(sqrt((x2).^2+(y2).^2)/(radius-0*dx2)*1);

[C, en, masl, Ntemp, Z2] = eigpiece(type_eigpiece);
    function [C, en, masl, Ntemp, Z2] = eigpiece(typep)
        exclude_is_in = circ(sqrt((x2).^2+(y2).^2)/Rex);
        is_in_circle3 = ~exclude_is_in & is_in_circle2; %это незадействованные ранее точки
%         C = zeros( length(x2(is_in_circle3)).^2,1); %это строка, cодержащая всевозможные строки изображения
%         en = zeros(Npie*Npie+1,1,'uint32');
        [Z2,~,~] = zerfit(x2,x2,y2,Rf,Deg); %тут мы сгенирировали фазу, которую
        
        
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
            C = zeros( en(end),1);
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
            C = zeros( en(end),1);
            
        end
        
        
    end

dat = da;
da = fix(da/1);
[x3,y3] = meshgrid(-1:1/(2*da):1-1/(2*da));
x3 = x3 + 1/(4*da); y3 = y3 + 1/(4*da);
exclude_new = circ(x3.^2+y3.^2);
da = dat;

A = A*0;   
P2 = tempz2;
if sq > 26
    sT = abs(E2);
else
    sT = abs(T);
end
stop_while = true;
G = zeros(l2,l2,250);
sh = 0;
da2 = l2/2;
% while stop_while && (449 > N) && any(P2(:)) 
% 
%     [mas_max(N),P2,cx,cy] = cut_PGA3(P2,sT,exclude_new,fix(da/1));
% 
%     if mas_max(1)/mas_max(N) > 2.5
%         stop_while = false;
%     end
%     A(N,1) = cx(1);
%     A(N,2) = cy(1);
%     %P1 - массив комплексных чисел, представляет собой вырезанный квадратик
%     %P2 - также массив комплексных чисел, содержит всевозможные вырезы
% 
%     %Квадратный вырез для окна l2
%     %Умножаем на простой Гаусс по ширине окна
%     S3(:,:) = T(cx+sh-da2:cx+sh+da2-1,cy+sh-da2:cy+sh+da2-1).*Ps2;
%     
%     G(:,:,N) = fftshift(fft2(fftshift(S3)));
%     N = N + 1;
% end

    [Esort, Ecr] = sort(sT(P2),'descend');
    N = 1;
    N2 = 1;
    szx = 184*2-td*2+1;
    szy = 448*2-td*2+1;
    cx_s = 513-184+td-1; cy_s = 513-448+td-1;
    while stop_while && (259 > N)                

        if Esort(1)/Esort(N2) > 2
            stop_while = false;
        end
        [cx,cy] = ind2sub([szx,szy],Ecr(N2));                
        cx = cx + cx_s;
        cy = cy + cy_s;
        if P2(cx,cy) == true
            A(N,1) = cx;
            A(N,2) = cy;                    
            P2(cx-da*2:cx+da*2-1,cy-da*2:cy+da*2-1) =...
            P2(cx-da*2:cx+da*2-1,cy-da*2:cy+da*2-1) & ~exclude_new;
            S3 = T(cx+sh-da2:cx+sh+da2-1,cy+sh-da2:cy+sh+da2-1).*Ps2;
            G(:,:,N) = fftshift(fft2(fftshift(S3)));
            N = N + 1;
        end
        N2 = N2+1;

    end

fprintf('Количество точек = %g \n',N);


sum_um = zeros(size(G(:,:,1)));
sum_vm = zeros(size(G(:,:,1)));

if is_eig == true    
    C1 = zeros(en(2),1);
    C2 = zeros(en(3)-en(2),1);
    C3 = zeros(en(4)-en(3),1);
    C4 = zeros(en(5)-en(4),1);
    C5 = zeros(en(6)-en(5),1);
    C6 = zeros(en(7)-en(6),1);
    C7 = zeros(en(8)-en(7),1);
    C8 = zeros(en(9)-en(8),1);
    C9 = zeros(en(10)-en(9),1);
    C10 = zeros(en(11)-en(10),1);
    C11 = zeros(en(12)-en(11),1);
    C12 = zeros(en(13)-en(12),1);
    C13 = zeros(en(14)-en(13),1);
    C14 = zeros(en(15)-en(14),1);
    C15 = zeros(en(16)-en(15),1);
    C16 = zeros(en(17)-en(16),1);
    C17 = zeros(en(18)-en(17),1);

    for kk = 1:N-1
%         cx = A(kk,1);
%         cy = A(kk,2);
%         S3 = T(cx+sh-da2:cx+sh+da2-1,cy+sh-da2:cy+sh+da2-1).*Ps2;
%         Gt = exp(1i*angle(fftshift(fft2(fftshift(S3)))));
        Gt = exp(1i*angle(G(:,:,kk)));

        xv = Gt(masl(:,:,1));
        Ct_ = xv*xv';
        C1 = C1 + Ct_(:);  

        xv = Gt(masl(:,:,2));
        Ct_ = xv*xv';    
        C2 = C2 + Ct_(:);

        xv = Gt(masl(:,:,3));
        Ct_ = xv*xv';    
        C3 = C3 + Ct_(:);

        xv = Gt(masl(:,:,4));
        Ct_ = xv*xv';    
        C4 = C4 + Ct_(:);

        xv = Gt(masl(:,:,5));
        Ct_ = xv*xv';    
        C5 = C5 + Ct_(:);

        xv = Gt(masl(:,:,6));
        Ct_ = xv*xv';    
        C6 = C6 + Ct_(:);

        xv = Gt(masl(:,:,7));
        Ct_ = xv*xv';    
        C7 = C7 + Ct_(:);

        xv = Gt(masl(:,:,8));
        Ct2 = xv*xv';    
        C8 = C8 + Ct2(:);

        xv = Gt(masl(:,:,9));
        Ct_ = xv*xv';    
        C9 = C9 + Ct_(:);

        xv = Gt(masl(:,:,10));
        Ct_ = xv*xv';    
        C10 = C10 + Ct_(:);

        xv = Gt(masl(:,:,11));
        Ct_ = xv*xv';    
        C11 = C11 + Ct_(:);

        xv = Gt(masl(:,:,12));
        Ct2 = xv*xv';    
        C12 = C12 + Ct2(:);

        xv = Gt(masl(:,:,13));
        Ct_ = xv*xv';    
        C13 = C13 + Ct_(:);

        xv = Gt(masl(:,:,14));
        Ct_ = xv*xv';    
        C14 = C14 + Ct_(:);

        xv = Gt(masl(:,:,15));
        Ct_ = xv*xv';    
        C15 = C15 + Ct_(:);

        xv = Gt(masl(:,:,16));
        Ct_ = xv*xv';    
        C16 = C16 + Ct_(:);

        xv = Gt(masl(:,:,17));
        Ct_ = xv*xv';    
        C17 = C17 + Ct_(:);

    %     sum_um = sum_um + Gt.*conj( circshift(Gt,[1 0]) );
    %     sum_vm = sum_vm + Gt.*conj( circshift(Gt,[0 1]) );
    end

    C(1:en(2)) = C1(:);
    C(en(2)+1:en(3)) = C2(:);
    C(en(3)+1:en(4)) = C3(:);
    C(en(4)+1:en(5)) = C4(:);
    C(en(5)+1:en(6)) = C5(:);
    C(en(6)+1:en(7)) = C6(:);
    C(en(7)+1:en(8)) = C7(:);
    C(en(8)+1:en(9)) = C8(:);
    C(en(9)+1:en(10)) = C9(:);
    C(en(10)+1:en(11)) = C10(:);
    C(en(11)+1:en(12)) = C11(:);
    C(en(12)+1:en(13)) = C12(:);
    C(en(13)+1:en(14)) = C13(:);
    C(en(14)+1:en(15)) = C14(:);
    C(en(15)+1:en(16)) = C15(:);
    C(en(16)+1:en(17)) = C16(:);
    C(en(17)+1:en(18)) = C17(:);
    
else
    for kk = 1:N-1
        Gt = exp(1i*angle(G(:,:,kk)));
        sum_um = sum_um + Gt.*conj( circshift(Gt,[1 0]) );
        sum_vm = sum_vm + Gt.*conj( circshift(Gt,[0 1]) );
    end
    aP_x = angle(sum_um).*is_in_circle2;
    aP_y = angle(sum_vm).*is_in_circle2;

    f3 = r_hand3(aP_x(1:end-1,:),aP_y(:,1:end-1));
end


P = zeros(l2,l2);
    
en2 = diff(en);
en2 = int32(sqrt(single(en2)));

if is_eig == true 
    for ii = 1:Ntemp %!!! ВОт тут происходит обратный сбор
        dsize = en2(ii);
        Ct = zeros(dsize,dsize);
        Ct(true(dsize,dsize)) = C(en(ii)+1:en(ii+1));
        [V,~] = eig(Ct);
        P(masl(:,:,ii)) = V(:,end);
    end
    
    aP = angle(P);

    [ masl_grad_x, masl_grad_y] = align_pieces( masl, is_in_circle2 );

    a_end = grad_svd3( aP, l2, Z2, is_in_circle2, is_in_circle4,...
        masl_grad_x, masl_grad_y );
    a_end(1:3) = 0;

    f(is_in_circle) = Z*a_end;
else
    a_end = f3(is_in_circle2)\Z2;
    a_end(1:3) = 0;
    f(is_in_circle) = Z*a_end';
    
end
s = a_end;

end




