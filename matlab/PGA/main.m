% тут будет программа для hNa со всеми ухищрениями

%Global constants
%%
addpath('HelpFun','helpPGA','-end')

N = 1024;
Global_constants = zeros(10,1);
lam=0.000841; Global_constants(1)=lam; %Длина волны
k=2*pi/(lam); Global_constants(2)=k; % Wave number in invers millimeters
fok=(140)/2; Global_constants(3)=fok; % Расстояние, на котором предполагается искажающая фаза
dx=20*10^-3; Global_constants(4)=dx; % it is pixel size in mm
L=dx*N; Global_constants(5)=L; % Length of image area
[x,y]=meshgrid(-L/2:dx:L/2-dx);
x = x + dx/2 ; % Это нужно, чтобы отцентрировать is_in_circle
y = y + dx/2;
kz=rphase3(dx,L,k); % This is phase for angular spectrum method
radius = 9.32; % 9.32; % Это радиус для любек-файлов
a=896; b=368; % размеры исходного изображения
Global_constants(10) = radius;
temp1024 = zeros(1024,1024);
temp1024(513-b/2:513+b/2-1,513-a/2:513+a/2-1) = true;

Global_constants(6) = 513; %xc
Global_constants(7) = 513; %yc
Global_constants(8) = b; %dfx
Global_constants(9) = a; %dfy тут нужно будет исправить

hNa_div2 = true;
hNa = true;

%Генерируем полиномы Цернике
%%
[Z5full,aZ5,is_in_circle3] = zerfit(zeros(N,N),x,y,radius/2,2);
Dl = 128; 
l3 = Dl*2;
Es = zeros(l3,l3);
dx3 = dx*1024/l3;
[x3,y3] = meshgrid(-L/2:dx3:L/2-dx3);
x3 = x3 + dx3/2; % Это нужно, чтобы отцентрировать is_in_circle
y3 = y3 + dx3/2;

ind = 4;
temp_little = false(l3,l3,zdim);
temp_little(1+ind:end-ind,1+ind:end-ind,:) = true;
Se = zeros(l,l);

[Z5,~,is_in_circle2_part] = zerfit(x3*0,x3,y3,radius/2,2);

l2 = 64;
dx2 = dx*l/l2;
[x2,y2] = meshgrid(-L/2:dx2:L/2-dx2);

x2 = x2 + dx2/2; % Это для центрального цикла
y2 = y2 + dx2/2;
Rd2 = radius/2; %Radius div 2
%Генерируем необходимые полиносы Цернике
[Z2,~,is_in_circle2] = zerfit(x2,x2,y2,Rd2,5); %тут мы сгенирировали фазу, которую
[Zmin,amin,~] = zerfit(x*0,x,y,radius/2,5); %Это для hNa/2
Deg_hNa = 4;
[Z,~,is_in_circle] = zerfit(x*0,x,y,radius,Deg_hNa);


%%
%Выбрать, какую часть изображения искать. Скорее всего будет удалено.
%%
typesq = 'full';
tempz = false(N,N);
switch typesq
    case char('piece128')
        td = 24;
        tempz(td:128-td,td:128-td) = true;
    case char('E256')
        td = 20;
        tempz(513-128+td:513+127-td,513-128+td:513+127-td) = true;
    case char('full')
        td = 0;
        tempz(513-184+td:513+184-td,513-448+td:513+448-td) = true;
    case char('full512')
        td = 18;
        tempz(513-184/2+td:513+184/2-td,513-448/2+td:513+448/2-td) = true;
    case char('USAF')        
        tempz(301:807,241:770) = true;

end
%%  

   
% Тут мы пытаемся восстановить Na/2
tic
if hNa_div2 == true
  
    %%Компенсируем дефокус    
    is_in_circle_r = LPF( is_in_circle3, 30 );
    U = ifft2( fft2(E).*fftshift(LPF( is_in_circle3, 30 )) );
    xc = 0; yc = 0;     
    
    Es(1+ind:end-ind,1+ind:end-ind,:) =...
          U(513+xc-Dl+ind:513+xc+Dl-1-ind,513+yc-Dl+ind:513+yc+Dl-1-ind,:);
    u = fft2(Es);       
        
    [ aZ5,~,~] =...
      coefficient_finding_defocus(u,temp_little,Z5,aZ5,is_in_circle2_part);    
    Se(is_in_circle3) = Z5full*aZ5;
    u = fft2(U).*fftshift(exp(-1i*Se));
    E2 = ifft2( u );    
    %%
    
    TypeWindow = 8;
    fft_method = 1;
    Deg = 5;
    fprintf('Энтропия исходного изображения = %g \n',...
                                                shann_entropy(E2(tempz)) );    
    N2 = 5;
    Gs_temp = zeros(l,l);
    Gs = repmat(Gs_temp,[1 1 N2+1]);
    Sh0 = zeros(N2,1);
    Q = [40,36,32,24,20,16,8];    

    for t2 = 1:N2
        fprintf('цикл %g \n',t2);
        fprintf('окно = %g \n',Q(t2));
        
%         Global_constants(10) = radius/2;                
%         [Se,amin] = PGA_find4_fft_help( E2,Q(t2),amin,...
%          0,Global_constants,TypeWindow,radius/2,Zmin, is_in_circle3, Deg );
        %% Тут вставляем главную функцию
            da = fix(Q(t2)/2);

            % Получаем множитель
            sig = 2;
            is_in_circle_t = circ(sqrt((x2).^2+(y2).^2)/(L/2*2*da/l2));
            Pst = exp( -1*( (x2).^2+(y2).^2 )/(2*sig^2) )/sqrt(2*pi)/sig;
            Psm = ifft2(fft2(fftshift(Pst)).*fft2(is_in_circle_t));
            m_psm = max(max(Psm));
            Ps2 = Psm/m_psm;
            
            %% Получаем разбиение центральной части
            
            rg_o = 4;
            tau_o = 4;
            C = zeros( length(x2(is_in_circle2)).^2,1);            
            
            masl = false(l2,l2);
            [tau,rn] = cart2pol(x2,y2);
            
            rt = Rd2/rg_o; %радиус центральной части
            dtau = 2*pi/tau_o; %шаг по tau
            temp_circ = rn <= rt;
            masl(:,:,1) = temp_circ;
            Ntemp = 1;
            clear en;
            en(Ntemp+1) = size(x2(temp_circ),1);
            for ii2 = 1:rg_o-1
                temp_circ2 = (rn <= Rd2*(1/rg_o+ii2/rg_o)) & ~temp_circ;
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
            %%
            
            % Получаем exclude_new
            dat = da;
            da = fix(da/1);
            [x3,y3] = meshgrid(-1:1/(2*da):1-1/(2*da));
            x3 = x3 + 1/(4*da); y3 = y3 + 1/(4*da);
            exclude_new = circ(x3.^2+y3.^2);
            da = dat;
            
            mas_max = zeros(200,1); %Это максимум для прекращения работы
            mas_max(1) = 1;
            N = 1;
            sT = abs(E2);
            A = zeros(200,2);
            td = 50;
            P2 = false(size(tempz));
            P2(513-184+td:513+184-td,513-448+td:513+448-td) = true;
            cx_s = 513-184+td-1; cy_s = 513-448+td-1; %Для другого вида индексирования
            stop_while = true;
            G = zeros(l2,l2,200);
            sh = 0;
            da2 = l2/2;
            
           
            [Esort, Ecr] = sort(sT(P2),'descend');
            N = 1;
            N2 = 1;
            szx = 184*2-td*2+1;
            szy = 448*2-td*2+1;
            while stop_while && (199 > N)                
                
                if Esort(1)/Esort(N2) > 2.5
                    stop_while = false;
                end
                [cx,cy] = ind2sub([szx,szy],Ecr(N2));                
                cx = cx + cx_s;
                cy = cy + cy_s;
                if P2(cx,cy) == true
                    A(N,1) = cx;
                    A(N,1) = cy;                    
                    P2(cx-da*2:cx+da*2-1,cy-da*2:cy+da*2-1) =...
                    P2(cx-da*2:cx+da*2-1,cy-da*2:cy+da*2-1) & ~exclude_new;
                    S3 = E2(cx+sh-da2:cx+sh+da2-1,cy+sh-da2:cy+sh+da2-1).*Ps2;
                    G(:,:,N) = fftshift(fft2(fftshift(S3)));
                    N = N + 1;
                end
                N2 = N2+1;
            
            end
            
            sum_um = zeros(size(G(:,:,1)));
            sum_vm = zeros(size(G(:,:,1)));
            sizeC = size(C);
            fprintf('Количество точек = %g \n',N-1);
            for kk = 1:N-1
                Gt = G(:,:,kk);
                for ii = 1:Ntemp                            
                    xv = Gt(masl(:,:,ii)); 
                    Ct = xv*xv';
                    C(en(ii)+1:en(ii+1)) = C(en(ii)+1:en(ii+1)) + Ct(:);
                end
            end            
            
            
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
            Se = E2*0;
            aP = angle(P);
            [ masl_grad_x, masl_grad_y] =...
                                       align_pieces( masl, is_in_circle2 );
            a_end = grad_svd3( aP, l2, Z2, is_in_circle2, is_in_circle2,...
                                                masl_grad_x, masl_grad_y );
            a_end(1:3) = 0;
            Se(is_in_circle3) = Zmin*a_end;
            

            
        
        %%
                                  
            
            
        %%
        
        Gs(:,:,t2) =  Se; % Записываем текущее поле
        
        E2 = ifft2( u.*ifftshift( ...
                  exp( -1i*(sum(Gs(:,:,1:t2),3)) )  )); % суммируем его
        Sh0(t2) = shann_entropy(E2(tempz));
        
        fprintf('Энтропия внутри цикла = %g \n',Sh0(t2));
        
    end
    [~,smin] = min(Sh0);
    Se = sum(Gs(:,:,1:smin),3);
    E2 = ifft2( u.*ifftshift( exp(-1i*Se)) );         
end


if hNa == true    
    Q = [44,36,32,28,24,20,16,10];        
    N2 = size(Q,2);
    Gst = repmat(E*0,[1 1 N2+1]);
    Sh0(1) = shann_entropy(E(tempz));
    u = fft2(E);
    Ev = E;
    for t2 = 1:N2
        
        fprintf('цикл %g \n',t2);
        fprintf('окно = %g \n',Q(t2));

        [Se,~] = PGA_find4_fft_hNa( Ev,Q(t2),tempz,...
            radius,Global_constants,TypeWindow,radius,...
            Z, is_in_circle, Deg_hNa, E2);
        %Se = Se.*seamless_border;
        Gst(:,:,t2) =  Se;
        Ev = ifft2( u.*ifftshift(exp(-1i*sum(Gst(:,:,1:t2),3))) ) ;
        Sh0(t2+1) = shann_entropy(Ev(tempz));   

        fprintf('Энтропия внутри цикла = %g \n',Sh0(t2+1));
    end
    [~,smin] = min(Sh0);
    Gs = sum(Gst(:,:,1:smin),3);
    Gk = ifft2( u.*ifftshift( exp(-1i*Gs)) ); 
    
end
toc

ap4 = Z\Gs(is_in_circle);
