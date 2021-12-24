addpath('C:\MatkivskiyV\!Актуальная наука\Статья 4\matlab\PGA_release_2\HelpFun\','-end')
addpath('C:\MatkivskiyV\!Актуальная наука\Статья 4\matlab\PGA_release_2\helpPGA\','-end')

NA_size = 2; %(1 - небольшая апертура, 2 - большая апертура)

zdim = 1;
readdata = true;
clear N_layer;
switch NA_size     
    case 1
        dik1 = fopen('c:\MatkivskiyV\!Актуальная наука\Статья 4\matlab\data\Retina_LowNA_firstVol_896x368x128_complex.raw');
        N_layer = 75:74+zdim;%10:9+zdim; %хороший слой, но можно взять любой другой /75 слой
        %N_layer = 36:35+zdim;
        
    case 2
        dik1=fopen('c:\MatkivskiyV\!Актуальная наука\Статья 4\matlab\data\Retina_HighNA_firstVol_896x368x128_complex.raw');
        N_layer = 65:64+zdim; %65й слой
        %nnn = 73;
        %N_layer = nnn:nnn-1+zdim; %65й слой
end


if readdata == true
D = single(zeros(896,368,128));
    for i=1:128      
        A = fread(dik1,[896*2,368],'single');
        ReA = A(1:2:896*2,1:368);
        ImA = A(2:2:896*2,1:368);
        B = ReA+1i*ImA;
        D(:,:,i) = B(:,:);
        
        %figure(1); imagesc((ReA.^2+ImA.^2).^0.5); colorbar ; colormap(gray);
             
          
        
    end
fclose(dik1);
end
    

%задается функция circ на основе реальных физ. параметров
l = 1024;
lam=0.000841;
k=2*pi/(lam);
dx=20*10^-3;
L=dx*l; 
[x,y]=meshgrid(-L/2:dx:L/2-dx);
kz=rphase3(dx,L,k); % This is phase for angular spectrum method
radius = 9.32; % Это радиус для любек-файлов
%radius = L/2;
%radius = 9.32;
circf = abs(sqrt(x.^2+y.^2)/radius) <= 1;

% circf = imgaussfilt(circ( sqrt(x.^2+y.^2)...
% /(radius - 7*dx*2)),7);

a=896; b=368; N = 128;% размеры исходного изображения

var = 2;
switch var
    case 1 %Изображение вставляется в 1024x512
        E = single(zeros(1024,512));
        E(28:995,72:439) = B(1:896,1:368);
    case 2 %Изображение вставляется в квадратный кадр NxN
        E = zeros(l,l,zdim);        
        E(l/2+1-a/2:l/2+a/2,l/2+1-b/2:l/2+b/2,1:zdim) = D(1:a,1:b,N_layer);
end
%figure(1); imagesc(abs(E)); colorbar ; colormap(gray);
%fE = fftshift(fft2(fftshift(E)));

%clear D; % Это удаляем когда намерены экономить память

switch NA_size
    case 1
        % %Двигаем поле в центр и режем апертурой. Эта история работала для файла 
        % %Retina_LowNA_firstVol_896x368x128_complex.raw
        for t1 = 1:zdim
            fE = fftshift(fft2(fftshift(E(:,:,t1))));
            fE = PL(-72.42-17,kz,fE);            
            fE = circshift(fE,[0 20]).*circf;
            E(:,:,t1) = double(ifftshift(ifft2(ifftshift(fE))));
            %E(:,:,t1) = rot90(E(:,:,t1)/3000/300,2);
        end
        
    case 2

        %Двигаем поле в центр и режем апертурой. Эта история работала для файла 
        %Retina_HighNA_firstVol_896x368x128_complex.raw
%         circf = LPF( circf, 20 );
        for t1 = 1:zdim
            fE = fftshift(fft2(fftshift(E(:,:,t1))));
            fE = PL(-153-17+170-140,kz,fE);
            fE = circshift(fE,[70 -15]).*circf;
            
            %figure(120); imagesc(abs(fE)); colorbar ; colormap(gray); caxis([0 10^5]);

            %fE = circshift(fE,[70 -10]).*circf;
            %fE = circshift(fE,[70 -10]).*circf;
            
            %fE = circshift(fE,[55 26]).*circf;            
            %fE = circshift(fE,[40 -10]).*circf; %[23 -46]
            E(:,:,t1) = double(ifftshift(ifft2(ifftshift(fE))));
            %E(:,:,t1) = rot90(E(:,:,t1)/3000/300,0);
        end
        
        %fE = PL(-80,kz,fE);
        %fE = circshift(fE,[50 -25]).*circf;                
end

%figure(11); imagesc(abs(fE)); colorbar ; colormap(gray);


tZ = zeros(l,l);
tZ(l/2+1-b/2:l/2+b/2,l/2+1-a/2:l/2+a/2) = 1;
for t1 = 1:zdim    
    E(:,:,t1) = flipud(rot90(E(:,:,t1))); % представляем данные в виде как в статье
    E(:,:,t1) = E(:,:,t1).*tZ;
end
%figure(23); imagesc(abs(fftshift(fft2(fftshift(E(:,:,1)))))); colorbar ; colormap(gray); caxis([0 15*10^5])


clear fE A B ReA ImA