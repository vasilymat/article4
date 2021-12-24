function [ a, Se, E] =...
                  coefficient_finding_defocus(u,temp,Z,a,is_in_circle)
%Тут еще можно серьезно дооптимизировать предварительно вычислив полиномы
%Ц. и уменьшив fft2 размер до 256х256

%Функция для поиска аберраций методом подбора коэффициентов
%Ev - субобъем
%subapertures - количество шагов ФНЧ. 1 - значит нет ФНЧ
%radius - радиус полной апертуры
%Deg - максимальная степень полиномов
   
% options = optimoptions('fminunc','Algorithm','quasi-newton','GradObj','on');
options.MaxFunEvals = 50;


[xdim, ydim, zdim] = size(u);
Se = zeros(xdim,ydim);

a(5) = fminunc(@sh_entr_fft, a(5),options );                

E = ifft2(  u.*fftshift(exp(-1i*Se))  ) ;

Sh0 = shann_entropy(E(temp));
fprintf('Энтропия Ev_redused = %g \n',Sh0);

    function [Sh,dS] = sh_entr_fft(a_local)
    %Подбираем коэффициенты, ориентируясь на эффект в нескольких плоскостях.
        a(5) = a_local;        
        Se(is_in_circle) = Z*a;        
        U = ifft2( u.*repmat( fftshift(exp(-1i*Se)), [1 1 zdim]) );
        
        I = U.*conj(U);     
        nI = sum(sum(I));
        I = I./nI;
        lI = log(I);
        lI(~temp) = 0;
        
        DFT_prod = ifft2( (1+lI).*conj(U/sqrt(nI)) );
        Im_prod = ifftshift( imag( DFT_prod.*u/sqrt(nI)...
                                             .*fftshift(exp(-1i*Se)) ) );
        dS = -2*sum( Im_prod(is_in_circle).*Z(:,5) );
        Sh = shann_entropy( U(temp) );
    end
end

