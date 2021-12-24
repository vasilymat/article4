function [ Gk, Gs ] = PGA_half_Na( E,Global_constants,TypeWindow,radius, Deg )
%������� ��������� �������� �����������, ���������� �� ����� � ���������
%�����
%subapertures - ���������� ����� �����������. 0 - ������ ���������� �������
%������
%radius - ������ ������ ��������
dx=Global_constants(4);
L=Global_constants(5);
[x,y]=meshgrid(-L/2:dx:L/2-dx);
x = x + dx/2;
y = y + dx/2;

Ev = E;
f = E(:,:,1)*0; %���� ����� ������������� ����
%Q = [40,36,32,24,20,16,8];

bs = 4;
N2 = bs+1;  

Q0 = 0:bs;
alfa = (40-16)/(bs^2);
Q = 40 - alfa*Q0.^2;
Q = round(Q/2)*2;
alfa2 = (0.05-0.05)/bs^2;
Q2 = 0.1+alfa2*Q0.^2;
Q = [40,36,32,24,20,16,8];

Gs = repmat(f,[1 1 N2+1]);
Gss = f*0;

[Z,a,~] = zerfit(x*0,x,y,radius,Deg);
Rex = 0;

tempz = false(1024,1024);%��� ����� �� ����� ��� �������� �������� �����������
tempz(513-184:513+183,513-448:513+447) = true;

is_in_circle = circ(sqrt(x.^2+y.^2)/radius);
Sh0(1) = shann_entropy(Ev(:,:,1));

disp(sprintf('�������� ��������� ����������� = %g',...
    shann_entropy(Ev(tempz)) ));


for t2 = 1:N2
        
    disp(sprintf('���� %g',t2));
    disp(sprintf('���� = %g',Q(t2)));
    
    [Se,a] = PGA_find4_fft_help( Ev,Q(t2),a,...
      Rex,Global_constants,TypeWindow,radius,Z, is_in_circle, Deg );    
    
    Ev = ifft2( fft2(Ev).*ifftshift((exp(-1i*Se))) );
    Sh0(t2+1) = shann_entropy(Ev(tempz));    
    
    disp(sprintf('�������� ������ ����� = %g',Sh0(t2+1)));
    
    % ���� ����� �������� ��������� ������� ����  
    Gs(:,:,t2+1) =  Se;
    
end
[~,smin] = min(Sh0);
Gst = sum(Gs(:,:,1:smin),3);
Gss = Gss + Gst;

Ev = ifft2( ifftshift( fftshift(fft2(E)).*exp(-1i*Gss) ) );

Gs = Gss;
Gk = Ev;

end

