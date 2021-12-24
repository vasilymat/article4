function f=rphase3(dx,L,k)
%������� ���������� ������ ��� �������� ���������
%  
%��� ��� ������� ��������
[kx,ky]=meshgrid(-1/(2*dx):1/L:1/(2*dx)-1/L);
KZ=2*pi*sqrt((k/(2*pi))^2-kx.^2-ky.^2); 
%�����

% %��� ��� ������ �� �������
% [x,y]=meshgrid(-L/2:dx:L/2-dx);
% KZ=x.^2+y.^2;
% %�����

f=KZ;
end