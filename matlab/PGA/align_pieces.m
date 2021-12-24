function [ masl_grad_x, masl_grad_y ] = align_pieces( masl, is_in_circle )
%Функция выравнивает константы между кусками апертурыр

N = size(masl,3);
l2 = size(masl,2);
S = zeros(l2,l2);

for ii = 1:N
    S(masl(:,:,ii)) = ii;
end


dS = diff(S,1,1);
is_in_circle2 = is_in_circle(2:end,:);
masl_grad_x = dS ~= 0 & is_in_circle2 ;


dS = diff(S,1,2);
is_in_circle2 = is_in_circle(:,2:end);
masl_grad_y = dS ~= 0 & is_in_circle2 ;

% dS = diff(S,1,1);
% is_in_circle2 = is_in_circle(2:end,:);
% masl_grad11 = dS == -1 & is_in_circle2 ;
% masl_grad12 = dS == 1 & is_in_circle2 ;
% 
% dS = diff(S,1,2);
% is_in_circle2 = is_in_circle(:,2:end);
% masl_grad21 = dS == -1 & is_in_circle2 ;
% masl_grad22 = dS == 1 & is_in_circle2 ;
% 
% masl_grad_x = masl_grad11 | masl_grad12;
% masl_grad_y = masl_grad21 | masl_grad22;
    

end

