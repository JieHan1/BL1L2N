%%% reference:
%%%                  Lin Dongfang, et al., Bias Reduction Method for Parameter Inversion of
%%%%                Ill-Posed Surveying Model,Journal of Surveying
%%%%                Engineering,  2020

%%%           The new formulae have been derived and the function need to
%%%           be updated. The updated date is 20230328
function  [t] = Bias_Reduced_Han(A, D, x_alpha,  lambda, p, sigma, y)


sigma0 = (y - A*x_alpha)'*(y - A*x_alpha)/abs(size(A,1)-size(A,2));
alpha   = p+ lambda*(1-sigma);


Dx      =[ diff(x_alpha, 1, 1); x_alpha(1,:) - x_alpha( end,:)] ;

[~, S, G] = svd((D'*D)\A'*A);
[H, V, U] = svd(D);


sim = min(size(S,1), size(S,2));
t = [];


for i = 1 : sim

%     Dxx = Dx'*H(:,i)*V(i,i)*(V(i,i)'*V(i,i))^(-1)*U(:,i)'*G(:,i); %% 方式1
% 
%     Tb_Covi(i) = alpha*alpha*(G(:,i)'*U(:,i)*U(:,i)'*G(:,i))/...
%         ((S(i,i) + alpha)^2) * sigma0 /(S(i,i)*V(i,i)*V(i,i));%
%     Tb_xTxi(i)  = Dxx*alpha*alpha/(S(i,i) + alpha)^2*Dxx';


%%% 审稿意见  %%%%%%%%%%%%%%%%%%%
    Dxx = x_alpha'*G(:,i); 

%     Dxx = x_alpha'*G(:,i);
    Tb_Covi(i) = (alpha*alpha* sigma0 + 2*alpha*sigma0*S(i,i))*(G(:,i)'*...
        U(:,i)*U(:,i)'*G(:,i))/((S(i,i) + alpha)^2)  /(S(i,i)*V(i,i)*V(i,i));%
    Tb_xTxi(i)  = Dxx*alpha*alpha/(S(i,i) + alpha)^2*Dxx';


  %%%%%%%%%%%%%%%%%%%%%%%%%%  
    if Tb_xTxi(i) >Tb_Covi(i)
        t = [t; i];
    end

end

end 