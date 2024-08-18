function [x_bc, bias, biasT] = L1L2_biasCorrection(x, A, D, p, lambda, sigma, zk, omigak, fx, t)
alpha = p + lambda*(1-sigma);
Dx      =[ diff(x, 1, 1); x(1,:) - x( end,:)] ;
aa = (zk - omigak);
cc =  [ diff(aa, 1, 1);  aa(1,:) - aa( end,:)] ;


% x_bc = x + (A'*A+ alpha*(D'*D))\(alpha*D'*Dx-p*cc);%
% bias  = (A'*A+ alpha*(D'*D))\(alpha*D'*Dx);%-p*cc
% 
% Dfx1  = [fx;zeros(size(A,2)-length(fx),1)];
% Dfx    = [ diff(Dfx1, 1, 1); Dfx1(1,:) - Dfx1( end,:)] ;
% biasT =  (A'*A+ alpha*(D'*D))\(alpha*D'*Dfx);%-p*cc


nt = length(t);
[~,m] = size(A);
ND = ((D'*D)\(A'*A) + alpha*eye(m));
[~, S, U] = svd(ND);
STS  = S(1:nt,1:nt);
bias    =  (U(:,1:nt)*inv(STS)*U(:,1:nt)')*((D'*D)*(alpha*D'*Dx-p*cc));%
x_bc   = x +(U(:,1:nt)*inv(STS)*U(:,1:nt)')*((D'*D)*((alpha*D'*Dx-p*cc)));

Dfx1  = [fx;zeros(size(A,2)-length(fx),1)];
Dfx    = [ diff(Dfx1, 1, 1); Dfx1(1,:) - Dfx1(end,:)] ;
biasT =  (U(:,1:nt)*inv(STS)*U(:,1:nt)')*((D'*D)*(alpha*D'*Dfx));

end