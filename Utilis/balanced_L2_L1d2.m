function [x_final, MSE, MSET, lambda,tn] = balanced_L2_L1d2(A1, y, lambda1,fx, D )

%%% objective function:
%%%                                phi(x) = min ||Ax-y||_2 + alpha*||x||_1
%%%                                + beta*||x||_2
%%%  Reference
%%%           Weibo Huo et al.,2022, balanced Tikhonov and Total Variation
%%%           Decovolution approach for Radar Forward-looking
%%%           Super-Resolution Imaging. GRSL,19:3505805.

% addpath('Addin-Toolbox\l1magic\Optimization\');
[n,m1] = size(A1);
A2     = zeros(n,min(n,m1));
A       = [A1, A2];
[~,m] = size(A);

sigma    = 0.5;             %% 过程中通过必要公式进行更新
x_alpha = (A'*A + lambda1*eye(m))\(A'*y);

%%% 梯度运算的矩阵形式
D = Construction_Toeplitz(x_alpha, [1,-1]);
% D           = eye(m);%%  Satellite imagery Positioning select the regularization matrix
p1          = 0.001;%%0.001
omiga   = zeros(m,1);    %% 过程中通过必要公式进行更新
lambda = lambda1;


sl = 10;
iter = 1;
while 1 

    p = sl*p1;
    inner_iter = 1;
    while p < 1e1
        
            zk = soft_thresholding(D, x_alpha, omiga, p, lambda,  sigma);

            xk = EstimateX(A, p, D, lambda, sigma, y, zk, omiga);
            [t] = Bias_Reduced_Han(A, D, x_alpha,  lambda, p, sigma, y);
            [xk,bias, biasT] = L1L2_biasCorrection(xk, A, D, p, lambda, sigma, zk, omiga,fx, t);
            
%             xk(xk<0) = 0;
            cc =[ diff(xk, 1, 1); xk(1,:) - xk( end,:)] ;% D*xk;%
            omigak = omiga + p*(cc - zk);

            sigmak  = EstimateSigma(D, xk);
%             sigmak = 0;
            %%% L-曲线更新lambda
            if inner_iter == 1
                lambda = MSEMin(A1, x_alpha(1:m1),  p, sigma, y); 
                if lambda == 0
                    lambda = LQuarve(A, p, D,  sigmak, y, zk, omigak);
                end
            end

            if norm(x_alpha - xk)/norm(xk) < 1e-5||  inner_iter >100
                x_alpha = xk;
                break;
            else
                p = sl*p;
                sigma    = sigmak;
                omiga    = omigak;
                x_alpha = xk;
                inner_iter = inner_iter + 1;
            end
    

    end
    iter = iter + 1;
    if iter >10

%         [x_alpha,bias, biasT] = L1L2_biasCorrection(x_alpha, A, D, p, lambda,...
%             zeros(size(sigma)), zeros(size(zk)), zeros(size(omiga)),fx, t);
        x_final  = x_alpha(1:m1);

        
        sigma0 = (y - A1*x_final)'*(y - A1*x_final)/abs(size(A1,1)-size(A1,2));
        D1         = Construction_Toeplitz(x_final, [1,-1]);
        ND       = A1'*A1 + (p+ lambda*(1-sigma))*(D1'*D1);
        MSE     = trace(sigma0*((ND)\(A1'*A1))/(ND)) + (bias(1:m1)'*bias(1:m1));
        MSET  = trace(sigma0*((ND)\(A1'*A1))/(ND)) + (biasT(1:m1)'*biasT(1:m1));

        
        break;
    end


end
tn = length(t);

end

function zk = soft_thresholding(D, x_alpha, omiga, p, lambda,  sigma)
% cc = [ diff(x_alpha, 1, 2), x_alpha(:, 1) - x_alpha(:, end)] ;
cc = D*x_alpha;
zk = omiga + cc;
zk(omiga + cc > lambda*sigma/p)  =zk(omiga + cc > lambda*sigma/p)- lambda*sigma/p;
zk(omiga + cc < -lambda*sigma/p) =zk(omiga + cc < -lambda*sigma/p)+ lambda*sigma/p;
zk(omiga + cc > lambda*sigma/p & omiga + cc < -lambda*sigma/p) = 0;
    
end

function xk = EstimateX(A, p, D, lambda, sigma, y, zk, omiga)

aa = (zk - omiga);
cc =  [ diff(aa, 1, 1);  aa(1,:) - aa( end,:)] ;
xk_up = A'*y + p*cc;
xk_dw = A'*A + (p+ lambda*(1-sigma))*(D'*D) ;
xk        = xk_dw\xk_up;

end

function  lambda = LQuarve(A, p, D, sigma, y, zk, omiga)

lambda_set = linspace(0, 1, 1000);

x_ord = zeros(length(lambda_set),1);
y_ord = zeros(length(lambda_set),1);

for i = 1 : length(lambda_set)

        x_alpha = EstimateX(A, p, D, lambda_set(i), sigma, y, zk, omiga);
        
%         cc = [ diff(x_alpha, 1, 1);  x_alpha(1,:) - x_alpha( end,:)] ;
        cc = D*x_alpha;
        x_ord(i) = 2*log((norm(y-A*x_alpha))^2);
        y_ord(i) = 2*log((1-sigma)*(norm(x_alpha)^2) + sigma*norm(cc,1));
end

dist_xy = sqrt((x_ord-min(x_ord)).^2 + (y_ord-min(y_ord)).^2);
lambda= lambda_set(dist_xy == min(dist_xy));

end


function sigmak  = EstimateSigma(D, xk)

step   = (abs(D*xk)).^(1/3);
step1 = step(:);
step2 = (step1 - min(step1))/(max(step1)-min(step1));
sigmak = sum(step2)/length(step2);

end