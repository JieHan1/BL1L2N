function  lambda = MSEMin(A, x_alpha,  p, sigma, y)

lambda_max = 10;
lambda_min = eps;


D                 = Construction_Toeplitz(x_alpha, [1,-1]);
alpha_max = (p+ lambda_max*(1-sigma));
alpha_min = (p+ lambda_min*(1-sigma));
sigma0        = (y - A*x_alpha)'*(y - A*x_alpha)/abs(size(A,1)-size(A,2));

MSE1       = GetDMSE(A, D, x_alpha, sigma0, alpha_min);

if MSE1 > 0 
    lambda = 0;
    return;
end

% while 1 
%     MSE=0;
%     MSE = MSE+GetDMSE(A, D, x_alpha, sigma0, alpha_max);
%     if MSE<0
%         lambda_max   = 2*alpha_max;
%         alpha_max= (p+ lambda_max*(1-sigma));
%     else
%         break;
%     end
%     x_alpha = (A'*A + alpha_max*(D'*D))\(A'*y);
%     sigma0  = (y - A*x_alpha)'*(y - A*x_alpha)/abs(size(A,1)-size(A,2));
% 
% end

while abs(lambda_min - lambda_max) > 1e-7

        Lambda_Median = (lambda_max +lambda_min)/2;
        alpha_Median      = (p+Lambda_Median*(1-sigma));
        MSE       = GetDMSE(A, D, x_alpha, sigma0, alpha_Median);
        if MSE < 0
            lambda_min = Lambda_Median;
        else
            lambda_max = Lambda_Median;
        end

end

lambda = Lambda_Median;
end

function MSE = GetDMSE(A, D, x_alpha, sigma0, alpha )

Nk = A'*A + alpha*(D'*D);
Sr   = D'*D;

MSE1 = -sigma0* trace( (A'*A)*inv(Nk)*(inv(Nk)*Sr + Sr*inv(Nk) )*inv(Nk));
MSE2 = 2*alpha*x_alpha'*Sr*inv(Nk)*inv(Nk)*Sr*x_alpha;
MSE3 = alpha*alpha*x_alpha'*Sr*inv(Nk)*(inv(Nk)*Sr + Sr*inv(Nk) )*inv(Nk)*Sr*x_alpha;

MSE   = MSE1 + MSE2 - MSE3;
end