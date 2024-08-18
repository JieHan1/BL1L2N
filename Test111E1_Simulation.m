clear ;
clc;
close all;
%%

A= [1.8812, -5.1186, 1.0129, 1.0806, -9.5331;
    -2.2202, 3.8944, 1.0656, -1.0268, 8.4156;
    -1.9014, 1.1472, 0.8832, -1.0990, 2.4498;
    -1.0519, 2.5056, 3.9539, -0.3660, 7.1488;
    -0.9673, 3.0783, 3.9738, -0.4710, 8.3454;
    1.0234,0.9959, -3.1213, 0.5479, 0.4053;
    3.0021, 6.8872, -3.1319, 1.6138, 12.6750;
    4.8996, -1.1349, -1.9069, 2.4316, -2.9337;
    3.9053, 1.9739, -1.9989, 1.8808, 2.9146;
    3.9626, 3.0953, -2.0645, 1.9927, 4.8799];

[U,S,D] = svd(A);
S(5,5) = 0.0000137;
A = U*S*D';
clear U S D



fx = [1;1;1;1;1];
%  z_yn = [-10.5120;10.4430;1.4485;11.9400;14.0850;...
%      -0.1535;21.1920;1.6535; 8.9494; 11.8650];
A1=A;
z_y = A1*fx;
for ii= 1  :1000
    z_yn = z_y + normrnd(0, 0.2, size(z_y));%sigma =2,1,0.5,0.2

    x = (A'*A)\A'*z_yn;



    %%%
    [m,n]   = size(A);
    alpha1  = 1e-10;
    alpha   = LQuarve(A, z_yn);
    D = Construction_Toeplitz(x, [1,-1]);%eye(n);%

    A = A1*(D'*D)^(-1/2);

    [x_L11, ~, ~,alphaL] = balanced_L2_L1d2(A1, z_yn, alpha, fx);


    NNMSE_L1L1(ii)  = (norm(fx - x_L11));

end


function  lambda = LQuarve(A,zk)

lambda_set = linspace(0, 1, 1000);

x_ord = zeros(length(lambda_set),1);
y_ord = zeros(length(lambda_set),1);

for i = 1 : length(lambda_set)

    x_alpha =(A'*A+lambda_set(i)*eye(size(A,2)))\A'*zk;

    x_ord(i) = 2*log((norm(zk-A*x_alpha))^2);
    y_ord(i) = 2*log((norm(x_alpha)^2));
end

dist_xy = sqrt((x_ord-min(x_ord)).^2 + (y_ord-min(y_ord)).^2);
lambda= lambda_set(dist_xy == min(dist_xy));

end
