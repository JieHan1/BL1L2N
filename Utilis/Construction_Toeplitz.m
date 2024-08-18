function [H] = Construction_Toeplitz(OneD_GT_Signal, psf_Signal)
% psf_Toeplitz_size --two elements matrix represents the size of Toeplitz
%                     first elements is row number
%                     second elements is column number
% psf_Signal        --column vector
psf_Toeplitz_size  = [length(OneD_GT_Signal), length(OneD_GT_Signal)];
L2 = length(psf_Signal);
H_row = zeros(psf_Toeplitz_size(2)+L2,1);%ÐÐ
H_col = zeros(psf_Toeplitz_size(1),1);%ÁÐ
H_col(1) = psf_Signal(1,1);
if size(psf_Signal,1) ==1
    H_row(1:size(psf_Signal,2)) = psf_Signal';
else
    H_row(1:size(psf_Signal,1)) = psf_Signal;
end

H1= toeplitz(H_row,H_col);

L = ceil((L2-1)/2);
H = H1(L+1:psf_Toeplitz_size(2)+L,:);
end