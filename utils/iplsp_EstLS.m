%�¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢�
% Estimate the signal x using least square
%
% y		: observation
% A	: Sensing matrix
% Written by Seokbeop Kwon
% Information Processing Lab., Korea Univ.
%�¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢�

function x_ls = iplsp_EstLS(y, A)

% 	x_ls=A\y;
	x_ls=pinv(A)*y;
end