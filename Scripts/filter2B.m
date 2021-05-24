function I_out=filter2B(h,I)

% Two-dimensional digital filter, based on FILTER2.m
%
%   Y = FILTER2B(B,X) filters the data in X with the 2-D FIR
%   filter in the matrix B.  The result, Y, is computed 
%   using 2-D correlation and is the same size as X. 
%
% The image X is periodically repeated to avoid edges problems.
%
%
%   FILTER2B uses FILTER2 to do most of the work.  


[Mx, My]=size(h);
if (rem(Mx,2)==0)||(rem(My,2)==0)
        error('h size must be odd');
end

Nx=(Mx-1)/2;
Ny=(My-1)/2;

% Expand the image
It=im_expand(I,Nx,Ny);
%Filter
I2=filter2(h,It);
%Original size
I_out=I2((Nx+1):end-Nx,(Ny+1):end-Ny);
