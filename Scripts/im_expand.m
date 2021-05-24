function I_out=im_expand(I,Nx,Ny)
%
% Simetric expansion of image I
%
% (Assumes a DCT-like simetry)
%
% Y=IM_EXPAND(I,Nx,Ny)
%
% Asumes a specular simetry of image I and replicate Nx columms
% and Ny lines.
%
% If size(I)= [M1, M2]
%    then size (I_out)=[M1+Ny,M2+Nx]
%

Ir=flipdim(I,1);
I_N=Ir(end-(Nx-1):end,:);
I_S=Ir(1:Nx,:);
Ir=flipdim(I,2);
I_E=Ir(:,1:Ny);
I_W=Ir(:,end-(Ny-1):end);
Ir=flipdim(Ir,1);
I_NE=Ir(end-(Nx-1):end,1:Ny);
I_NO=Ir(end-(Nx-1):end,end-(Ny-1):end);
I_SE=Ir(1:Nx,1:Ny);
I_SO=Ir(1:Nx,end-(Ny-1):end);

I_out=[I_NO I_N I_NE
    I_W  I   I_E
    I_SO I_S I_SE];
