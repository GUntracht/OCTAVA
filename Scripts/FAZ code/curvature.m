function K = curvature(I)

[m,n]=size(I);

phix = [I(1:m,2)-I(1:m,1) 1/2*(I(1:m,3:n)-I(1:m,1:n-2)) I(1:m,n)-I(1:m,n-1)];
phiy = [I(2,1:n)-I(1,1:n);1/2*(I(3:m,1:n)-I(1:m-2,1:n));I(m,1:n)-I(m-1,1:n)];

phiyy = [I(3,1:n)-2*I(2,1:n)+I(1,1:n)
         I(1:m-2,1:n)-2*I(2:m-1,1:n)+I(3:m,1:n)
         I(m,1:n)-2*I(m-1,1:n)+I(m-2,1:n)];

phixx = [I(1:m,3)-2*I(1:m,2)+I(1:m,1) I(1:m,1:n-2)-2*I(1:m,2:n-1)+I(1:m,3:n) I(1:m,n)-2*I(1:m,n-1)+I(1:m,n-2)];

     
%phi_xy
xyu = (I(2,3:n)-I(1,3:n)-I(2,1:n-2)+I(1,1:n-2))/2; %upper
xyb = (I(m,3:n)-I(m-1,3:n)-I(m,1:n-2)+I(m-1,1:n-2))/2; %bottom
xyl = (I(3:m,2)-I(1:m-2,2)-I(3:m,1)+I(1:m-2,1))/2; %left
xyr = (I(3:m,n)-I(1:m-2,n)-I(3:m,n-1)+I(1:m-2,n-1))/2; %right
xyulc = I(2,2)-I(2,1)-I(1,2)+I(1,1); %upper left corner
xyurc = I(2,n)-I(1,n)-I(2,n-1)+I(1,n-1); %Upper right corner
xyllc = I(m,2)-I(m-1,2)-I(m,1)+I(m-1,1); %Lower left corner
xylrc = I(m,n)-I(m-1,n)-I(m,n-1)+I(m-1,n-1); %Lower rigth corner

%Putting them together
xyu = [xyulc xyu xyurc]; 
xyb = [xyllc xyb xylrc]; 


     
temp = [xyl (I(2:m-1,2:n-1)-I(1:m-2,3:n)-I(3:m,1:n-2)-I(1:m-2,1:n-2))/4 xyr];
phixy = [xyu ; temp ; xyb];

K = (phixx.*phiy.^2-2*phix.*(phiy.*phixy)+phiyy.*phix.^2)./(phix.^2+phiy.^2+eps).^(3/2);
K = K./max(max(abs(K))); % Normalize

end