function Y=fuzzy_thresholding(X,nth,w)
% Reference:
%
%       Aja-Fernández, S., A. Hernán Curiale, and G. Vegas-Sánchez-Ferrero, A local fuzzy thresholding methodology for multiregion image segmentation", 
%       Knowledge-Based Systems, vol. 83, pp. 1-12, 07/2015.
%
X=double(X);

sigmaX=std(X(:));              
muX=mean(X(:));
norma=max(X((X>(muX-2.*sigmaX))&(X<(muX+2.*sigmaX))));
X=X./norma;
h=fspecial('gaussian', 11, 1.5);
gX= smooth_function(h,X);

% fuzzy c-means clustering
Xfcm=fcm(gX(:),nth,[2,100,1e-5,0]);    	
max_gX=max(gX(:));
min_gX=min(gX(:));
Xfcm=sort(Xfcm);

for imax=1:nth
    if imax==1 %First set
        a1=min_gX;
        a2=min_gX;
        a3=Xfcm(1);
        a4=Xfcm(2);
    elseif imax==nth %Last set
        a1=Xfcm(nth-1);
        a2=Xfcm(nth);
        a3=max_gX;
        a4=max_gX;
    else 
        a1=Xfcm(imax-1);
        a2=Xfcm(imax);
        a3=Xfcm(imax);
        a4=Xfcm(imax+1);
    end

    Xtrap=trapezoid_fuzzy_memberhip(X,a1,a2,a3,a4);
    Xtrap3(:,:,imax)=Xtrap; 
end 
   
[~, ~, n3]=size(Xtrap3);
for ii=1:n3
    medXtrap3(:,:,ii)=medfilt2(Xtrap3(:,:,ii),[w,w]);
end
[~, Y]=max(medXtrap3,[],3);
  
%Correct any possible side effect of median operator
Y=Y.*(Y>0)+(Y==0);

end

function Y=smooth_function(h,X)
    [Mx, My]=size(h);

    nx=(Mx-1)/2; ny=(My-1)/2;

    Xe=iexpansion(X,nx,ny);
    fXe=filter2(h,Xe);
    Y=fXe((nx+1):end-nx,(ny+1):end-ny);
end

function Y= iexpansion(X,nx,ny)
Xr=flipdim(X,1);
Xn=Xr(end-(nx-1):end,:);
Xs=Xr(1:nx,:);
Xr=flipdim(X,2);
Xe=Xr(:,1:ny);
Xw=Xr(:,end-(ny-1):end);
Xr=flipdim(Xr,1);
Xne=Xr(end-(nx-1):end,1:ny);
Xno=Xr(end-(nx-1):end,end-(ny-1):end);
Xse=Xr(1:nx,1:ny);
Xso=Xr(1:nx,end-(ny-1):end);

Y=[Xno Xn Xne
    Xw X Xe 
    Xso Xs Xse];
end

function y=trapezoid_fuzzy_memberhip(x,a,b,c,d)
x=double(x);
if (b<a)||(c<b)||(d<c)
	error('Bad params')
end

if (a==b)
	y=((x-c).*(-1)./(d-c)+1).*(x>c).*(x<d)+(x<=c);
elseif (b==c)
	y=(x-a)./(b-a).*(x<=b).*(x>a)+((x-c).*(-1)./(d-c)+1).*(x>c).*(x<d);
elseif (c==d)
	y=(x-a)./(b-a).*(x>a).*(x<b)+(x>=b);
else
	y=(x-a)./(b-a).*(x>a).*(x<b)+(x>=b).*(x<=c)+((x-c).*(-1)./(d-c)+1).*(x>c).*(x<d);
end

end
