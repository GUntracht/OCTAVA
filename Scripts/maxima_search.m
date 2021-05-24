function Nmax=maxima_search(I,Th)

% Search of the maxima in the histogram of I
% 
% Nmax=maxima_search(I,Th)
%
% INPUT:
%           I:      input image
%           Th:     pruning threshold



dim = ndims(I);

if dim==2  
    window = fspecial('gaussian', 11, 1.5);
    I2=filter2B(window,I);
elseif dim>2
    SM=1;
    if SM==1
        % Gaussian smooth
        siz=[5 5 5];% 11 in each dimension
        sig = smooth.*[1 1 1]; % 1.5 in each dimension
        siz   = (siz-1)/2;
        [x,y,z] = ndgrid(-siz(1):siz(1),-siz(2):siz(2),-siz(3):siz(3));
        h = exp(-(x.*x/2/sig(1)^2 + y.*y/2/sig(2)^2 + z.*z/2/sig(3)^2));
        window = h/sum(h(:));    
        I2 = convn(I,window,'valid');
     else
         window = fspecial('gaussian', 11, 1.5);
         for ii=1:Mz
            I2(:,:,ii)=filter2B(window,I(:,:,ii));
         end
     end
end

S1=std(I2(:));              
M1=mean(I2(:));
Mw=max(I2((I2>(M1-2.*S1))&(I2<(M1+2.*S1))));

I2=I2./Mw;

[h1,x]=hist(I2(:),100);
h2=conv(ones([7,1])/7,h1);
h=h2(4:(end-3));
h=(h./sum(h))./(x(2)-x(1));

Dr=diff(h);

%MAX search----------
%Find the points when the function (dDr) is 0.
M=[Dr(1:end-1).*Dr(2:end)];
Maxt=find(double(M<0).*(Dr(1:end-1)>0))+1;
Nmax=length(Maxt); %(Number of fuzzy sets)		


%Prune of maxima with height smaller than Th% of maximum height

Tp=h(Maxt);
Tp=Tp./max(Tp);
Ns2=sum(Tp>Th);
if Ns2<(Nmax-1)
    Nth2=Ns2+1;
    %Maxt2=Maxt(Tp>Th);
    %Max2=floor(mean([Maxt(Nth2),length(x)]));
    %Maxt=[Maxt2 Max2];
    Nmax=Nth2;   
end
