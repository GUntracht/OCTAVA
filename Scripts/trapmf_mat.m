function y=trapmf_mat(x,a,b,c,d)

% Construction of trapezoid-shaped fuzzy membership functions
% (Matrix Version)
%
%y=trapmf_mat(x,a,b,c,d)
%
% INPUT:
%       x:          input value (number of matrix)
%       [a,b,c,d]:  main points of the trapezoid
%           If b=c: triangle
%           If a=b: initial trapezoid inicial (a and b <x)
%           If c=d: final trapezoid (c and d >x)
%
% OUTPUT: 
%       y: evaluation of x in the membership function defined by [a,b,c,d]

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

