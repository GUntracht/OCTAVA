function A = get_area(x,y)
%input: the x and y coordinates of the points in contour
%output: the area of the contour
%
  [a,n]=size(x);
  i=1;
  sum=0;
  while i<n
      sum=sum+ x(i)*y(i+1)-x(i+1)*y(i);
      i=i+1;
  end
    
  A=1/2*sum;
  
end