function [Cx,Cy] = get_central(x,y)
%input: the x coordinates and the y coordinates of the contour
%output: the centoid's coordinates of 


  [a,n]=size(x);
  %To make the contour close
   x(n+1)=x(1);y(n+1)=y(1);
  A=get_area(x,y);%the area of the contour
 
  sum_x=0;
  sum_y=0;
  i=1;
  while i<n+1
      sum_x=sum_x+ (x(i)*y(i+1)-x(i+1)*y(i))*(x(i)+x(i+1));
      sum_y=sum_y+ (x(i)*y(i+1)-x(i+1)*y(i))*(y(i)+y(i+1));
      i=i+1;
  end
    
  Cx=sum_x/(6*A);
  Cy=sum_y/(6*A);
  
end