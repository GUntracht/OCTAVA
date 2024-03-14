function [s,index] = get_lines(c)
%s.v:the value of the isoline
%s.x:the x coordinate(array) of a certain value
%s.y:the y coordinate(array) of a certain value
%index: the number of the maximum length contour

    sz = size(c,2);     % Size of the contour matrix c
    ii = 1;             % Index to keep track of current location
    jj = 1;             % Counter to keep track of # of contour lines
    maxi= 0;
    index= 1;
    while ii < sz       
        n = c(2,ii);    
        s(jj).v = c(1,ii);        % Value of the contour
        s(jj).x = c(1,ii+1:ii+n); % X coordinates
        s(jj).y = c(2,ii+1:ii+n); % Y coordinates        
         if(n>maxi)
            maxi=n;
            index=jj;
        end
        ii = ii + n + 1;          % Skip the next contour line
        jj = jj + 1;              % Increment number of contours
    end

end