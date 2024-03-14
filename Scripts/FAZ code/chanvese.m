function [phi,s] = chanvese(I,n,my,lambda1,lambda2,maxsize)
% Input
% I is the Image 
% n is the number of iterations
% my, lambda1 and lambda2 are parameters for Chan-Vese
% maxsize is the max number of pixels on the long edge on the output image


% Things for the stop criteria
S = stop({'Stop me on the ok button:'}) ;


% Rescale image/pixel values
temp = double(I)./255;


% Rescaling the size of the image
scale = maxsize/min(size(temp));
   
if scale<=1   
    temp = imresize(temp,scale);
    I = imresize(I,scale);
end



% Choose initial points
figure, imshow(temp,'Initial','fit'); 
title('Select initial points, end by pressing "enter" on the key board'); hold on; 
[x1,y1] = ginput;
hold off;
  
% the levelset function
[x,y] = meshgrid(1:size(temp,2),1:size(temp,1));
m = inpolygon(x,y,x1,y1);
phi = bwdist(1-m)-bwdist(m)-m;

% get the initial contour
imshow(temp); hold on;
title('Initial contour');
contour(phi, [0 0], 'r', 'LineWidth',2); drawnow; hold off;
figure,


% Iterativ method 
dt = 1;

%Chan-Vese method
for i=1:n
    %get the c- and c+
    heavi = heaviside(phi);
    cp = sum(sum(heavi.*temp))/length(find(phi>=0)); %inside mean gray level
    cm = sum(sum((1-heavi).*temp))/length(find(phi<0)); %outside mean gray level
    % Curvature
    kappa = curvature(phi);

    % Iterative step
    phi_t = diracdelta(phi).*(my*kappa - lambda1*(temp-cp).^2 + lambda2*(temp-cm).^2);
    phi = phi + dt.*(phi_t./(max(max(abs(phi_t))))); 

    % Print contour
    if(mod(i-1,10) == 0) 
        imshow(temp); hold on;
        contour(phi, [0 0], 'r', 'LineWidth',2); drawnow; hold off;
        
        xlabel('The iteration can be stopped by pressing the OK button')
    end
    
    % Stop iteration
    if S.Stop();
        % Print Orginal image with the final contour
        figure,
        imshow(256-I); hold on;
        title(['Finally contour on original image after ' num2str(i) ' iterations']);
        [c,h]=contour(phi, [0 0], 'k', 'LineWidth',0.02); drawnow; %hold off;
        %To remove the small circles
        [s,index]=get_lines(c);
        %s(index).x,s(index).y are the co
        plot(s(index).x,s(index).y,'r','LineWidth',2);
        %get the center of the avascular zone
        [Cx,Cy]=get_central(s(index).x,s(index).y);
        plot(Cx,Cy,'r*');
%         % Plot level set
         figure,
         mesh(double(phi)); hold on;
         mesh(double(0*phi)); hold off;
        return;
    end

    % Reinitialize the level set equation phi
    phi = reinitialization(phi, 0.5);
    
end
end

function Y = heaviside(X)
epsilon = 0.3;
Y=1/2*(1+2/pi*atan(X/epsilon));
end

function d = diracdelta(X)
epsilon = 0.3;
d = (pi*(epsilon^2+X.^2)).\epsilon;
end
