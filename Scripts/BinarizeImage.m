%Function for binarizing image in apptest 2
function [Im]= BinarizeImage(app,Im)
%%Load info 
FrangiBool = app.FrangiCheckBox.Value;  %check if Frangi is on 
Method = app.BinarizationMethodDropDown.Value; %Check selected binarization method 


%% %% Make sure vessels are white and the rest is black 
[y] = max(Im);
[~,cc] = max(y);
[~,rr] = max(Im(:,cc));

%% Binarize 
if FrangiBool
    opts.sigmarange = [1 app.FrangiSigma_MaxSlider.Value];
    opts.sigmastepsize=2;
    opts.correctionconst1=0.8;
    opts.correctionconst2=15;
    Im=frangi_2Dfilter(Im,opts);
end

if strcmp(Method,'Fuzzy Means')
    %Implements Fuzzy thrshold binarization 
    %Number of sets;
    nth=2;
    Im=fuzzy_thresholding(Im,nth,3)-1;  
    
elseif strcmp(Method,'Adaptive Thresholding')
    %Implements adaptive threshold binarization
    Im= double(adaptivethresholding(Im,70));
    %Maybe make  the number here edjustable? 
else
    warndlg('Incorrect Segmentation method') ;
end

Im(Im==Im(rr,cc)) = 255;
Im(Im~=Im(rr,cc)) = 0;




