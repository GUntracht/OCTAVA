%% This funtion allows the user to load a single OCT MIP into the app

[Filename, Pathname] = uigetfile({'*.*'}, 'Select MIP');

if Filename == 0
    warndlg('You did not choose an image file') ; 
else
    a = fullfile(Pathname, Filename);
    a = imread(a);
    app.inputImage = double(int16(squeeze(a(:,:,1))));


    imagesc( app.inputImage, 'Parent', app.inputImageAxes);
end

inds = strfind(Filename,'.');
app.FileLabelName = Filename(1:inds(end)-1);
app.FileLabel.Text = app.FileLabelName;



