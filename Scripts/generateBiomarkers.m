if Batch
    dir = app.batchSaveFile;
    app.biomarkers.name = app.FileLabelName;

else
    dir = app.singleSaveFile;
end

app.biomarkers.name = app.FileLabelName;
app.csvFile = char(app.csvFile);
if ~(app.csvFile(end) == 'v')
    app.csvFile = app.csvFile(1:end-1);
end
    
inds= strfind(app.csvFile,'Table');
inds2= strfind(app.csvFile,'.');
imgnum = app.csvFile(inds+6:inds2-1);
if ~isempty(imgnum)
    imgnum = ['-' imgnum ];
end


[Nx, ~] = size(app.BinaryImage);
pix_to_mm = app.RealImageSizemm.Value/Nx; %This assumes images are square - we may want to update this! 



%LoadDiameterTable - Loads the diameter, length, and tortuosity
%of each identified vessel in the image 
DiamTableFile = [dir '\CSVs\' app.csvFile];
DiamTable = xlsread(DiamTableFile); 
app.biomarkers.Diams = rmmissing(DiamTable(:,3)*pix_to_mm*1000); %Diameters (in microns)
app.biomarkers.Tort = rmmissing(DiamTable(:,21)); %Tortuosity
app.biomarkers.Lengths = rmmissing(DiamTable(:,19))*pix_to_mm; %In mm

%Read ROIS to count nodes
%Still not 100% sure this is correct - need to double check
%junctions vs nodes 
% ROI_folder = [dir '\ROIs\' app.FileLabelName '_Network-RoiSet' imgnum '.zip'];
% f = unzip(ROI_folder);
%     app.biomarkers.nodes = 0;
% for i = 1:length(f)
%     tmp = [strfind(f{i},'unction') strfind(f{i},'node')];
%     if ~isempty(tmp)
%         app.biomarkers.nodes = app.biomarkers.nodes+1;
%     end
%     delete(f{i})
% end

%% Calulate Biomarkers

%Vessel Area Fraction 
app.biomarkers.VAD = sum(sum(app.BinaryImage))/255/length(app.BinaryImage(:));

%Vessel Length Density
app.biomarkers.totalLength = sum(app.biomarkers.Lengths); %In mm
app.biomarkers.VLD = app.biomarkers.totalLength/app.RealImageSizemm.Value^2; %in 1/mm

%Mean Diameter
app.biomarkers.MeanDiam = mean(app.biomarkers.Diams); %In microns
app.biomarkers.MedianDiam = median(app.biomarkers.Diams); 

%Branchpoint density
app.biomarkers.BD =  app.biomarkers.nodes/app.biomarkers.totalLength; %in 1/mm

%Display results - caluclulated biomarkers
app.VesselAreaDensityLabel.Text = ['Vessel Area Density = ' num2str(app.biomarkers.VAD)]; 
app.VesselLengthDensityLabel.Text = ['Vessel Length Density = ' num2str(app.biomarkers.VLD) ' (1/mm)'];
app.MeanDiameterLabel.Text = ['Mean Diameter = ' num2str(app.biomarkers.MeanDiam) ' (' char(181) 'm)'];
app.MedianDiameterLabel.Text = ['Median Diameter = ' num2str(app.biomarkers.MedianDiam) ' (' char(181) 'm)'];
app.BranchpointDensityLabel.Text = ['Branchpoint Density = ' num2str(app.biomarkers.BD) ' (1/mm)'];

%Generate histograms of diameter, tortuosity 
%cens_diams = linspace(0,300,31);  
%dHist = hist(Diams,cens_diams);
%bar(cens_diams,dHist,'Parent',app.DiamsAxes)
%Will using histogram vs hist cause an issue with non-uniform
%bin size between datasets? Need to check this! 
if ~Batch
    histogram(app.biomarkers.Diams,'BinWidth',10,'Parent',app.DiamsAxes);%Make Bin size adjustable?
    histogram(app.biomarkers.Tort,'BinWidth',0.1,'Parent', app.TortsAxes);
    histogram(app.biomarkers.Lengths,'BinWidth',0.1,'Parent', app.LengthsAxes);
end

app.biomarkers_counter = app.biomarkers_counter+1;
app.biomarkers_cumulative(app.biomarkers_counter).name = app.biomarkers.name;
app.biomarkers_cumulative(app.biomarkers_counter).nodes = app.biomarkers.nodes;
app.biomarkers_cumulative(app.biomarkers_counter).VAD = app.biomarkers.VAD;
app.biomarkers_cumulative(app.biomarkers_counter).totalLength = app.biomarkers.totalLength;
app.biomarkers_cumulative(app.biomarkers_counter).VLD = app.biomarkers.VLD;
app.biomarkers_cumulative(app.biomarkers_counter).MeanDiam = app.biomarkers.MeanDiam;
app.biomarkers_cumulative(app.biomarkers_counter).MedianDiam = app.biomarkers.MedianDiam;
app.biomarkers_cumulative(app.biomarkers_counter).BD = app.biomarkers.BD;


