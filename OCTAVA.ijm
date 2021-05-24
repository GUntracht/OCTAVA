/***************************************************************************************
****************************************************************************************
GENERAL INFO
Developer: Rolando Szilveszter Matos and Gavrielle Untracht
Contact: r.matos@surrey.ac.uk
Campagnolo Lab and Biophotonics Lab:
	University of Surrey - FHMS - Bioscience Department
	Lead Investigators:
		Campagnolo, Paola Dr (Sch of Biosci & Med) <p.campagnolo@surrey.ac.uk>
		Heiss, Christian Prof (Sch of Biosci & Med) <c.heiss@surrey.ac.uk>
		McVey, John Prof (Sch of Biosci & Med) <j.mcvey@surrey.ac.uk>
		Sampson, Danuta Dr (Elec Electronic Eng) <danuta.sampson@surrey.ac.uk>

For More information consult README.MD file
This version has been modified for use with the MATLAB MOVA Software
***************************************************************************************
/***************************************************************************************
//---------------------        Version Information      ----------------------
//Last Updated: 12 March 2021
Version: MatLab_v3
State: Under Development
**************************************************************************************/
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//---------------------        Baseline Variables      ----------------------
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//<sectionStart

//****************************Function Universal Variables


run("Colors...", "foreground=white background=black selection=yellow");
run("Options...", "iterations=1 count=1 black");


//Angiogenesis Analyser Table name
//Changed all occurrences of "Stat Results Table" to angStatsTable to make the name change easier
var angStatsTable = "Stat Results Table";
var statsHeading = "\\Headings:Image Name\tAnalysed area\tNb extrem.\tNb nodes\tNb Junctions\tNb master junction\tNb master segments\tTot. master segments lenght\tNb meshes\tTot.meshes area\tNb peaces\tNb segments\tNb branches\tNb isol. seg.\tTot. lenght\tTot. branching lenght\tTot. segments lenght\tTot. branches lenght\tTot. isol. branches lenght\tBranching interval\tMesh index\tMean Mesh Size\tPath\tTotal Area\tMask Area Fraction\tSegment Area Fraction\tMean Diameter\tMean Len\tMean Endpoint Distance\tMean Tortuosity \tID-Mask Diamater\tMask Diameter_Mean\tMask Diameter_StDev\tMask Diameter_Mode\tMask Diameter_Min\tMask Diameter_Max\tMask Diameter_Perim\tMask Diameter_IntDen\tMask Diameter_Median\tMask Diameter_RawIntDen\tUse Auto Threshold\tAuto Threshold Value\tTwig Size\tUse Spackle\tSpackle Radious\tPixel Width\tUnit\tOriginal Width\tNew Width";

//My Tables
var myStatsTable = "Summary of Results";
var statsLine = "";
var myRestuls = "Complete list of Results";
var myResultsHead = "\\Headings:#\tLabel\tAverage Thickness\tStdDev\tMode\tMin\tMax\tPerim.\tBX\tBY\tWidth\tHeight\tMedian\tMinThr\tMaxThr\tImage ID\tComponent (code)\tComponent\tAdjusted Len\tEndpoint Distance\tTortuosity";

//All the type of values measured
var componentKeys = newArray(
			"Branch",
			"Extremity",
			"Extremity-Edge",
			"Isolated-Element",
			"Isolated-Twig",
			"Junction",
			"Master-Junction",
			"Master-Segment",
			"Mesh",
			"Node",
			"Segment",
			"Twigs"
		);

//Other variables
var isCMP = 0;

//sectionEnd>

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//----------        Angiogenesis Analyser U-Variables ORG      --------------
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//<sectionStart
//Was copied from --> Angiogenesis Original Variable Poz
var mapJunctionID=0, mapSegmentID=0, mapBranchID=0, mapMasterJunctionID=0;
var initid=0,workid=0,workingima="",mapExtremaID=0,mapNodesID=0,starttime,q3id,excludeObjectsSize=10,seconds=0;
var twigSize=25,nbExtrema=0,nbNodes=0,nbJunctions=0,nbLoop=0,nbTwigs,nbSegments,nbIsolated,nbTotalPeaces,loopSize=850, seuilIsolFinal=25,showLoopInBlue=0;
var nbMesh=0, totalMeshSize=0, initSkeleton=0;
var totalLenght=0, totalSegmentsLenght=0, totalTwigsLenght=0, totalBranchingLenght=0,totalIsolatedLenght=0,analyzedArea=0,branchingIndex=0,meshIndex=0,meanMeshSize=0;
var findloopOption=1,pass=0,supressTwig =1,supressIsolated=0,showSupressIsolated=1,nbIt=3;
var TwigID,showOnlyFinal=1,showLoops=0,showMesh=1,showNodes=1,showExtrema=1,showTwig=1,showMasterSegments=1,showSegments=1,showIsolated=1,showPassChoice=3,rmspeckles=1,Speckelradius=5;
//check - change speckle radius to 10?
var supressVerySmallSegments=1,smoothEachLimbing=1; // state = 1 required for ImageJ skeletonize function
var supressSmallSegments=1, sizeExcludeSegments=30;
var finalmapBigNodesID=0, smootMapID=0;
var recordStepOption=0, step=1, stepByStep=0;
var analyseMasterTree=1, masterTreeID=0, mapMasterTreeID=0, nbMasterJunctions=0, nbMasterJunctionsFinal=0, nbMasterSgment=0, totalMasterSegmentLenght=0;
var singleAnaMapExtID=0, singleAnaMapNodeID=0, singleShowMaps=0;
// variable for tab part
var windowTabName=angStatsTable,nameOfStatTab="["+windowTabName+"]",label="",undoErease="",batchChoice=0,imageFolder="",nameInitImageorFolder="";
var objectsChoices= newArray ("Analyze HUVEC Phase Contrast","Analyze HUVEC Fluo"), objectsChoice="Analyze HUVEC Phase Contrast", objects=0, excludeIlandsSize = 20;
var menuChoices= newArray ("-","Analyze Binary Tree","-","Save Current Analysis" , "Close Current Analysis","Save Current Analysis and Close","-","Blink Overlay [b]","Hide Overlay [h]","Show Overlay [s]");

var TreatedFolderSuffix ="=AN/", treatedImageSuffix="", ext1=".tif",ext2=".TIF";
var batchStatus="[Batch Progress Window]", batchStatusWindow="Batch Progress Window",countBatchTreated=0,countBatch=0,starttime=0,fileProcessed="";
var ecranX=screenWidth, ecranY=screenHeight,onFoldeByImage=0,folderTreatedpath;
var sufficAn="-tr";
var segmentCoul = "#f00efc", branchCoul = "#23e500", extremityCoul= "#f70133", extremityEdgeCoul= "#f7eb00", junctionCoul="#0035ff", twigColour="#19e7fb", isolatedTwigCoul="#00f7fa";
var nodesCoul="#ff1303", masterJunctionCoul="#ff2d00", masterSegmentCoul="#ffd300", isolatedElementCoul="#0005fe", meshCoul="#06adfd";
var rebuiltMapId=0;
var	reBuiltElements = newArray ("segment","branch","extremity","extremityEdge","junction","twig","isolatedTwig","nodes","masterJunction","masterSegment","isolatedElement","mesh");
var	reBuiltChoices = newArray (12), reBuiltChoicesNames = newArray (12);

var errorNetMessage ="Error: ";
var urllist = "http://image.bio.methods.free.fr/ij/ijupdatetest/ListOfMacros.txt";// to check the internet access
var onlinedoclink = "http://image.bio.methods.free.fr/ImageJ/?Angiogenesis-Analyzer-for-ImageJ.html";

var demoimagelink1 = "http://image.bio.methods.free.fr/ij/ijmacro/Angiogenesis/HUVEC-Pseudo-Phase-Contrast.tif.zip";
var demoimagelink2 = "http://image.bio.methods.free.fr/ij/ijmacro/Angiogenesis/HUVEC-Fluo.tif.zip";

var demoimagename1 = "HUVEC-Pseudo-Phase-Contrast.tif";
var demoimagename2 = "HUVEC-Fluo.tif";

var xx = requires147d(); // check version at install time
function requires147d() {requires("1.47d"); return 0; }

var concatMenuChoices=Array.concat(objectsChoices,menuChoices);
var menuanalysis=0;
var dCmds = newMenu("Network Analysis Menu Tool",concatMenuChoices);
//sectionEnd>

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//-------------------        MOVA Tools      ----------------------
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



//****************************Complete Process Macro Tool
//Variables
	//Input Settings
	var title="No-Title";
	var typeOptions = newArray("Phase Contrast", "OCT","Binary Mask", "Original Masking");//Supported image inputs - I cna probebly remove more here
	var imgType=typeOptions[2]; //Type of image used for masking
	var plotList;
	var dir = "C:/Users/rm00955/AppData/Local/Temp/",imgID,roiFolder = "ROIs",imgFolder = "Images",stackFolder = "Stacks",csvFolder = "CSVs",summaryFolder = "Summary",graphFolder = "Graphs";

	//Dialog
//<<<<<<< HEAD
	var isComprase=1,	compraseRatio=512, isAutoTh=1, imgSave=1, roiSave=1;

	//Key variables
	var thKeys = newArray("Junction","Master-Junction","Node");
	var thKeysDefault = newArray(0,0,0,0,0,1,1,0,0,1,0,0);
	var interestKeys = newArray("Branch","Isolated-Element","Master-Segment");
	var interestKeysDefaults = newArray(1,0,0,1,0,0,0,1,0,0,0,0);

	//Output options
	var resultDialog= newArray("Generate Result Stacks","Generate Individual Results Tables","Generate Summary Table","Generate Complete Results Table","Generate Skater Plot");
	var resultOptions = newArray("Stacks","Results Tables","Summary Table","Complete Table","Skater Plot");
	var optionsGenList, optionsSaveList;

	//Individual File options
	var fileTypes = newArray("jpg","tif");
	var possibleOutputs = newArray("CompleteMask","ProcessPlots","LocalThickness","Skeleton","Network-RoiSet","Overlay","SegmentMask","ThicknessDistribution");
	var possibleOutputsDefaults = newArray(1,0,1,1,1,1,0,1,1);
	var possibleOutputsClose = newArray(0,1,0,0,0,0,1,0,1);
	var outputSaveList,outputTypeList,outputCloseList;



	//****************************Macro
	macro "MOVA" {
		isCMP = 0;//if true cancels exit batch mode on ang analyser.
		//Settings
		setOption("ExpandableArrays", true);
		run("Input/Output...", "jpeg=85 gif=-1 file=.csv use_file copy_row save_column save_row");

		//Variables
		scaterPlot = newArray();
		imgTitleOnlyList = newArray();
		yLimits = newArray();
		plotList = "";

		//Start
		//Read settings sent from MATLAB
		args=getArgument();
		fname = substring(args, (indexOf(args,"fname=")+6), indexOf(args,"Savedir="));
		dir = substring(args, (indexOf(args,"Savedir=")+8), indexOf(args,"Manual"));
		ManCur = substring(args, (indexOf(args,"Manual=")+7), indexOf(args,"TwigSize="));
		twigSize = substring(args, (indexOf(args,"TwigSize=")+9));
		checkResultFolders(dir);
		fnameOCT = fname + "-original";

		//Process image
			//System Settings
			setOption("BlackBackground", 1);
			imgID = 1;

			//Image Settings
			orgImg = fname;
			//Initiate
				selectWindow(orgImg);

				getPixelSize(orgUnit, orgPixelWidth, pixelHeight);
				run("Set Scale...", "distance=0 known=0 unit=pixel");

				//IMG info - I can probably streamline a bit more here!
				title = orgImg;
				imgTitleOnlyList[1] = title;
				getDimensions(orgWidth, height, channels, slices, frames);

				//Check if binary mask
				if (is("binary")) {
					//GU Assuming images are already compressed images --> can you actually tell me what the picture sizes are? You might want to reduce the size that would significantly speedup the network processing.
						compressedOrg = orgImg;

					//Create Duplicate
					selectWindow(orgImg);
					run("Duplicate...", "title=["+orgImg+"-CompleteMask]");
					maskingImg = getTitle();
					// //Check Background Colour
					// selectWindow(maskingImg);
					// if (getValue("color.background")!=0) {
					// 		//run("Invert LUT");
					// 		//run("Invert");
					// 		print("TRUE");
					// }
					//Get Relative Area
					run("Clear Results");
					selectWindow(maskingImg);
						run("Set Measurements...", "area area_fraction redirect=None decimal=3");
						run("Measure");
						selectWindow("Results");
							totalArea = getResult("Area", 0);
							maskAreaFunction = getResult("%Area", 0);
							run("Close");
					//Save Image
						maskingImg = checkAndSave(dir+imgFolder+"/",maskingImg,".tif",1,0);
				} else exit(title+"\nNot Binary Picture");

			//Images Based on Final Mask
			//<sectionStart


//Images Based on Final Mask
			selectWindow(maskingImg);
			run("Set Scale...", "distance=0 known=0 unit=pixel");

			//Make selection Based on Mask
			roiManager("reset");
			selectWindow(maskingImg);
			run("Create Selection");
			roiManager("Add");
			roiManager("Select", 0);
			roiManager("Rename", "CompleteMask");
			roiManager("Select", 0);
			maskROI = title+"_MaskROI";
			maskROI = checkAndSave(dir+roiFolder+"/",maskROI,".roi",1,roiSave);

				//Mask Accuracy
				selectWindow(maskingImg);
				setBatchMode(1);
				overlayImg = title+"-MaskAccuracy";
				overlayImg = overlayMask(compressedOrg,overlayImg,1);
				overlayImg = checkAndSave(dir+imgFolder+"/",overlayImg,".jpg",1,0);
				setBatchMode("exit and display");


				//Thickness
				selectWindow(maskingImg);
				Roi.remove
				thicknessImg = maskToDiameter(maskingImg,imgType); //Check
				thicknessImg = checkAndSave(dir+imgFolder+"/",thicknessImg,".tif",1,imgSave);


					//Measue Average Thickness based on ROI
					roiManager("Select", 0);
					run("Set Measurements...", "mean standard modal min perimeter integrated median redirect=None decimal=3");
					roiManager("Measure");
					//Copy Table
					setOption("CopyHeaders", 0);
					String.copyResults;
					maskStats = String.paste;
					//Close results table
					selectWindow("Results");
					run("Close");



				//Skeletonize
				skeletonImg = makeSkeleton(maskingImg);
				skeletonImg = checkAndSave(dir+imgFolder+"/",skeletonImg,".tif",1, imgSave);
				//Run Angiogenesis Analyzer
				networkImg = analyzeTree(skeletonImg);

				//Get Components to ROI and extract complete list
				componentList = networkToROI(networkImg,ManCur);


				setBatchMode(1);
				overlayImg2 = title+"-Overlay";
				overlayImg2 = overlayMask(fnameOCT,overlayImg2,1);
				overlayImg2 = checkAndSave(dir+imgFolder+"/",overlayImg2,".jpg",1,imgSave);
				setBatchMode("exit and display");

				//Save ROI (And update components list???)
				roiSetSave = title+"_Network-RoiSet";
				roiSetSave = checkAndSave(dir+roiFolder+"/",roiSetSave,".zip",1,roiSave);
				close(networkImg);

			//Threshold based on nodes
				setBatchMode(1);
				segmentMask = keyBasedTH(thicknessImg,componentList,thKeys);
				thKeyValue = segmentMask[1];
				segmentMask = segmentMask[0];//
				segmentMask = checkAndSave(dir+imgFolder+"/",segmentMask,".tif",1,0);
					//Get Relative Area
					run("Clear Results");
					selectWindow(segmentMask);
						run("Set Measurements...", "area area_fraction redirect=None decimal=3");
						run("Measure");
						selectWindow("Results");
							segmentMaskAreaFunction = getResult("%Area", 0);
							run("Close");
				setBatchMode("exit and display");

			//Make Results table
				thTableTitle = title+"_Local-Diameter_Table";
				thTable = measureROIs(thicknessImg,componentList,interestKeys,thTableTitle,"mean standard modal min perimeter bounding median limit");

			//Get Average values for Summary
				selectWindow(thTableTitle);
				getAverages = newArray("Average Thickness","Adjusted Len","Endpoint Distance","Tortuosity");
				valuesString = "";
				for (k = 0; k < getAverages.length; k++) {
					myAverageArray = newArray();
					l = 0;
					averageValue = Table.getColumn(getAverages[k]);
					for (j = 0; j < averageValue.length; j++) {
						if (!isNaN(averageValue[j])) {
							myAverageArray[l] = averageValue[j];
							l++;
						}
					}
					Array.getStatistics(myAverageArray, min, max, myMean, stdDev);
					valuesString = valuesString+"\t"+myMean;
				}



			//Close Part Table Tables
				thTableTitle = checkAndSave(dir+csvFolder+"/",thTableTitle,".csv",1,1);
				print(thTableTitle); //Maybe and error here! Check!


				selectWindow(angStatsTable);
					run("Close");


	}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//---------------------------        Functions      -------------------------
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//Make Stats
function measureROIs(img,list,keys,table,measurs) {
	//Reset Results
	run("Clear Results");
	selectWindow(img);
	//Variables
	categoryColumn = newArray();

	//Measure ROIs
		for (i = 0; i < keys.length; i++) {
			//Measure ROI
			cROI = listToArray(list,keys[i]);//Get Key Specific ROI
			roiManager("Select",cROI);
			if(cROI.length>0) {
				run("Set Measurements...", measurs+" display redirect=["+img+"] decimal=3");
				roiManager("Measure");
			}

			//Key position in components list
			for (j = 0; j < componentKeys.length; j++) {
				if (componentKeys[j] == keys[i]) categoryValue = j;
			}

			//Categorise ROI
			categoryArray = newArray(cROI.length);
				Array.fill(categoryArray, categoryValue);
				categoryColumn = Array.concat(categoryColumn,categoryArray);
		}

		//Rename Category values
		categoryColmnName = newArray();
		for (i = 0; i < categoryColumn.length; i++) {
			categoryColmnName[i] = componentKeys[categoryColumn[i]];
		}

		//Creat New Length
			myLenArray = newArray();
			myStraightLine = newArray();
			myTortuosity = newArray();

			//Get Measurements
			selectWindow("Results");
			headers = Table.headings;
			isLen = (indexOf(headers, "Length")!=-1);
				perimeterArray = Table.getColumn("Perim.");
				widthArray = Table.getColumn("Width");
				heightArray = Table.getColumn("Height");
				if (isLen) lengthArray = Table.getColumn("Length");

			//Calculate New Values
			for (i = 0; i < perimeterArray.length; i++) {
				if (isLen){
				 	if(perimeterArray[i]==lengthArray[i]) {
						myLenArray[i] = lengthArray[i];
					} else myLenArray[i] = perimeterArray[i]/2;
				} else myLenArray[i] = perimeterArray[i]/2;
				//Angiogenesis analyser uses this value to calculate the total lengths
				myStraightLine[i] = Math.sqrt(Math.sqr(widthArray[i])+Math.sqr(heightArray[i]));
				myTortuosity[i] = myLenArray[i]/myStraightLine[i];
			}

		//Set new Table with all the components and category
		selectWindow("Results");
			imgIDColumn = newArray(categoryColumn.length);
			Array.fill(imgIDColumn,imgID);
			if (isLen) Table.deleteColumn("Length");
			Table.setColumn("Image ID", imgIDColumn);
			Table.setColumn("Component (code)", categoryColumn);
			Table.setColumn("Component", categoryColmnName);
			Table.setColumn("Adjusted Len", myLenArray);
			Table.setColumn("Endpoint Distance", myStraightLine);
			Table.setColumn("Tortuosity", myTortuosity);
			Table.renameColumn("Mean", "Average Thickness");
				Table.update();

			//Copy Table
			setOption("CopyHeaders", 0);
			String.copyResults;
				resultString = String.paste;

			//Rename Table
			Table.rename("Results", table);

	return resultString;
}

//make plot
function generatePlot(plotTitle,table,x,y,sortBy,type,highRes) {
	//Variables
		tmpPlot = "tmpPlots";
		rstTable = plotTitle+"_Values";
	//Get values
	selectWindow(table);
	yValues = Table.getColumn(y);
	if (type!="Single") {
		xValues = Table.getColumn(x);
		run("Plots...", "width=600 height=340 font=14 draw_ticks list minimum=0 maximum=0 interpolate");
	} else {
		xValues = newArray(yValues.length);
		Array.fill(xValues,1);
		run("Plots...", "width=200 height=340 font=14 draw_ticks list minimum=0 maximum=0 interpolate");
	}
	//X max
	Array.getStatistics(xValues, xMin, xMax, mean, stdDev);

	//Get Categories
	catergories = Table.getColumn(sortBy);
	Array.sort(catergories);
	unique = newArray();
	j = 1;
	unique[0] = catergories[0];

	for (i = 1; i < catergories.length; i++) {
		if (catergories[i]!=catergories[i-1]) {
			unique[j]  = catergories[i];
			j++;
		}
	}

	//Create Plot
	selectWindow(table);
	catergories = Table.getColumn(sortBy);
	legend = "";
	Plot.create(plotTitle, "Components", y);
	colours = newArray("balck","green","red","blue","purple");
	maxArray = newArray();

	for (i = 0; i < unique.length; i++) {
		yPlot = newArray();
		xPlot = newArray();
		k = 0;
		for (j = 0; j < yValues.length; j++) {
			if (catergories[j]==unique[i]) {
				yPlot[k] = yValues[j];
				xPlot[k] = xValues[j];
				k++;
			}
		}
		if (type=="Single") {
			legend = legend+unique[i]+" n="+xPlot.length+"\n";
		} else {
			legend = legend+unique[i]+"\n";
		}
		Plot.add(unique[i], xPlot, yPlot);
		Plot.setStyle(i, colours[i]+",none,2.0,X");
		maxArray = Array.concat(maxArray,yPlot);
	}
	Array.getStatistics(maxArray, min, yMax, mean, stdDev);
	Plot.addLegend(legend, "Auto");
	Plot.setLimits(xMin-0.5, xMax+0.5, 0, yMax*1.1);
	Plot.setLimitsToFit();
	Plot.show();

	selectWindow("Plot Values");
		Table.rename("Plot Values",rstTable);
		selectWindow(rstTable);
		headingsString = Table.headings;
			headings = split(headingsString,"\t");

			evenHeadings = headings.length==unique.length*2;
			if(evenHeadings) {
				k=0;
				for (i = 0; i < unique.length; i++) {
					if (indexOf(headingsString,headings[k])!=-1) {
						Table.renameColumn(headings[k], "X-"+unique[i]);
					} else print("ERROR: No column "+headings[k]+" found in 'Headings' of Table: "+rstTable);
					k++;
					if (indexOf(headingsString,headings[k])!=-1) {
						Table.renameColumn(headings[k], "Y-"+unique[i]);
					} else print("ERROR: No column "+headings[k]+" found in 'Headings' of Table: "+rstTable);
					k++;
				}
			}

	//Create High Res image
	if (highRes) {
		selectWindow(plotTitle);
			rename(tmpPlot);
		Plot.makeHighResolution(tmpPlot,4.0,"disable");
			rename(plotTitle);
		selectWindow(tmpPlot);
			run("Close");
	}

	return rstTable;
}

//Find Unique values
function findUniqueValues(array) {
	Array.sort(array);
	unique = newArray();
	j = 1;
	unique[0] = array[0];

	for (i = 1; i < array.length; i++) {
		if (array[i]!=array[i-1]) {
			unique[j]  = array[i];
			j++;
		}
	}

	return unique;
}

//Check File in folder
function getSaveName(dir, titleWoExt, ext) {
	resultName = titleWoExt;
	file = dir+titleWoExt+ext;
		fileExists = File.exists(file);
		replicate = 0;

	//AlterFilename until there is no overwrite
		while (fileExists) {
			replicate += 1;
			resultName = titleWoExt+"-"+replicate;
			file = dir+resultName+ext;
			fileExists = File.exists(file);
		}

 	return resultName;
}

//Check not to overwrite and Save
function checkAndSave(saveFolder,windowName,ext,check,doSave) {
	//Save File
	saveName = windowName;

	if (doSave) {
		if (check) saveName = getSaveName(saveFolder, saveName,ext);

		if (ext==".zip" || ext==".roi") {
			roiManager("Save", saveFolder+saveName+ext);
		} else {
			selectWindow(windowName);
			if (ext==".csv") {
				saveAs("Results", saveFolder+saveName+ext);
				saveName = getInfo("window.title");
				selectWindow(saveName);
				Table.setLocationAndSize(0,0,500,500,saveName);
			} else if (ext==".jpg") {
				saveAs("Jpeg", saveFolder+saveName+ext);
				saveName = getTitle;
				selectWindow(saveName);
				setLocation(0, 0, 500, 500);
			} else {
				saveAs("Tiff", saveFolder+saveName+ext);
				saveName = getTitle;
				selectWindow(saveName);
				setLocation(0, 0, 500, 500);
			}
		}
	}
	return saveName;
}

//Compress image - Check - I am pretty sure I can delete this
function compressImg(img,newTitle,ref,useImg) {
	selectWindow(img);
		if (!useImg) run("Duplicate...", "title=["+newTitle+"]");
			run("Size...", "width="+ref+" depth=1 constrain average interpolation=Bilinear");

	return newTitle;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//-----------------------        CMP Functions    ---------------------------
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//--------------------        General Functions      ------------------------
//Functions just to help with the processing in general

//****************************Background Functions
//<sectionStart
//Array to String
	function arrayToString(array,delimeter) {
		n = array.length;
		string = ""+array[0];

		if (n>1) {
			for (i = 1; i < n; i++) {
				string += delimeter+array[i];
			}
		}

		return string;
	}

//Close Image-array
	//Closes all image in input array
	function closeImgs(imgArray) {
		for (i = 0; i < imgArray.length; i++) {
			close(imgArray[i]);
		}
	}

//Get Array from list
	//Converts extracts an array from a list assigned to a key.
	function listToArray(list,key) {
		List.setList(list);
		myArray = List.get(key);

		//Check if components had values
		isArray = lastIndexOf(myArray,",");

		//Convert component string to array
		if(isArray!=-1) {
			myArray = substring(myArray, 0, isArray);
			myArray = split(myArray,",");
		} else myArray = newArray();

		return myArray;
	}

//Get Value from list
	function listToValue(list,key) {
		List.setList(list);
		myValue = List.get(key);

		return myValue;
	}

//Get Dimensions
	//Gets the longest side of the image and the pixel width and unit - returns array
	function getmyDim(img) {
		selectWindow(img);
		getDimensions(width, height, channels, slices, frames);
			if (width>height) {
				maxSide=width;
			} else maxSide=height;
		getPixelSize(unit, pixelWidth, pixelHeight);

		return newArray(width,height,maxSide,pixelWidth,unit);
	}

//Check Result Folders
	//Check if subfolders are there
	function checkResultFolders(dir) {
		folders = newArray(roiFolder,imgFolder,stackFolder,csvFolder,summaryFolder,graphFolder);
		//Create Subfolders
			for (i = 0; i < folders.length; i++) {
				if (!File.isDirectory(dir+folders[i])) File.makeDirectory(dir+folders[i]);
			}
	}
//sectionEnd>

//****************************Image Processing Functions
//<sectionStart
//Mask to diameter map
	function maskToDiameter(img,maskType) {
		rstImg = title+"-LocalThickness";
		tmpImgs = newArray("tmp");

		selectWindow(img);
			run("Duplicate...", "title=tmp");
			if (maskType==typeOptions[3]) run("Invert");
			run("Local Thickness (masked, calibrated, silent)");
			setMinAndMax(0, 255);
			call("ij.ImagePlus.setDefault16bitRange", 8);
			rename(rstImg);

		//Close tmpImgs
		closeImgs(tmpImgs);

		return rstImg;
	}

//Skeletonize
	//converts the mask image in to a skelton
	function makeSkeleton(img) {
		rstImg = title+"-Skeleton";
		//Make Skeleton
		selectWindow(img);
			run("Duplicate...", "title=["+rstImg+"]");
			run("Skeletonize");
			run("Convert to Mask");
			run("Invert LUT");
			//run("Invert");

	return rstImg;
	}

//Mask Overlay
	//Combines the original image woth the mask to indicate accuracy of the process
	function overlayMask(img,overlay,roi) {
		rstImg = overlay;
		tmpImgs = newArray("Outline","Combined-8bit","Composite");

		//Processing

		selectWindow(img);
			run("Duplicate...", "title=Combined-8bit");
			run("RGB Color");
			run("8-bit");



		//Overlay ROI
		if (roi) {
			selectWindow(img);
				roiManager("Show All without labels");
				run("Flatten");
				rename("withROI");

			selectWindow("withROI");
				rename(rstImg);
		}

		//Close tmpImgs
		closeImgs(tmpImgs);

		return rstImg;
	}

//sectionEnd>

//------------------        Network Analysis Functions      ------------------
//Combining CMP macros with the angiogenesis analyser and translating information

//****************************Network Processing
//<sectionStart GU - this is maybe the part causing an eror? Maybe comment out, but test!
//Complete Angiogenesis Analysis
	//Runs the complete process of angiogenesis analyser
	//It assumes that the picture is a RGB picture
	function runAngAnalysis(img) {
		rstImg = title+"-AngAnalyser"

		//Run Angiogenesis Analyser single network analysis
		objects = 1;
		menuanalysis=1;
		singleNetworkAnalysis();

	return rstImg;
	}

//Tree analysis
	//Analyse the skeleton image based on the angiogenesis analyser
	//returns image with overlay
	function analyzeTree(img) {
		//My Results name
		rstImg = title+"-Network";
		//Angiogenesis analyser name
		if (lastIndexOf(img, ".") > 1) {
			angName = substring(img,0,lastIndexOf(img, "."));
			angName = angName+"-tr";
		} else {
			angName = img+"-tr";
		}

		//Run Angiogenesi Analyser
		selectWindow(img);
		analyseFromATree();
			rename(rstImg);

		//Change Input image name back to original
		selectWindow(angName);
			rename(img);

	return rstImg;
	}

//Extract Network
	//Generate ROI out of the image with
	//colour coded overlay based on the angiogenesis analyser.
	function networkToROI(img,userInput) {
		var segmentCount=0, branchCount=0, extremityCount=0, extremityCountEdge=0, junctionCount=0,
			twigCount=0, nodeCount=0, masterJunctionCount=0, masterSegmentCount=0,
			isolatedElementCount=0, meshCount=0, isolatedTwigCount=0;

		segmentArray=""; branchArray=""; extremityArray=""; extremityEdgeArray=""; junctionArray="";
		twigArray=""; nodeArray=""; masterJunctionArray=""; masterSegmentArray="";
		isolatedElementArray=""; meshArray=""; isolatedTwigArray="";

		//Get Overlay from image
		selectWindow(img);
			run("To ROI Manager");
			nROI = roiManager("count");
			if (nROI==0) exit("No ROI");

		//Assign new names and categories to each ROI based on colour code
		for (i=0; i< nROI; i++) {
			roiManager("select",i);
			colorSelection = getInfo("selection.color");
			cathegory="";
			if (toString(colorSelection) == toString(extremityCoul)) {
				extremityArray = extremityArray+roiManager("index")+",";
				extremityCount++;
				cathegory="extremity";
				name="extremity-"+extremityCount;
				roiManager("Rename",name);
			}
			if (toString(colorSelection) == toString(extremityEdgeCoul)) {
				extremityEdgeArray = extremityEdgeArray+roiManager("index")+",";
				extremityCountEdge++;
				cathegory="extremityEdge";
				name="extremityEdge-"+extremityCountEdge;
				roiManager("Rename",name)
			}
			if (toString(colorSelection) == toString(junctionCoul)) {
				junctionArray = junctionArray+roiManager("index")+",";
				junctionCount++;
				cathegory="junction";
				name="junction-"+junctionCount;
				roiManager("Rename",name);
			}
			if (toString(colorSelection) == toString(twigColour)) {
				twigArray = twigArray+roiManager("index")+",";
				twigCount++;
				cathegory="twig";
				name="twig-"+twigCount;
				roiManager("Rename",name);
			}
			if (toString(colorSelection) == toString(segmentCoul)) {
				segmentArray = segmentArray+roiManager("index")+",";
				segmentCount++;
				cathegory="segment";
				name="segment-"+segmentCount;
				roiManager("Rename",name);
			}
			if (toString(colorSelection) == toString(branchCoul)) {
				branchArray = branchArray+roiManager("index")+",";
				branchCount++;
				cathegory="branch";
				name="branch-"+branchCount;
				roiManager("Rename",name);
			}
			if (toString(colorSelection) == toString(nodesCoul)) {
				nodeArray = nodeArray+roiManager("index")+",";
				nodeCount++;
				cathegory="nodes";
				name="nodes-"+nodeCount;
				roiManager("Rename",name);
			}
			if (toString(colorSelection) == toString(masterJunctionCoul)) {
				masterJunctionArray = masterJunctionArray+roiManager("index")+",";
				masterJunctionCount++;
				cathegory="masterJunction";
				name="masterJunction-"+masterJunctionCount;
				roiManager("Rename",name);
			}
			if (toString(colorSelection) == toString(masterSegmentCoul)) {
				masterSegmentArray = masterSegmentArray+roiManager("index")+",";
				masterSegmentCount++;
				cathegory="masterSegment";
				name="masterSegment-"+masterSegmentCount;
				roiManager("Rename",name);
			}
			if (toString(colorSelection) == toString(isolatedElementCoul)) {
				isolatedElementArray = isolatedElementArray+roiManager("index")+",";
				isolatedElementCount++;
				cathegory="isolatedElement";
				name="isolatedElement-"+isolatedElementCount;
				roiManager("Rename",name);
			}
			if (toString(colorSelection) == toString(meshCoul)) {
				meshArray = meshArray+roiManager("index")+",";
				meshCount++;
				cathegory="mesh";
				name="mesh-"+meshCount;
				roiManager("Rename",name);
			}
			if (toString(colorSelection) == toString(isolatedTwigCoul)) {
				isolatedTwigArray = isolatedTwigArray+roiManager("index")+",";
				isolatedTwigCount++;
				cathegory="isolatedTwig";
				name="isolatedTwig-"+isolatedTwigCount;
				roiManager("Rename",name);
			}
			name = "";
		}

		//User Curraction Process
		if (userInput){
			nRoiUserStart = roiManager("count"); // ROI number before Start of curation

			//User Curation
			selectWindow(img);
				roiManager("Show All without labels");
				setTool("freeline");
			waitForUser("Manual Curation", "Happy with ROIs?\nPress 't' to add new selection to ROI list\nRemove unwanted ROIs in ROI manager.");

		//Rename User defined ROIs
			nRoiUserEnd = roiManager("count"); // ROI number after Start of curation

			for (i = nRoiUserStart+1; i<nRoiUserEnd; i++) {
				roiManager("select",i);
				branchArray = branchArray+roiManager("index")+",";
				branchCount++;
				cathegory="branch";
				name="branch-"+branchCount+"-USER";
				roiManager("Rename",name);
			}

		//Redefine Branches Array
			segmentCount=0; branchCount=0; extremityCount=0; extremityCountEdge=0; junctionCount=0;
				twigCount=0; nodeCount=0; masterJunctionCount=0; masterSegmentCount=0;
				isolatedElementCount=0; meshCount=0; isolatedTwigCount=0;

			segmentArray=""; branchArray=""; extremityArray=""; extremityEdgeArray=""; junctionArray="";
			twigArray=""; nodeArray=""; masterJunctionArray=""; masterSegmentArray="";
			isolatedElementArray=""; meshArray=""; isolatedTwigArray="";

			for (i=0; i< nRoiUserEnd; i++) {
					roiManager("select",i);
					roiTitle = Roi.getName;
					cathegory="";
					if (lastIndexOf(roiTitle,"extremity")!=-1) {
					extremityArray = extremityArray+roiManager("index")+",";
					extremityCount++;
					cathegory="extremity";
					name="extremity-"+extremityCount;
					roiManager("Rename",name);
					} else if (lastIndexOf(roiTitle,"extremityEdge")!=-1) {
					extremityEdgeArray = extremityEdgeArray+roiManager("index")+",";
					extremityCountEdge++;
					cathegory="extremityEdge";
					name="extremityEdge-"+extremityCountEdge;
					roiManager("Rename",name)
					} else if (lastIndexOf(roiTitle,"junction")!=-1) {
						junctionArray = junctionArray+roiManager("index")+",";
						junctionCount++;
						cathegory="junction";
						name="junction-"+junctionCount;
						roiManager("Rename",name);
					} else if (lastIndexOf(roiTitle,"twig")!=-1) {
						twigArray = twigArray+roiManager("index")+",";
						twigCount++;
						cathegory="twig";
						name="twig-"+twigCount;
						roiManager("Rename",name);
					} else if (lastIndexOf(roiTitle,"segment")!=-1) {
						segmentArray = segmentArray+roiManager("index")+",";
						segmentCount++;
						cathegory="segment";
						name="segment-"+segmentCount;
						roiManager("Rename",name);
					} else if (lastIndexOf(roiTitle,"branch")!=-1) {
						branchArray = branchArray+roiManager("index")+",";
						branchCount++;
						cathegory="branch";
						name="branch-"+branchCount;
						roiManager("Rename",name);
					} else if (lastIndexOf(roiTitle,"nodes")!=-1) {
						nodeArray = nodeArray+roiManager("index")+",";
						nodeCount++;
						cathegory="nodes";
						name="nodes-"+nodeCount;
						roiManager("Rename",name);
					} else if (lastIndexOf(roiTitle,"masterJunction")!=-1) {
						masterJunctionArray = masterJunctionArray+roiManager("index")+",";
						masterJunctionCount++;
						cathegory="masterJunction";
						name="masterJunction-"+masterJunctionCount;
						roiManager("Rename",name);
					}else if (lastIndexOf(roiTitle,"masterSegment")!=-1) {
						masterSegmentArray = masterSegmentArray+roiManager("index")+",";
						masterSegmentCount++;
						cathegory="masterSegment";
						name="masterSegment-"+masterSegmentCount;
						roiManager("Rename",name);
					}else if (lastIndexOf(roiTitle,"isolatedElement")!=-1) {
						isolatedElementArray = isolatedElementArray+roiManager("index")+",";
						isolatedElementCount++;
						cathegory="isolatedElement";
						name="isolatedElement-"+isolatedElementCount;
						roiManager("Rename",name);
					}else if (lastIndexOf(roiTitle,"mesh")!=-1) {
						meshArray = meshArray+roiManager("index")+",";
						meshCount++;
						cathegory="mesh";
						name="mesh-"+meshCount;
						roiManager("Rename",name);
					} else if (lastIndexOf(roiTitle,"isolatedTwig")!=-1) {
						isolatedTwigArray = isolatedTwigArray+roiManager("index")+",";
						isolatedTwigCount++;
						cathegory="isolatedTwig";
						name="isolatedTwig-"+isolatedTwigCount;
						roiManager("Rename",name);
					} else {
						branchArray = branchArray+roiManager("index")+",";
						branchCount++;
						cathegory="branch";
						name="branch-"+branchCount+"-USER";
						roiManager("Rename",name);
					}
					name = "";
				}
		}


		componentsArray = newArray(
			branchArray,
			extremityArray,
			extremityEdgeArray,
			isolatedElementArray,
			isolatedTwigArray,
			junctionArray,
			masterJunctionArray,
			masterSegmentArray,
			meshArray,
			nodeArray,
			segmentArray,
			twigArray
		);

		List.fromArrays(componentKeys, componentsArray);
			networkResults = List.getList;

			numberOfNodes = junctionCount;
			print(numberOfNodes);

		return networkResults;
	}

//Set Threshold based on input keys
	//Used to set threshold based on
	function keyBasedTH(img,list,keys) {
		//Input descriptions
			//list - list containing the ROIs for components as keys.
			//img - image to perform analysis on
			//keys - list of keys that have to be included in the threshold calibration

		//Reset Results
			run("Clear Results");

		//Variables
			rstImg = title+"-SegmentMask";
			allROIs = newArray();

		//Generate ROI array
			for (i = 0; i < keys.length; i++) {
				cROI = listToArray(list,keys[i]);
				allROIs = Array.concat(allROIs,cROI);
			}

		//Measure selected Only eg.
			selectWindow(img);
			roiManager("select", allROIs);
			run("Set Measurements...", "mean display nan redirect=["+img+"] decimal=3");
					roiManager("Measure");

		//Get Values from Mean Column to measure the threshold
			selectWindow("Results");
				//
				keyMeansArray = newArray();
				l = 0;
				keyMeans = Table.getColumn("Mean");
				for (j = 0; j < keyMeans.length; j++) {
					if (!isNaN(keyMeans[j])) {
						keyMeansArray[l] = keyMeans[j];
						l++;
					}
				}
				//
				Array.getStatistics(keyMeansArray, min, max, keyAvg, stdDev);
				run("Close");

		//Set Threshold Based on the Nodes
			selectWindow(img);
			getStatistics(area, mean, thMin, max, std, histogram);
			if (isAutoTh) setThreshold(thMin, keyAvg);

		//Convert To Mask
			selectWindow(img);
			run("Duplicate...", "title=["+rstImg+"]");
			setThreshold(thMin, keyAvg);
				run("Convert to Mask");
				setAutoThreshold("Default");

		return newArray(rstImg,keyAvg);
	}

//sectionEnd>


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//--------        Angiogenesis Analyser Original Plugin      ----------------
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//<sectionStart
// Angiogenesis Analyzer
// Author : Gilles Carpentier
// Faculte des Sciences et Technologie,
// Universite Paris Est Creteil Val de Marne, France.

//**Angiogenesis Original Variable Poz**

macro "Network Analysis Menu Tool - C000D00D01D02D03D08D09D0aD0bD0cD0dD0eD0fD10D11D12D13D14D1aD1bD1cD1dD1eD20D21D22D23D24D2dD30D32D33D40D41D42D48D49D50D51D56D57D58D59D5fD66D67D68D6fD74D75D76D77D84D85D86D8cD90D96D9bD9cDa0DaaDabDacDadDb9DbaDbbDbcDbdDc0Dc1Dc2Dc3Dc8Dc9DcaDcbDccDceDd0Dd1Dd2Dd3Dd8Dd9DdaDdbDe0De1De2De8De9DeaDebDf0Df1Df2Df8Df9DfaDfbC002D69D83C000D31C573Dc6C212De7C021Dc7Cc84Da7C005Da9C001D18D5aDbeC0c8D5eC763D29Cdb0DafCa78D5dC004D07Dc4DcfC001D19D65D78D7fD87C695D7dC009D4fC430Dd7Cbf1D6bC015Db8C111DecC399D17C729D26Ccf0D4eCfffDddDdeDdfDedDeeDefDfdDfeDffC003D25D7cC000Db2C085D2eC007D94C050Db1Cbd5D82C006Da5C6a8D6cC269D97Cdb0Df4Cab6D06C012D95C6b6D93C02aDf5Ca91Df6Cfd2D89C017D2bD34C111DcdC1a9D46Cf38Db5Cee0D70C002D04D4aC685D81C034Dd6Ccc1D16C006D43D60D9aDe3C597D9fC216D2aCf87D71C686Da8C029D55DaeC640DbfCad4D6dC016D15C07bDf3Cf27D5cCdf0D3bC003D8bC000DfcC4f4D38C058D4bC070Db0Cfd2D98C8f6D39C53aDe6Ccc1D3fCac5D7eC024Da1Caf4D3aD3cC05aDf7C8c1D28Cff2D7aC017D80C6c8D8aCf69D27Cdf0D63C286D35C044D5bC005D2cD61C4d8D7bC555DdcC8a3Da4C029Db3Caf2D92C025D8dC598D53Ce4aD4cC395D6eC037Da6Cdd3Da3C9b6De4C26bD88Cca8Db6C014D47C8a6D9eC14aD79Cff1Dc5C0faD3dCf47D4dC576D62C136D91Cec2D45C697D99C349D6aCe8cD37C996Db4Cce3D8eC025D1fC29fDb7Cf27De5C004D52D9dC5f4D3eC168D64Cfd2D54C8f7D2fC36aDa2Ccf6Dd5C054Dd4C259D73Cef3D36C027D8fC8a8D44Cf2eD72C7b9D05"{
	cmd = getArgument();
	if (cmd!="-" && cmd != "Close Current Analysis" && cmd != "Save Current Analysis and Close" && cmd != "Save Current Analysis" && cmd != "Analyze Binary Tree" && cmd != "Blink Overlay [b]" && cmd != "Hide Overlay [h]" && cmd != "Show Overlay [s]") {
		objectsChoice = cmd;
		menuanalysis=1;
		if (objectsChoice == objectsChoices[0]) objects = 1;
		if (objectsChoice == objectsChoices[1]) objects = 2;
		singleNetworkAnalysis ();
	}
	if (cmd!="-" && cmd == "Close Current Analysis") {closeAll ();}
	if (cmd!="-" && cmd == "Save Current Analysis") {saveAndCloseAnalysis (0);}
	if (cmd!="-" && cmd == "Save Current Analysis and Close") {saveAndCloseAnalysis (1);}
	if (cmd!="-" && cmd == "Analyze Binary Tree") {analyseFromATree ();}
	if (cmd!="-" && cmd == "Blink Overlay [b]") {blinkOverlay ();}
	if (cmd!="-" && cmd == "Hide Overlay [h]") {run("Hide Overlay");}
	if (cmd!="-" && cmd == "Show Overlay [s]") {run("Show Overlay");}
	menuanalysis=0;
}

macro "Blurred Mask Action Tool - C40bD00D01D02D0aD0bD0cD0dD0eD0fD10D11D1bD1cD1dD1eD1fD20D2cD2dD2eD2fD30D3dD3eD3fD4eD4fD5fD6fD70D7fD80D90Da0Da1Db0Db1Dc0Dc1Dc2Dd0Dd1Dd2Dd3DdfDe0De1De2De3De4DefDf0Df1Df2Df3Df4Df5DfeDffCe30D63D8cCc04D7cCf50D5aD9cC40cD03D08D09D12D1aD21D2bD3cD40D4dD50D5eD60D6eD8fD91D9fDafDb2DbfDc3DcfDd4De5DeeDf6Df7DfcDfdC60dD13D22D2aD41D4cD51D61D92Dc4DceC50cD04D07D31D5dD81Da2DdeDe6DedDf8DfbCfc0D47D96Da7Cc04D25DccCb07D24Cf90D48D6aDa6Db7C70eD82Da3DaeC90bD42Cfe0D56D68Da9Ce40D38D53D6bDa5Cb06D5bCf80D7bC60dD18D3bD6dD8eD9eDbeDddDe7DecCfe0D55D65C80dDc5Ca0aD28Cfb0D85D9bDabC80dD14D17D29D32D5cD7dDb4Dd6De8DebC90cD3aD4bDcdDe9DeaCff3D88Cd12D73Cc06D39Dd8DdbCf70D37Dc8Cfd0D69D8aC50dD05D06D19D71D7eDb3Dd5Df9DfaCa09D9dCfa0D54D59D64D8bDbbCa0bD62D93Cff2D66D67D77D78D87D98D99Cf50D34Db6DbcCb07D33D83Db5Dc6Cf90D36D44D74D95Dc9DcaCfe0D57D79D86D97D9aDa8DaaCa0bD52D8dDd7DdcCfc0D45D7aDb8Cd21D43D94Cc05D4aCf60D49D84DacDcbCfd0D46D58D75Db9DbaC90cD15D16D23D72Cb08DadCff1D76D89Cf80D35Cd13D26Dd9DdaCa09D6cDbdCe31Dc7Cc05D27Cb08Da4"{
	setBatchMode(true);
	imageTomask=getImageID();
	if (selectionType() == -1) {
		showMessage ("This tool requires a selection");
		setTool("ellipse");
		makeEllipse(30, 30, 100, 100, 0.4);
		exit;
	}
	run("Select None");
	run("Duplicate...", "title=masqueflou");
	masqueflouID = getImageID();
	run("Restore Selection");
	run("Gaussian Blur...", "sigma=4");
	run("Enlarge...", "enlarge=10");
	run("Gaussian Blur...", "sigma=2");
	run("Enlarge...", "enlarge=10");
	run("Gaussian Blur...", "sigma=1");
	run("Select None");
	run("Select All");
	run("Copy");
	close();
	selectImage(imageTomask);
	setupUndo();
	run("Paste");run("Select None");
}

var settingsMenuChoices= newArray ("Find a Tree","Find & Remove Loops","Find Extremities","Find Nodes and Junctions","Find Nodes and Branches","Record the Steps of Limbing","-","Get Maps of Selections","-","Blink Overlay [b]","Hide Overlay [h]","Show Overlay [s]","-","Settings");
var dCmds = newMenu("Tuning Functions Menu Tool",settingsMenuChoices);
macro "Tuning Functions Menu Tool - C000D19D6eCfffD00D01D02D03D04D05D06D07D08D0bD0cD0dD0eD0fD10D11D12D13D14D15D1cD1dD1eD1fD20D21D22D23D24D2bD2cD2dD2eD2fD30D31D32D33D3aD3bD3cD3dD3eD3fD40D41D42D43D49D4aD4bD4cD4dD4eD4fD50D51D52D53D58D59D5aD5bD5cD60D61D62D63D69D6aD6bD6cD6fD70D71D72D73D79D7aD7fD80D81D82D83D8fD90D91D92D9eD9fDa0Da1Da2DaeDafDb0Db1DbdDbeDbfDc0Dc1DccDcdDceDcfDd0Dd9DdaDdbDdcDddDdeDdfDe0De8De9DeaDebDecDedDeeDefDf7Df8Df9DfaDfbDfcDfdDfeDffC99aD97C777D55CcccDc9Dd5De3C666D67CbbbD64Da3Da7Da8DbbDc8C998D44CdddD93C445D36CaaaD54C888D46D87CccdD2aDd4De4Df2C777D78CbbbDabC999D74DacDbaDc7De1De2CfffD09D39D48D5fDcbDd8De7C111D26D7dCaaaDb6Dc2C778D25D65CcccD16D1bDadDbcC667D37D75D77CbbbDa9DaaDb7Db8Dc4C999D56D66D95D9dDb9Df4CeeeD5dD68D7bDb2C555D8dCaabD9cDa6Db5Dc5C989D38D76CdddDf0C777D85D99Da4Dc3Dd2CbbcD8eDf3C111D18D1aCaaaDc6Dd7De5C778D9aC666D29D5eDf5C999D84Db4CeeeDd1Df6C555D45D88CaaaD34Da5Dd3C889D86D98C767D47D57C333D28C888D94Dd6Df1C767D8aC666D17D7cD8cCdddD0aDcaC99aD96C777Db3De6C555D6dD8bC223D35C101D27C444D7eD89Ca9aD9b"{
	cmd = getArgument();
	if (cmd!="-" && cmd == "Find a Tree") {findTree ();}
	if (cmd!="-" && cmd == "Find & Remove Loops") {FindAnfRemoveLoops ();}
	if (cmd!="-" && cmd == "Find Extremities") {findExtremities ();}
	if (cmd!="-" && cmd == "Find Nodes and Junctions") {findTheNodes ();}
	if (cmd!="-" && cmd == "Find Nodes and Branches") {
		findNodesAndBranches ();
		setBatchMode("exit and display");
		selectImage (workid);
	}
	if (cmd!="-" && cmd == "Get Maps of Selections") {getMapsOFSelections ();}
	if (cmd!="-" && cmd == "Record the Steps of Limbing") {recordTheStepOfLimbing ();}
	if (cmd!="-" && cmd == "Blink Overlay [b]") {blinkOverlay ();}
	if (cmd!="-" && cmd == "Hide Overlay [h]") {run("Hide Overlay");}
	if (cmd!="-" && cmd == "Show Overlay [s]") {run("Show Overlay");}
	if (cmd!="-" && cmd == "Settings") {globalSettings ();}
}

macro "Batch Image Treatment Action Tool - C578D1fD2fD3fD4fD5fD6fD7fD8fD9fDafDbfDcfDdfDefC9bdD23D33D43D53D89D99Da9Dc8Dd8De9C69cD1eDeeCbccD12D62C69bD3eD4eD5eD6eD7eD8eD9eDceCabdDe6C8acD24D34D44D54CbdeD95C689D05D06D07D08D09D64D74D84D94Da4Db4Dc4Dd4De4Df5Df6Df7Df8Df9CabdD88C7acD2cDccCdddD00D01D02D03D10D11D20D21D30D31D40D41D50D51D60D61D70D71D72D73D80D81D82D83D90D91D92D93Da0Da1Da2Da3Db0Db1Db2Db3Dc0Dc1Dc2Dc3Dd0Dd1Dd2Dd3De0De1De2De3Df0Df1Df2Df3C79aD13D22D32D42D52D63CacdD76C8adDebCcdfD15D25D35C689D0eDfcDfdDfeC9bdDe8C7acD8cD9cDacDbcDdcC79bD2eDaeDbeDdeCacdD37D47D57D86D96C8bcD2bD3bD4bD5bD8aD9aDaaDbaDcaDdaCcdeD45D65D75D85C9bdDe7C8acD6bD7bD8bD9bDabDbbDcbDdbC8abD04Df4CbceD16Dd5C8bdD1aCcdfD55C9bdD29D39D49D59D69D78D79D98Da7Da8Db7Db8Dc7Dd6Dd7C7acD2dD3dD4dD5dD6dD7dD8dD9dDadDbdDcdDddCacdD17D18D27D38D48D58D67D77D87Da6CacdD28D68Db6Dc6C7acD1cD3cD4cD5cD6cD7cDecCbcdD26D36D46D56D66C9bdD2aD3aD4aD5aD6aDb9Dc9Dd9C689D0cD0dC7acD1dDedC9cdD19CbceDe5C9cdD97C8bdD7aDeaC689D0fDffCbdeDa5Db5Dc5C689D0aD0bD14DfaDfbC8bdD1b"{
	batchChoice=1;
	batchProcessing ();
}

var ftabCmds = newMenu("Measurement Table Manager Menu Tool",newArray("Save a Stat Results Table","Save a Stat Results Table at the Image Location","-","Remove Last Table Line","Undo Remove Last Table Line", "-","Open a \"Stat Results Table\" Table File"));
macro "Measurement Table Manager Menu Tool - C112D01D03CcacD5cDbcDccDecDfcC35aD26D28D2aD2eD32D62Da2Dd2CfffDdeDdfDeeDefDfeDffC137D11D1aD21D41D51D52D71D81D91Db1Dc1De1Df1CeeeD00D10D20D30D40D44D47D49D4bD4dD4fD50D54D57D59D5bD5dD5fD60D70D77D79D7bD7dD7fD80D84D87D89D8bD8dD8fD90D94D97D99D9bD9dD9fDa0Db0Db4Db7Db9DbbDbdDbfDc0Dc4Dc7Dc9DcbDcdDcfDd0De0De4De7De9DebDedDf0Df4Df7Df9DfbDfdCbbbD46D48D4aD4eD56D58D5aD5eD64D69D6bD6dD6fD76D78D7aD7eD86D88D8aD8eD96D98D9aD9eDa4Da9DabDadDafDb6Db8DbaDbeDc6Dc8DcaDceDd4Dd9DdbDddDe6De8DeaDf6Df8DfaC444D63Da3Da5Dd3Dd5CcdfD45D55D74C47eD19D42C359D1cD82Dc2Df2C235D02D04D07D09D0bD0dD0fD13D23C46bD14D17D2cDb2De2C248D16D18D25C666D43D53D73D75D83D85D93D95Db3Db5Dc3Dc5De3De5Df3Df5C79cD65C777D66D68D6aD6cD6eDa6Da8DaaDacDaeDd6Dd8DdaDdcC124D06D08D0aD0cD0eD31D61Da1Dd1C46aD37C346D15D1eC58eD12D22D24D29D2bD2dD2fD72D92C579D34D39D3bD3dD3fC47dD1bD1dD1fD27C357D35D36D38D3aD3cD3eCabcD8cD9cC234D33C9abD67Da7Dd7CaccD4cD7cC013D05" {
	cmd = getArgument();
	if (cmd !="-" && cmd == "Save a Stat Results Table") {saveTab ("--",windowTabName,"");}
	if (cmd !="-" && cmd == "Save a Stat Results Table at the Image Location") {saveTab (imageFolder,windowTabName,workingima);}
	if (cmd !="-" && cmd == "Remove Last Table Line") {rmLastLine ();}
	if (cmd !="-" && cmd == "Undo Remove Last Table Line") {undormLastLine ();}
	if (cmd !="-" && cmd == "Open a \"Stat Results Table\" Table File") {openStatResultTable ();}
	if (cmd !="-" && cmd == "Open a \"Stat Results Table\" Table File") {openStatResultTable ();}
	if (cmd !="-" && cmd == "Open a \"Stat Results Table\" Table File") {openStatResultTable ();}
}

var dCmds = newMenu("Documentation & Demo Menu Tool",newArray ("Download a Pseudo Phase Contrast HUVEC Image","Download a Fluo HUVEC Image", "-", "Open On Line Documentation"));
macro "Documentation & Demo Menu Tool - C000D89D8dD96D9dDa5Da6Db3Db4Db5Db6DbdC06fD12D13D14D15D16D17D18D19D1aD1bD1cD1dD1eC0f3D32D33D34D35D36D37D38D39D3aD3bD3cD3dD3eD3fCeeeD80C444DbeDcfC73fD01Cfd0D42D43D44D45D46D47D48D49D4aD4bD4cD4dD4eCfffD7cD7fD8fD91D98Dc1DceDd1Dd3DdeDdfDe0De1De2De3De4De5De6De7De8De9DeaDebDecDedDeeDefDf0Df1Df2Df3Df4Df5Df6Df7Df8Df9DfaDfbDfcDfdDfeDffC222C0feD22D23D24D25D26D27D28D29D2aD2bD2cD2dD2eC888D88CfffD71D90DcbDd0Dd2Dd4Dd5Dd6DdcCf74D50CcccD9fCbbbD97Da0Dc2C000Da4Db2C06fD11D1fC0f4D31CeeeD70D75D81D84D92D93Dc0Dc8Dd7DdaDdbCf40D51D52D53D54D55D56D57D58D59D5aD5bD5cD5dD5eD5fCa7fD00Cfd0D41D4fC444D95DbaDc9C0ffD21D2fC4f6D30CfffD7bCf88D62D63D64D65D6bD6cCeeeD72D73D74DabDbbDccDd8Cf89D6dC666Da7C84fD02D03D04D05D06D07D08D09D0aD0bD0cD0dD0eC333D86D87D8aD99D9cDb9DbcCaaaD8eDa1Dc3Dc4Dc5Dc6CeddD76Cf89D61D66D6aD6eD6fC222DadDcdC48fD10Cfd4D40C4ffD20CaaaDacDcaDddCfaaD60C555Da2C73fD0fC333D9bDa3Da9Db1C999D7dCdddD7aD83Da8Dc7CbbbDb8C111D9aDb7CaaaD8bCf99D67D68D69C777D8cDaaCdeeD82D85C999DafCdddDd9C766D79C555DaeC999D78CcccD7eD94C888D9eC777DbfC666Db0Ca99D77"{
	cmd = getArgument();
	if (cmd!="-" && cmd == "Open On Line Documentation") {
		doc ();
	}
	if (cmd!="-" && cmd == "Download a Pseudo Phase Contrast HUVEC Image") {
		OpenImageLink(demoimagelink1,demoimagename1,1);
	}
	if (cmd!="-" && cmd == "Download a Fluo HUVEC Image") {
		OpenImageLink(demoimagelink2,demoimagename2,1);
	}
}

var dCmds = newMenu("Version and Update Information Menu Tool",newArray ("Version and Update Information", "About \"Angiogenesis Analyzer\""));
macro "Version and Update Information Menu Tool - C006DaeC3cfD91C139D8dCfedD12D40D46Dd1Df3C049D8aCeddD02D11D30D85Db9Dd0De1Df2C78aD4bCfffD48D49D55D56D65D68D69D74D9bDb4Db5DbaDc5Dc6Dc7Dc8Dc9DcaDcbDcfDd7Dd8Dd9C029DeaCdddD00D01D0eD0fD10D1fD20D95De0DefDf0Df1DfeDffC57aDadCffeD0bD47D66D75Da5Dd6C07cD97CdeeD60D90C99cDafCfffD5aD9cDabDacDbbC008D8eCdddDc4C16aDe9CffeD0cC06aD36D87CeddD21D73Da3Dc0C79cD08C03aD19D6cD6dD7bD7cD7dCddeD0aC1acDa2CfffD57D59D84D94Da4C09cDc3CeeeD05D1eD4fDdfDfbCaacD5fC007D6eDdcC7deDb1C327D8fCfeeDe2C059D77C78bD2dC048D18C29cD82CffeD1dD2eD38D3fD4aD64DfcC08cD53CeeeD2fD31D39D50DfdC9abDb8C018D4dCbeeDd2C47aD37CffeD58DdeDedC06bD98CeedD0dC6ceD23CdeeD22C2acDe4C09eD43CabcD1cC007D5eD9eDbeDcdC4dfDa1C159Da9CfedDb0C04aD89C69bDf9C039D5cC48bDe7C08bD44CeeeD13DbfDc1C9abD45D79C017D4eCcdeDf5C369D3aCfeeDb6C06aD27C8acD07C0bfD52C0adD24D62Db2Dd3CabcDb7C007D3dCadeD41C249D8cCfeeD04D78D83D93Df4C05aD17C89aD88C039D1aD3bD7aC59bDa6C08cD16CeeeDbcC9bcDb3C028DebCcdeD70D80De3C46aDddC07bDe6CeedD03Da0DaaC8bdDf6C2cfD81CcccD5bC3dfD71C149D1bD8bC04aD2aC78bDceC029D5dC48bDf8C07cD35C9acD3eC008D7eCdddD67C35aDccC7abD63C1adD15C09cD72CabcDfaC7dfD51C338D7fCfeeDeeC89bDecC049D6bC3acD92C07cD26CaabD6aC018DbdC36aDdbC7cdD14C2beDc2CbccD54C5dfD61C04aD99Da8C69bD76C039D2bD4cC48bDe8C08bDd4C018D2cCdddDdaC36aD29C06bD86D96C9abD9aC1beD33CabcDd5CbddD06C459D9fC05aDa7C99aD9dC59cDf7C08dD25C9bcD09C029D3cC56aD6fC07bDe5C7ceD32C3cfD42CdccD28C09dD34" {
	cmd = getArgument();
	if (cmd == "Version and Update Information") {
		VersionInfos ();
	}
	if (cmd == "About \"Angiogenesis Analyzer\"") {
		aboutTheTools ();
	}
}

//macro "About \"Angiogenesis Analyzer\" Action Tool - C05dD9aC678DdeDedC06dD5aCcbaD3fC06fD38D39D3aD3bD3cD4cD5cD5dDabDacDadDbaDbbDbcDbdDc8Dc9DcaDcbDccDd7Dd8Dd9DdaDdbDe6De7De8De9DeaDebC6afD22C26dDecCcceDa0Db0C05fD3dD7eD8eD9eC9abDf7C16fD6cD9cCdddD00D01D0fD10Dc0De0DefDf0Df1DfeDffC9bfD91C29fDc4CeedD02D0eD1fD20C06eD58D59DaaDb8Ca98DfcC06fD9bCbcdDd1C8bfD13Dc1C37cD3eCbdfD8cC16fD4eC9abDf5C07fD75CedcD0cC9cfD51C59eD98CefeD88C05eD68C99aD2eDbfC06eD57CdcbDdfC6bfD84C37fD17D18D19CcdeD60D70C05fD27D28D29D2aD2bD2cD4dD5eD6dD6eD9dDaeDc7DcdDdcC6aeDa7C26fD7dD8dC9cfD33D61D71C49fD95CfedD11C05fDb7DbeC99aD9fDafC18fD46CddcD05D06C8cfDb1C37fD16CddeD40CaceD87C08fDc5CeddDd0CadfD52D62C5aeD56CfeeD03C06eD67Da9C899DfbC07dD47CccbDf3C5bfD34C38eD8aCcceD50D80C9abD5fD6fD7fC8cfD41Da2Db2C29fD35Dd4CaaaD1dC06fD26D37D4bD85Dc6Dd6De5CddcD07D0dD2fC8bfD32D53D74Dc2C37eD1bCcefD54CabbDf4C17fDe4CedcDe1C9dfD92C49fDd2CeefD64C26cDceC89aDf8Df9C06fD4aCcccD0aC8bfD78C28fD25CbdfD7cC8aeD65C49fD23CeeeD30C06fD5bDb9CaaaD4fC18fD66CdedD97C8cfD31C38eD7aCdefD93CaceD79D89C08fD36Dd5CfedD55CaefD72C5afD44D45Db3CeffD83C05eD6aDa8C06dD48D49CdcbDfdC5afDb4C9cfD81Da1Da3Cba9DcfDeeCddcD08CcdfD63D94C9abDf6C58fD14C06eD6bC07eD76CdcbD1eC7bfD43Db5C7beD7bC49fDb6C37fD15CdeeD12D21C9ceD96CeedD04Da6CffeDa5C89aDfaC26eD2dCccbD0bC6bfDc3C27fD1aC58dD1cCabdDe2CedcDf2CadfD82C4afD24CeffD73C16eD99DddC17eD77CdccD09C8ceD8bC4afDd3C16eD69C99bD8fC48fDe3CaefD42CfffDa4CccdD90C17eD86" {
//	aboutTheTools ();
//}

function globalSettings () {

	archtwigSize=twigSize;
	archexcludeObjectsSize=excludeObjectsSize;
	archloopSize = loopSize;
	archseuilIsolFinal = seuilIsolFinal;
	archsizeExcludeSegments = sizeExcludeSegments;
	archnbIt = nbIt;
	archshowPassChoice = showPassChoice;
	alertMessage="";
	listCheckBoxe1=newArray ("Show maps of elements \(single analysis\)","Show nodes and junctions","Show extremities","Show loops","Show meshes","Show branches","Show segments","Analyze master tree","Show master segments","Remove small master segments","Supress isolated elements","Show supressed isolated elements");
	listCheckBoxesStatus=newArray ( lengthOf(listCheckBoxe1));
	listCheckBoxesStatus[0]=singleShowMaps;
	listCheckBoxesStatus[1]=showNodes;
	listCheckBoxesStatus[2]=showExtrema;
	listCheckBoxesStatus[3]=showLoops;
	listCheckBoxesStatus[4]=showMesh;
	listCheckBoxesStatus[5]=showTwig;
	listCheckBoxesStatus[6]=showSegments;
	listCheckBoxesStatus[7]=analyseMasterTree;
	listCheckBoxesStatus[8]=showMasterSegments;
	listCheckBoxesStatus[9]=supressSmallSegments;
	listCheckBoxesStatus[10]=supressIsolated;
	listCheckBoxesStatus[11]=showSupressIsolated;
	testOK = 0;
	while (testOK == 0) {
		Dialog.create("Angiogenis Analyzer Settings");
		Dialog.addMessage(alertMessage);
		Dialog.addCheckboxGroup(4,4, listCheckBoxe1,listCheckBoxesStatus);
		Dialog.addMessage("\n");
		Dialog.addNumber("Minimum object size", excludeObjectsSize, 0, 3, "pixels");
		Dialog.addNumber("Minimum branch size", twigSize, 0, 3, "pixels");
		Dialog.addNumber("Artifactual loop size", loopSize, 0, 3, "pixels");
		Dialog.addNumber("Isolated element size threshold", seuilIsolFinal, 0, 3, "pixels");
		Dialog.addNumber("Master segment size threshold", sizeExcludeSegments, 0, 3, "pixels");
		Dialog.addNumber("Iteration number \(advised 2 to 5\)", nbIt, 0, 3, "iteration\(s\)");
		Dialog.addNumber("Show iteration \(for single analysis\)", showPassChoice, 0, 3, "iteration\(s\)");
		html="http://image.bio.methods.free.fr/ImageJ/?Angiogenesis-Analyzer-for-ImageJ";
		Dialog.addHelp(html);
		Dialog.show();
		alertMessage="";
		singleShowMaps = Dialog.getCheckbox();
		showNodes = Dialog.getCheckbox();
		showExtrema = Dialog.getCheckbox();
		showLoops = Dialog.getCheckbox();
		showMesh = Dialog.getCheckbox();
		showTwig = Dialog.getCheckbox();
		showSegments = Dialog.getCheckbox();
		analyseMasterTree = Dialog.getCheckbox();
		showMasterSegments = Dialog.getCheckbox();
		supressSmallSegments = Dialog.getCheckbox();
		supressIsolated = Dialog.getCheckbox();
		showSupressIsolated =Dialog.getCheckbox();
		excludeObjectsSize = Dialog.getNumber();
		twigSize = Dialog.getNumber();
		loopSize = Dialog.getNumber();
		seuilIsolFinal = Dialog.getNumber();
		sizeExcludeSegments = Dialog.getNumber();
		nbIt = Dialog.getNumber();
		showPassChoice = Dialog.getNumber();
		if (supressIsolated == 0) showSupressIsolated =1;
		testOK=1;
		if (excludeObjectsSize <0 ||  isNaN(excludeObjectsSize)) {
			alertMessage = alertMessage + "! Bad settings\: " + "\"Minimum object size\" must be > 0"  + ". Replaced by the previous value.\n";
			excludeObjectsSize = archexcludeObjectsSize;
			testOK=0;
		}
		if (twigSize <0 ||  isNaN(twigSize)) {
			alertMessage = alertMessage + "! Bad settings\: " + "\"Minimum branch size\" must be \> 0"  + ". Replaced by the previous value.\n";
			twigSize = archtwigSize;
			testOK=0;
		}
		if (loopSize <0 ||  isNaN(loopSize)) {
			alertMessage = alertMessage + "! Bad settings\: " + "\"Artifactual loop size\" must be \> 0"  + ". Replaced by the previous value.\n";
			loopSize = archloopSize ;
			testOK=0;
		}
		if (seuilIsolFinal <0 ||  isNaN(seuilIsolFinal)) {
			alertMessage = alertMessage + "! Bad settings\: " + "\"Isolated element size threshold\" must be \> 0"  + ". Replaced by the previous value.\n";
			seuilIsolFinal = archseuilIsolFinal ;
			testOK=0;
		}
		if (sizeExcludeSegments <0 ||  isNaN(sizeExcludeSegments)) {
			alertMessage = alertMessage + "! Bad settings\: " + "\"Master segment size threshold\" must be \> 0"  + ". Replaced by the previous value.\n";
			sizeExcludeSegments = archsizeExcludeSegments;
			testOK=0;
		}
		if (nbIt <1 ||  isNaN(nbIt)) {
			alertMessage = alertMessage + "! Bad settings\: " + "\"Iteration number\" must be \>= 1"  + ". Replaced by the previous value.\n";
			nbIt = archseuilIsolFinal ;
			testOK=0;
		}
		if (showPassChoice > nbIt) showPassChoice = nbIt;
		if (showPassChoice < -1 || showPassChoice > nbIt ||  isNaN(showPassChoice)) {
			alertMessage = alertMessage + "! Bad settings\: " + "\"Show iteration\" must be \>= -1 and" + "<= to iteration number, currently =" + nbIt  + ". Replaced by the previous value.\n";
			showPassChoice = archshowPassChoice ;
			testOK=0;
		}
	}
}

function singleNetworkAnalysis () {
	if (is("binary")) exit ("This function requires an 8, 16 or 24 bits image.");
	getStatistics(area, mean, min, max, std, histogram);
	testDepth=pow(2,bitDepth())-1;
	if (bitDepth() == 24) testDepth = 255;
	if (mean == 0 || mean == testDepth) exit ("Empty image");
	batchChoice=0;
	workingImagTemp="";
	windowTabName=angStatsTable; nameOfStatTab="["+windowTabName+"]"; //changed stats table name to variable org name "Stat Results Table"
	nameInitImageorFolder="";
	starttime=getTime;
	imageFolder=getDirectory("image");
	findNodes (1);
	if (isOpen (workid)) selectImage (workid);
	nomdimage = getTitle;
	if (lastIndexOf(nomdimage, sufficAn) > 1) {
		workingima = substring (nomdimage,0,lastIndexOf(nomdimage, sufficAn));
	} else {workingima=nomdimage;}
	nameInitImageorFolder=workingima;
	if (isOpen (q3id)) selectImage (q3id);
	rename(workingima+"-Final Tree (q3)");
	if (isOpen (mapExtremaID)) {
		selectImage(mapExtremaID);
		run("Duplicate...", "title=Map_of_Extremities");
		singleAnaMapExtID=getImageID();
		selectImage (mapExtremaID); close ();
	}
	if (isOpen (mapNodesID)) {
		selectImage(mapNodesID);
		run("Duplicate...", "title=Map_of_Nodes");
		singleAnaMapNodeID=getImageID();
		selectImage (mapNodesID); close ();
	}
	if (singleShowMaps == 0) {
		if (isOpen (singleAnaMapExtID)) {selectImage (singleAnaMapExtID); close ();}
		if (isOpen (singleAnaMapNodeID)) {selectImage (singleAnaMapNodeID); close ();}
		if (isOpen (q3id)) {selectImage (q3id); close ();}
		if (isOpen(mapMasterTreeID)) {selectImage (mapMasterTreeID); close();}
		if (isOpen(masterTreeID)) {selectImage (masterTreeID); close();}
	}
	seconds = (getTime-starttime)/1000;
	setBatchMode("exit and display");
	if (isOpen (workid)) selectImage (workid);
	menuanalysis=0;
}

function reBuiltMap (couleur) {
	currentID=getImageID();
	selectImage (rebuiltMapId);
	run("Restore Selection");
	run("Add Selection...", "stroke=&couleur width=1 ");
	run("Select None");
	selectImage(currentID);
}

function findTree () {
	getStatistics(area, mean, min, max, std, histogram);
	testDepth=pow(2,bitDepth())-1;
	if (bitDepth() == 24) testDepth = 255;
	if (mean == 0 || mean == testDepth) exit ("Empty image");
	batchChoice=0;
	starttime=getTime;
	settings ();
	tree ();
	if (isOpen (q3id)) {
		selectImage (q3id);
		run("Invert");
		run("Select All");
		run("Copy");run("Select None");
	} else exit ("No tree found");
	if (isOpen (workid)) {
		selectImage (workid);
		if (bitDepth() == 8 || bitDepth() == 16) run("RGB Color");
		run("Duplicate...", "title=tempBlink");
		tempBlinkid=getImageID();
		setPasteMode("Add");
		run("Paste");
		run("Select None");
		setPasteMode("Copy");
	}
	seconds = (getTime-starttime)/1000;
	selectImage (q3id);
	run("Invert");
	selectImage (workid);
	if (bitDepth() == 8 || bitDepth() == 16) run("RGB Color");
	run("Add Image...", "image=tempBlink x=0 y=0 opacity=100");
	if (isOpen (tempBlinkid)) {selectImage(tempBlinkid); close();}
	selectImage (q3id);
	rename(workingima+"-Tree (q3)");
	setBatchMode("exit and display");
	if (isOpen (workid)) {selectImage(workid);}
	rename(workingima+" & Tree (overlay)");
	//print ("timing="+seconds +" sec");
}

function FindAnfRemoveLoops () {
	if (!is("binary")) exit ("This function requires a binary image.");
	getStatistics(area, mean, min, max, std, histogram);
	if (mean == 0 || mean == (pow(2,bitDepth())-1)) exit ("Empty image");
	batchChoice=0;showLoopInBlue=1;
	initLoop=getImageID();
	setBatchMode(true);
	run("Duplicate...", "title=lookloop");
	workid=getImageID();
	toAnalyse=getImageID();
	archloop=showLoops;
	showLoops=1;
	findLoop (toAnalyse,0);
	nomdimage = getTitle;
	showLoops=archloop;showLoopInBlue=0;
	selectImage(workid);
	if (getInfo("overlay") != "") {
		Overlay.copy;
		selectImage(initLoop);
		Overlay.paste
	}
	if (isOpen (initLoop)) selectImage (initLoop);
	nomdimage = getTitle;
	run("Duplicate...", "title=lookloop");
	if (lastIndexOf(nomdimage, ".") > 1) {
		workingima = substring (nomdimage,0,lastIndexOf(nomdimage, "."));
	} else {workingima=nomdimage;}
	rename(workingima+"-with Loops");
	selectImage (initLoop);
	run("Remove Overlay");
	if (isOpen (workid)) selectImage (workid);
	rename(workingima+"-with Loops Removed (q3)");
	run("Select None");
	setBatchMode("exit and display");
}

function findExtremities () {
	if (!is("binary")) exit ("This function requires a binary image.");
	getStatistics(area, mean, min, max, std, histogram);
	if (mean == 0 || mean == (pow(2,bitDepth())-1)) exit ("Empty image");
	batchChoice=0;
	initSkeleton=getImageID();
	choix=showOnlyFinal;
	showOnlyFinal=0;
	nomdimage = getTitle;
	initid=getImageID();
	imageFolder=getDirectory("image");
	sufix= ".";
	if (lastIndexOf(nomdimage, sufix) > 1) {
		workingima = substring (nomdimage,0,lastIndexOf(nomdimage, sufix));
	} else {workingima=nomdimage;}
	setBatchMode(true);
	run("Make Binary");
	selectImage(initSkeleton);
	run("Duplicate...", "title=q3Single");
	rename(workingima);
	workid=getImageID();
	run("Duplicate...", "title=q3 bis");
	q3id=getImageID();
	mapExtrema (q3id,q3id,1,1);
	selectImage (q3id);
	rename(workingima+"-Extremities (q3)");
	showOnlyFinal=choix;
	if (isOpen(workid)) {selectImage (workid); close();}
	setBatchMode("exit and display");
	selectImage (q3id);
}

function findTheNodes () {
	if (!is("binary")) exit ("This function requires a binary image.");
	getStatistics(area, mean, min, max, std, histogram);
	if (mean == 0 || mean == (pow(2,bitDepth())-1)) exit ("Empty image");
	batchChoice=0;
	initSkeleton=getImageID();
	choix=showOnlyFinal;
	showOnlyFinal=0;
	nomdimage = getTitle;
	initid=getImageID();
	imageFolder=getDirectory("image");
	sufix= ".";
	if (lastIndexOf(nomdimage, sufix) > 1) {
		workingima = substring (nomdimage,0,lastIndexOf(nomdimage, sufix));
	} else {workingima=nomdimage;}
	setBatchMode(true);
	run("Make Binary");
	selectImage(initSkeleton);
	run("Duplicate...", "title=q3Single");
	rename(workingima);
	workid=getImageID();
	run("Duplicate...", "title=q3 bis");
	q3id=getImageID();
	mapNodes (q3id,q3id,1,1,"Nodes");
	selectImage (q3id);
	rename(workingima+"-Nodes (q3)");
	showOnlyFinal=choix;
	if (isOpen(workid)) {selectImage (workid); close();}
	setBatchMode("exit and display");
	selectImage (q3id);
}

function findNodesAndBranches () {
	if (!is("binary")) exit ("This function requires a binary image.");
	getStatistics(area, mean, min, max, std, histogram);
	if (mean == 0 || mean == (pow(2,bitDepth())-1)) exit ("Empty image");
	batchChoice=0;
	mapExtremaID=0;mapNodesID=0;mapJunctionID=0; mapSegmentID=0; mapBranchID=0; mapMasterTreeID; masterTreeID; mapSegment2Id=0;
	initSkeleton=getImageID();
	choix=showOnlyFinal;
	showOnlyFinal=0;
	nomdimage = getTitle;
	initid=getImageID();
	imageFolder=getDirectory("image");
	sufix= ".";
	if (lastIndexOf(nomdimage, sufix) > 1) {
		workingima = substring (nomdimage,0,lastIndexOf(nomdimage, sufix));
	} else {workingima=nomdimage;}
	setBatchMode(true);
	run("Make Binary");
	selectImage(initSkeleton);
	run("Duplicate...", "title=q3Single");
	rename(workingima);
	workid=getImageID();
	run("Duplicate...", "title=q3 bis");
	q3id=getImageID();
	simpleNode ();
	selectImage (workid);
	rename(workingima+"-Twig Map");
	selectImage (q3id);
	rename(workingima+"-Nodes & Extremities (q3)");
	showOnlyFinal=choix;
}
//Note Here
function analyseFromATree () {
	if (!is("binary")) exit ("This function requires a binary image.");
	getStatistics(area, mean, min, max, std, histogram);
	if (mean == 0 || mean == (pow(2,bitDepth())-1)) exit ("Empty image");
	batchChoice=0;
	workingImagTemp="";
	windowTabName=angStatsTable; nameOfStatTab="["+windowTabName+"]";
	nameInitImageorFolder="";
	initSkeletonID=getImageID();
	mapExtremaID=0;mapNodesID=0;mapJunctionID=0; mapSegmentID=0; mapBranchID=0; mapMasterTreeID; masterTreeID; mapSegment2Id=0;
	nomdimage = getTitle;
	workingima=nomdimage;
	nameInitImageorFolder=workingima;
	initid=getImageID();
	imageFolder=getDirectory("image");
	imageFolderSave=imageFolder;
	sufix= ".";
	if (lastIndexOf(nomdimage, sufix) > 1) {
		workingima = substring (nomdimage,0,lastIndexOf(nomdimage, sufix));
	} else {workingima=nomdimage;}
 	run("Select None");
 	titre=workingima+sufficAn;
 	workingima=workingima+sufficAn;
	setBatchMode(true);
	run("Make Binary");
	selectImage(initSkeletonID);
	run("Duplicate...", "title=q3Single");
	workid=getImageID();
	rename(workingima);
	run("Duplicate...", "title=q3 bis");
	q3id=getImageID();
	starttime=getTime;
	findNodes (0);
	imageFolder=imageFolderSave;
	if (isOpen (mapExtremaID)) {
		selectImage(mapExtremaID);
		run("Duplicate...", "title=Map_of_Extremities");
		singleAnaMapExtID=getImageID();
		selectImage (mapExtremaID); close ();
	}
	if (isOpen (mapNodesID)) {
		selectImage(mapNodesID);
		run("Duplicate...", "title=Map_of_Nodes");
		singleAnaMapNodeID=getImageID();
		selectImage (mapNodesID); close ();
	}
	if (isOpen (initSkeletonID)) {
		selectImage (initSkeletonID);
		close ();
	}
	if (isOpen (q3id)) {
		selectImage (q3id);
		rename(workingima);
	}
	if (singleShowMaps == 0) {
		if (isOpen (singleAnaMapExtID)) {selectImage (singleAnaMapExtID); close ();}
		if (isOpen (singleAnaMapNodeID)) {selectImage (singleAnaMapNodeID); close ();}
		if (isOpen(mapMasterTreeID)) {selectImage (mapMasterTreeID); close();}
		if (isOpen(masterTreeID)) {selectImage (masterTreeID); close();}
	}
	seconds = (getTime-starttime)/1000;
	if (!isCMP) setBatchMode("exit and display");//Cancel exit mode if processed from cmp macro
	if (isOpen (workid)) selectImage (workid);
	menuanalysis=0;
}


function getMapsOFSelections () {
	setBatchMode(true);
	var segmentCount=0, branchCount=0, extremityCount=0, extremityCountEdge=0, junctionCount=0, twigCount=0, nodeCount=0, masterJunctionCount=0, masterSegmentCount=0;
	var isolatedElementCount=0, meshCount=0, isolatedTwigCount=0;
	run("Duplicate...", "title=tempOverlay");
	tempOverlayID = getImageID();
	run("To ROI Manager");
	getDimensions(width, height, channels, slices, frames);
	newImage("Re-Built Map", "8-bit White", width, height, 1);
	rebuiltMapId = getImageID();
	Array.fill(reBuiltChoicesNames,0);
	if(roiManager("count")<= 0) exit ("nothing in the ROI Manager");
	nbselections=roiManager("count");
	Dialog.create("Re-Built Map Choices");
	Dialog.addMessage("Choose the elements of the map");
	for (i=0; i< lengthOf(reBuiltElements); i++) {
		Dialog.addCheckbox(reBuiltElements[i], reBuiltChoices [i]) ;
	}
	Dialog.show();
	for (i=0; i< lengthOf(reBuiltElements); i++) {
		reBuiltChoices [i] = Dialog.getCheckbox();
		if (reBuiltChoices [i] == 1) reBuiltChoicesNames [i] = reBuiltElements[i];
	}
	selectImage (tempOverlayID);
	for (i=0; i< nbselections; i++) {
		roiManager("select",i);
		colorSelection = getInfo("selection.color");
		cathegory="";
		if (toString(colorSelection) == toString(extremityCoul)) {extremityCount++; cathegory="extremity"; name="extremity-"+extremityCount; roiManager("Rename",name);}
		if (toString(colorSelection) == toString(extremityEdgeCoul)) {extremityCountEdge++; cathegory="extremityEdge"; name="extremityEdge-"+extremityCountEdge; roiManager("Rename",name);}
		if (toString(colorSelection) == toString(junctionCoul)) {junctionCount++; cathegory="junction"; name="junction-"+junctionCount; roiManager("Rename",name);}
		if (toString(colorSelection) == toString(twigColour)) {twigCount++; cathegory="twig";name="twig-"+twigCount; roiManager("Rename",name);}
		if (toString(colorSelection) == toString(segmentCoul)) {segmentCount++; cathegory="segment";name="segment-"+segmentCount; roiManager("Rename",name);}
		if (toString(colorSelection) == toString(branchCoul)) {branchCount++; cathegory="branch";name="branch-"+branchCount; roiManager("Rename",name);}
		if (toString(colorSelection) == toString(nodesCoul)) {nodeCount++; cathegory="nodes";name="nodes-"+nodeCount; roiManager("Rename",name);}
		if (toString(colorSelection) == toString(masterJunctionCoul)) {masterJunctionCount++; cathegory="masterJunction";name="masterJunction-"+masterJunctionCount; roiManager("Rename",name);}
		if (toString(colorSelection) == toString(masterSegmentCoul)) {masterSegmentCount++; cathegory="masterSegment";name="masterSegment-"+masterSegmentCount; roiManager("Rename",name);}
		if (toString(colorSelection) == toString(isolatedElementCoul)) {isolatedElementCount++; cathegory="isolatedElement";name="isolatedElement-"+isolatedElementCount; roiManager("Rename",name);}
		if (toString(colorSelection) == toString(meshCoul)) {meshCount++;cathegory="mesh"; name="mesh-"+meshCount; roiManager("Rename",name);}
		if (toString(colorSelection) == toString(isolatedTwigCoul)) {isolatedTwigCount++; cathegory="isolatedTwig";name="isolatedTwig-"+isolatedTwigCount; roiManager("Rename",name);}
		for (j=0; j< lengthOf(reBuiltElements); j++) {
			if (cathegory == reBuiltChoicesNames[j]) { reBuiltMap (colorSelection); }
		}
		name = "";
	}
	if (isOpen (tempOverlayID)) {selectImage (tempOverlayID); close();}
	setBatchMode("exit and display");
}

//Here
function findNodes (treeOption) {
	if (batchChoice == 0 && (treeOption == 1) && (menuanalysis==0)) settings ();
	meanBgGreenValue=0;
	workid=getImageID();
	getStatistics(area, mean, min, max, std, histogram);
	testDepth=pow(2,bitDepth())-1;
	if (bitDepth() == 24) testDepth = 255;
	if (mean != 0 && mean != testDepth) {
		//step=0;
		if (treeOption == 1) tree ();
		run("Options...", "iterations=1 count=1 edm=Overwrite do=Nothing");
		node ();
		selectImage (workid);
		metaResults="Results for vascular branching analysis:\n";
		metaResults=metaResults+ "Analyzed area="+analyzedArea+ " timing=" +seconds +" sec\n";
		metaResults=metaResults+"nb Extrema="+ nbExtrema +" nb Nodes="+ nbNodes + " nb Junctions=" + nbJunctions + "\n";
		metaResults=metaResults+"Nb master segments="+ nbMasterSgment + " Tot. master segments lenght="+  totalMasterSegmentLenght + "\n";
		metaResults=metaResults+"Nb meshes="+ nbMesh + " Tot.meshes area="+  totalMeshSize + "\n";
		metaResults=metaResults+"nb total peaces="+ nbTotalPeaces +" nb segments="+ nbSegments + " nb branches="+ nbTwigs + " nb isolated segments="+ nbIsolated+"\n";
		metaResults=metaResults+"Total lenght="+totalLenght+ " total branching lenght="+totalBranchingLenght +"\n";
		metaResults=metaResults+"Total segments lenght="+totalSegmentsLenght+" total branches lenght="+ totalTwigsLenght + " total isolated branches lenght="+totalIsolatedLenght+"\n";
		metaResults=metaResults+"Branching interval=" + branchingIndex + " Mesh Index=" + meshIndex + " Mean Mesh Size = " + meanMeshSize + "\n";
		metaResults=metaResults+"Path of analyzed file=" + imageFolder ;
		selectImage (workid);
		setMetadata("Info", metaResults);
	} else {
		workingima=getTitle();
		if (lastIndexOf(workingima, ".") > 1) {
			workingima = substring (workingima,0,lastIndexOf(workingima, "."));
		}
		rename (workingima);
		run("Clear Results");
		analyzedArea=0;  nbExtrema=0; nbNodes=0; nbJunctions=0; nbMasterJunctionsFinal=0;  nbMasterSgment=0; totalMasterSegmentLenght=0; nbMesh=0; totalMeshSize=0; nbTotalPeaces=0; nbSegments=0; nbTwigs=0; nbIsolated=0; totalLenght=0; totalBranchingLenght=0; totalSegmentsLenght=0; totalTwigsLenght=0; totalIsolatedLenght=0; branchingIndex=0; meshIndex=0;  meanMeshSize=0;
	}
	if (recordStepOption == 1) {
		workingimaArch=workingima;
		workingima = workingima + "-step=" + step;
		TabWindow (1);
		workingima = workingimaArch;
		titreStep=" pass="+stepByStep+"-step"+ (step++);
		recordSteps (workid,workingima,titreStep);
	} else {
		TabWindow (1);
	}
	return workid;
}

function tree () {
	requires("1.43f");
	setBatchMode(true);
	run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel");
	depth = bitDepth; nbslice = getSliceNumber();
	if (depth != 24 && (objects == 1)) exit('Image must be 8 bit RGB encoded');
	if (bitDepth() == 24 && objects == 2) exit('Image must be 8 or 16 bits encoded');
	nomdimage = getTitle;
	initid=getImageID();
	imageFolder=getDirectory("image");
	sufix= ".";
	if (lastIndexOf(nomdimage, sufix) > 1) {
		workingima = substring (nomdimage,0,lastIndexOf(nomdimage, sufix));
	} else {workingima=nomdimage;}
 	run("Select None");
 	titre=workingima+sufficAn;
 	workingima=workingima+sufficAn;
  	run("Duplicate...", "title=&titre");
	workid = getImageID();
	if (objects == 1) q3id = treeHUVECPhase(workid); // HUVEC Phase contrast
	if (objects == 2) q3id = treeHUVECFluo(workid); // HUVEC fluo
	selectImage (q3id);
	setForegroundColor (255,255 ,255);
	getDimensions(width, height, channels, slices, frames);
	setLineWidth(1);drawRect(0, 0, width, height);
	run("Options...", "iterations=1 count=1 edm=Overwrite do=Nothing");
}

function treeHUVECPhase(workid){
	run("Duplicate...", "title=initRGB");
	imaid=getImageID();
	if (rmspeckles==1) {
		// find small bright objects
		selectImage (imaid);
		run("RGB Stack");
		setSlice(2);
		run("Duplicate...", "temp");
		tempIDee=getImageID();
		selectImage (imaid);
		run("RGB Color");
		selectImage (tempIDee);
		deflouageRelief2(tempIDee);
		histoAnalyse(1,"Mean");
		segmentSpec ();
		specID=getImageID();
		selectImage (tempIDee); close ();
		selectImage (specID);rename ("speckles");
		run("Dilate");
	}
	selectImage (imaid);
	run("RGB Stack");
	setSlice(2);
	run("Duplicate...", "tempfinale");
	run("Select All");
	run("Copy"); close ();
	selectImage (imaid);
	run("RGB Color");
	run("8-bit");
	run("Paste");
	run("Select None");
	deflouageRelief(imaid);
	run("Maximum...", "radius=0");
	run("Gaussian Blur...", "sigma=0.1");
	histoAnalyse(1,"MinError");
	if (rmspeckles==1) {
		imageCalculator("Subtract ", imaid,specID);
		if (isOpen (specID)) {selectImage  (specID); close ();}
		selectImage (imaid);
	}
	segMaskMeas=getTitle(); segMaskMeasid=getImageID();
	rename("q3"); q3Name=getTitle;
	q3id = getImageID();
	run("Options...", "iterations=2 count=3 pad edm=Overwrite do=Erode");
	run("Skeletonize");
	run("Options...", "iterations=2 count=2 pad edm=Overwrite do=Dilate");
	run("Skeletonize");
	selectImage(q3id);
	return q3id;
}

function deflouageRelief(tempDEid) {
	selectImage(tempDEid);
	run("Bandpass Filter...", "filter_large=10 filter_small=0 suppress=None tolerance=5 autoscale");
	run("Duplicate...", "title=tempdefloue");
	tempdefloueID=getImageID();
	run("Enhance Contrast", "saturated=0 normalize ");
	selectImage(tempDEid);
	run("Enhance Contrast", "saturated=0 normalize "); // a verifier
	run("Find Edges");
	imageCalculator("Average", tempDEid,tempdefloueID);
	selectImage(tempdefloueID); close ();
	selectImage (tempDEid);
}

function deflouageRelief2(tempDEid) {
	selectImage(tempDEid);
	run("Bandpass Filter...", "filter_large=10 filter_small=0 suppress=None tolerance=5 autoscale");
	run("Find Edges");
}

function treeHUVECFluo (workid) {
	selectImage(workid);
	run("Duplicate...", "title=initRGB");
	setAutoThreshold("Percentile dark");
	run("Convert to Mask");
	run("Options...", "iterations=1 count=1 pad edm=Overwrite");
	run("Dilate");
	run("Skeletonize");
	rename("q3"); q3Name=getTitle;
	q3id = getImageID();
	return q3id;
}


function segmentSpec () {
	//Speckelradius=2;
	limiteSpec=floor((Speckelradius*Speckelradius)*3);
	run("Set Measurements...", "area min perimeter shape redirect=None decimal=2");
	run("Analyze Particles...", "size=2-"+(limiteSpec)+" circularity=0-1.00 show=Masks clear");
}

function findLoop (treeID,showBigLoops) {
	if (! is("Batch Mode")) {setBatchMode(true);}
	nbLoop=1;
	nbMesh=0; totalMeshSize=0;
	if (isOpen (treeID)) selectImage (treeID); else exit ("no tree image available");
	nom="loopTemp";
	run("Duplicate...", "title=&nom");
	loopTempID = getImageID();
	run("Options...", "iterations=1 count=1 edm=Overwrite do=Dilate");
	nom="q3-loops";
	run("Duplicate...", "title=&nom");
	loopsID = getImageID();
	getDimensions(width, height, channels, slices, frames);
	selectImage (loopTempID);
	getStatistics(area, mean, min, max, std, histogram);
	loopTemp=histogram[255];
	while (loopTemp != 0) {
		nbLoop=0;
		setThreshold(255, 255, "none");
		run("Set Measurements...", "area min perimeter redirect=None decimal=2");
		run("Analyze Particles...", "size=0-Infinity circularity=0.00-1.00 show=Nothing clear include record");
		for (i=0; i< nResults; i++) {
			selectImage (loopTempID);
			setThreshold(255, 255, "red");
			doWand(getResult("XStart",i), getResult("YStart",i),0, "8-connected");
			getStatistics(area, mean, min, max, std, histogram);
			// erease objects without loop
			if (histogram[0] == 0) {
				selectImage (loopTempID);
				run("Clear", "slice");
				run("Select None");
			}
			// invert object containing almost on level of loop
			if (histogram[0] != 0) {
				run("Invert");run("Select None");
				nbLoop ++;
			}
		}
		run("Select None");
		// remove loops from loopTemp, add loops in q3-loops - loops containing onbject kept in loopTemp
		if (nbLoop > 0) {
			selectImage (loopTempID);
			run("Analyze Particles...", "size=0-Infinity circularity=0.00-1.00 show=Nothing clear include record");
			selectImage (loopTempID);
			if (nResults > 0) {
				for (i=0; i< nResults; i++) {
			 		selectImage (loopTempID);
			 		doWand(getResult("XStart",i), getResult("YStart",i),0, "8-connected");
					getStatistics(area, mean, min, max, std, histogram);
			 		if (histogram[0] == 0) {
						run("Clear", "slice");run("Select None");
						if (area < loopSize) {
							selectImage (loopsID);
							run("Restore Selection");
							run("Invert");
							run("Select None");
							if (showLoops == 1 && isOpen(workid)) {
								selectImage (workid);
								run("Restore Selection");
								run("Enlarge...", "enlarge=1");
								if (showLoopInBlue==1) {
									run("Add Selection...", "stroke=blue width=1 ");
								} else {
									run("Add Selection...", "stroke=yellow width=1 ");
								}
								run("Select None");
							}
						} else {
							if (showBigLoops == 1) {
								// detect meshes
								selectImage (workid);
								run("Restore Selection");
								nbMesh ++;
								getStatistics(area, mean, min, max, std, histogram);
								totalMeshSize = totalMeshSize + area;
								run("Enlarge...", "enlarge=-8");
								run("Add Selection...", "stroke=&meshCoul width=2 ");
								run("Select None");
							}
						}
					} else {
						run("Invert");run("Select None");
						if (showBigLoops == 1) {
							// detect meshe included in another meshe
							selectImage (workid);
							run("Restore Selection");
							nbMesh ++;
							getStatistics(area, mean, min, max, std, histogram);
							totalMeshSize = totalMeshSize + area;
							run("Enlarge...", "enlarge=-8");
							run("Add Selection...", "stroke=&meshCoul width=2 ");
							run("Select None");
						}
					}
				}
			}
		}
		selectImage (loopTempID);
		getStatistics(area, mean, min, max, std, histogram);
		loopTemp=histogram[255];
	} // end while
	run("Select None");
	selectImage (loopTempID);
	getStatistics(area, mean, min, max, std, histogram);
	loopTemp=histogram[255];
	if (isOpen (loopTempID)) {selectImage (loopTempID); close ();}
	if (isOpen (loopsID)) {
		selectImage (loopsID);
		run("Skeletonize");
		run("Select All");
		run("Copy"); close();
		selectImage (treeID);
		setPasteMode("Copy");
		run("Paste");
	}
	if (totalMeshSize >0 && nbMesh > 0) {meanMeshSize=d2s((totalMeshSize/nbMesh),1);} else {meanMeshSize=0;}
	setForegroundColor (255,255 ,255);
	getDimensions(width, height, channels, slices, frames);
	setLineWidth(1);drawRect(0, 0, width, height);
}

function simpleNode () {
	// note: q3id contains the tree
	if (isOpen (q3id)) selectImage (q3id); else exit ("no tree image available");
	getStatistics(area, mean, min, max, std, histogram);
	analyzedArea=area;
	totalLenght=0; totalSegmentsLenght=0; totalTwigsLenght=0; totalBranchingLenght=0; branchingIndex=0; meshIndex=0;
	mapExtremaID= mapExtrema (q3id,workid,1,0);
	mapNodesID= mapNodes (q3id,workid,1,0,"Nodes");
	mapTwigID=twigMap (q3id,workid,mapExtremaID,mapNodesID,0);
	mapExtremaID= mapExtrema (q3id,q3id,1,0);
	mapNodesID= mapNodes (q3id,q3id,1,0,"Nodes");
	if (isOpen(mapTwigID)) {//selectImage(mapTwigID);
			close(mapTwigID);}
	if (isOpen (mapExtremaID)) {selectImage (mapExtremaID); close ();}
	if (isOpen (mapNodesID)) {selectImage (mapNodesID); close ();}
	if (isOpen (mapJunctionID)) {selectImage (mapJunctionID); close ();}
	if (isOpen (mapSegmentID)) {selectImage (mapSegmentID); close ();}
	if (isOpen (mapBranchID)) {selectImage (mapBranchID); close ();}
}

function node () {
	// note: q3id contains the tree
	if (isOpen (q3id)) selectImage (q3id); else exit ("no tree image available");
	getStatistics(area, mean, min, max, std, histogram);
	analyzedArea=area;
	totalLenght=0; totalSegmentsLenght=0; totalTwigsLenght=0; totalBranchingLenght=0; branchingIndex=0; meshIndex=0;
	pass = -1; // initial branching
	if (showPassChoice <0) {
		mapExtremaID= mapExtrema (q3id,workid,1,0);
		mapNodesID= mapNodes (q3id,workid,1,0,"Nodes");
	}
	if (showPassChoice >=0) {
		mapExtremaID= mapExtrema (q3id,workid,0,0);
		mapNodesID= mapNodes (q3id,workid,0,0,"Nodes");
	}
	for (pass=0; pass<= nbIt; pass++) {
		if (showPassChoice >= pass) {
			if (pass <2 && findloopOption == 1) {
				findLoop (q3id,0); // remove loops created by the first main simplification
				mapExtremaID= mapExtrema (q3id,workid,0,0);
				mapNodesID= mapNodes (q3id,workid,0,0,"Nodes");
			}
			cutTwig=1;
			if (pass==showPassChoice) cutTwig=0;
			mapTwigID=twigMap(q3id,workid,mapExtremaID,mapNodesID,cutTwig);
			if (isOpen(mapTwigID)) {//selectImage (mapTwigID);
				close(mapTwigID);}
			if (isOpen (mapExtremaID)) {selectImage (mapExtremaID); close ();}
			if (isOpen (mapNodesID)) {selectImage (mapNodesID); close ();}
			if (isOpen (q3id) && smoothEachLimbing ==1) {
				// to smooth the zig-zag aspect after limbing and to avoid edge effects
				selectImage (q3id);
				run("Smooth");
				run("Maximum...", "radius=3");
				run("Make Binary");
				run("Skeletonize");
				setForegroundColor (255,255 ,255);
				getDimensions(width, height, channels, slices, frames);
				setLineWidth(8);drawRect(0, 0, width, height);
				setLineWidth(1);
				setThreshold(255, 255, "none");
			}
			//mapExtremaID= mapExtrema (q3id,workid,1,0);
			mapExtremaID= mapExtrema (q3id,workid,1,1);
			mapNodesID= mapNodes (q3id,workid,1,1,"Nodes");
			if (pass==showPassChoice && showMesh==1 && analyseMasterTree==0) {
				findLoop (q3id,1);
			}
			if (isOpen (mapJunctionID)) {selectImage (mapJunctionID); close ();}
			if (isOpen (mapSegmentID)) {selectImage (mapSegmentID); close ();}
			if (isOpen (mapBranchID)) {selectImage (mapBranchID); close ();}
		}
	}
}

function twigMap (q3twigid,destination,extremeID,nodesID,rmTwigy) {
	if (mapSegmentID ==1) {selectImage (mapSegmentID); close();}
	if ( !isOpen (extremeID) || !isOpen (nodesID) || !isOpen (q3twigid)) exit ("No maps available");
	getDimensions(width, height, channels, slices, frames);
	newImage("map Segments", "8-bit White", width, height, 1);
	mapSegmentID = getImageID();
	if (isOpen(mapBranchID)) {selectImage(mapBranchID); close();}
	getDimensions(width, height, channels, slices, frames);
	newImage("map Twigs", "8-bit White", width, height, 1);
	mapBranchID = getImageID();
	// to remove small segments
	if (isOpen (nodesID) && supressVerySmallSegments == 1) {
		// to avoide too small segments betxeen two near nodes
		selectImage (nodesID);
		nom="finalmapBigNodes";
		run("Duplicate...", "title=&nom");
		run("Options...", "iterations=1 count=1 edm=Overwrite do=Dilate");
		finalmapBigNodesID= getImageID();
		imageCalculator("Subtract create", q3twigid,finalmapBigNodesID);
		mapTwigID= getImageID(); // contains segments separated by nodes
	} else {
		imageCalculator("Subtract create", q3twigid,nodesID);
		mapTwigID= getImageID(); // contains segments separated by nodes
	}
	selectImage (mapTwigID);
	rename ("twigDetection and Analysis");
	run("Analyze Particles...", "size=0-Infinity circularity=0.00-1.00 show=Nothing clear include record");
	nbSegments=0;nbTwigs=0;nbIsolated=0;nbTotalPeaces=0;
	totalLenght=0; totalTwigsLenght=0; totalSegmentsLenght=0; totalBranchingLenght=0; totalIsolatedLenght=0;
	segmentcount=0;
	if (nResults > 0) {
		for (i=0; i< nResults; i++) {
			selectImage (mapTwigID);
			longueurSegment= floor (getResult ("Perim.",i)/2);
			doWand(getResult("XStart",i), getResult("YStart",i),0, "8-connected");
			xstartSelect=getResult("XStart",i);
			ystartSelect=getResult("YStart",i);
			selectImage (extremeID);
			run("Restore Selection");
			getStatistics(area, mean, min, max, std, histogram);
			extremHisto=histogram;
			rmCol=0;supress=0;
			if (extremHisto [255] ==1 ) {
				if (area < twigSize) {
					supress=1;
					showSelection=0;
					if (showOnlyFinal == 1 && pass == showPassChoice && showTwig==1) showSelection=1;
					if (showOnlyFinal == 0 && showTwig==1) showSelection=1;
					if (showSelection == 1 && isOpen (destination) && rmCol==0) {
						selectImage (destination);
						run("Restore Selection");
						run("Add Selection...", "stroke=&twigColour width=1"); // show too short twigs as overlay
						selectImage (mapBranchID);
						run("Restore Selection");
						run("Add Selection...", "stroke=&twigColour width=1"); // show too short twigs as overlay
						setForegroundColor(0, 0, 0);
						run("Fill", "slice");
						run("Select None");
						nbTwigs++; totalTwigsLenght=totalTwigsLenght+longueurSegment;
					}
				} else {
					if (rmCol==0)  {
						nbTwigs++; totalTwigsLenght=totalTwigsLenght+longueurSegment;
						// place to put the mean green value of twig for futur option of color sorting
					}
					run("Select None");
					showSelection=0;
					if (showOnlyFinal == 1 && pass == showPassChoice && showTwig==1) showSelection=1;
					if (showOnlyFinal == 0 && showTwig==1) showSelection=1;
					if (showSelection == 1 && isOpen (destination) && rmCol==0) {
						// futur color test here
						selectImage (destination);
						run("Restore Selection");
						run("Add Selection...", "stroke=&branchCoul width=1"); // show twigs as overlay
						selectImage (mapBranchID);
						run("Restore Selection");
						run("Add Selection...", "stroke=&branchCoul width=1"); // show nodes as overlay
						setForegroundColor(0, 0, 0);
						run("Fill", "slice");
						run("Select None");
					}
				}
				if (supress == 1 && supressTwig ==1 && rmTwigy && rmCol==0) {
					selectImage (mapTwigID);
					run("Restore Selection");
					run("Clear", "slice");
					run("Select None");
				}
			}
			if (extremHisto [255] == 2 ) {
				if (supressIsolated == 0 && rmCol==0 && longueurSegment > seuilIsolFinal) {totalIsolatedLenght=totalIsolatedLenght+longueurSegment; nbIsolated ++;}
				showSelection=0;
				if (showOnlyFinal == 1 && pass == showPassChoice && showTwig==1) showSelection=1;
				if (showOnlyFinal == 0 && showTwig==1) showSelection=1;
				if (showSelection == 1 && isOpen (destination) && rmCol==0) {
					selectImage (destination);
					run("Restore Selection");
					if (supressIsolated == 0 && longueurSegment > seuilIsolFinal) {run("Add Selection...", "stroke=&isolatedElementCoul width=1");}
					if (supressIsolated == 1 && showSupressIsolated == 1) {run("Add Selection...", "stroke=&isolatedTwigCoul width=1"); }
					if (showSupressIsolated == 1 && longueurSegment <= seuilIsolFinal) {run("Add Selection...", "stroke=&isolatedTwigCoul width=1");}
					run("Select None");
					if (supressIsolated == 1 || longueurSegment < seuilIsolFinal && rmCol==0) {
						selectImage (mapTwigID);
						run("Restore Selection");
						run("Clear", "slice");
						run("Select None");
					}
				}
			}
			if (extremHisto [255] == 0) {
				if (rmCol==0) {nbSegments ++; totalSegmentsLenght=totalSegmentsLenght+longueurSegment;}
				showSelection=0;
				if (showOnlyFinal == 1 && pass == showPassChoice && showSegments==1) showSelection=1;
				if (showOnlyFinal == 0 && showTwig==1) showSelection=1;
				if (showSelection == 1 && isOpen (destination) && rmCol==0) {
					// futur test for color here
					selectImage (destination);
					run("Restore Selection");
					run("Add Selection...", "stroke=&segmentCoul  width=1"); // show segment as overlay
					selectImage (mapSegmentID);
					run("Restore Selection");
					setForegroundColor(0, 0, 0);
					run("Fill", "slice");
					run("Add Selection...", "stroke=&segmentCoul width=1"); // show segment as overlay
					run("Select None");
					selectImage (mapSegmentID);
					run("Select None");
				}
			}
		}
		selectImage (mapTwigID);
		run("Analyze Particles...", "size=0-Infinity circularity=0.00-1.00 show=Nothing clear include record");
	}
	run("Select None");
	if (supressVerySmallSegments == 0) imageCalculator("Add ", mapTwigID, nodesID);
	if (isOpen (nodesID) && isOpen (finalmapBigNodesID) && supressVerySmallSegments == 1 ) {
		if (isOpen(finalmapBigNodesID)) {
			selectImage (finalmapBigNodesID);
			getStatistics(area, mean, min, max, std, histogram);
		}
		if (isOpen(finalmapBigNodesID) && max!=0) {
			selectImage (finalmapBigNodesID);
			run("Options...", "iterations=1 count=1 edm=Overwrite do=Dilate");
			imageCalculator("Add ", mapTwigID, finalmapBigNodesID);
			selectImage (mapTwigID);
		}
	}
	selectImage (mapTwigID);
	run("Select None");
	rename ("twig resolve");
	run("Skeletonize");
	// remove smal objects
	setForegroundColor (255,255 ,255);
	getDimensions(width, height, channels, slices, frames);
	setLineWidth(1);drawRect(0, 0, width, height);
	setThreshold(255, 255, "none");
	run("Analyze Particles...", "size=["+(excludeObjectsSize)+"] -Infinity circularity=0.00-1.00 show=Masks clear record");
	masqueid=getImageID();
	if (nResults >0) {
		rename ("mask"+pass);
		run("Select All");
		run("Copy");
		selectImage (mapTwigID);
		setPasteMode("Copy");
		run("Paste");
		run("Select None");
	}
	selectImage (masqueid); close ();
	selectImage (mapTwigID);
	run("Select All");
	run("Copy");
	selectImage (q3twigid);
	run("Paste");
	run("Select None");
	if (isOpen(finalmapBigNodesID)) {
		selectImage (finalmapBigNodesID);
		getStatistics(area, mean, min, max, std, histogram);
		if (max ==0) close ();
	}
	// option analyseMasterTree
	if (analyseMasterTree == 1 && (nbSegments >1) && pass == showPassChoice) {
		masterJunstionMap (mapSegmentID, destination);
	}
	if (isOpen(finalmapBigNodesID)) {selectImage(finalmapBigNodesID); close ();}
	if (isOpen (smootMapID)) {selectImage (smootMapID); close ();}
	// end option analyseMasterTree
	selectImage (q3twigid);
	run("Select None");
	selectImage (mapTwigID);
	run("Select None");
	nbTotalPeaces= nbIsolated+nbTwigs+nbSegments;
	totalTwigsLenght=floor(totalTwigsLenght);
	totalIsolatedLenght=floor(totalIsolatedLenght);
	totalSegmentsLenght=floor(totalSegmentsLenght);
	totalBranchingLenght=totalSegmentsLenght+totalTwigsLenght;
	totalLenght=totalBranchingLenght+totalIsolatedLenght;
	if (totalSegmentsLenght >=1 && nbTwigs >0)   { branchingIndex=d2s((totalSegmentsLenght/(nbTwigs)),3);} else {branchingIndex =0;}
	selectImage(mapTwigID);
		mapTwigTitle = getTitle();

	return mapTwigTitle;
}

function masterJunstionMap (mapSegment2Id, destination2ID) {
	setForegroundColor (255, 255, 255);
	nbMasterSgment=0; totalMasterSegmentLenght=0;
	nom="masterTreeID";
	if (isOpen(finalmapBigNodesID) && supressVerySmallSegments == 1) {
		// method without the very small segments
		selectImage (mapSegment2Id);
		run("Invert");
		imageCalculator("Add create", mapSegment2Id,finalmapBigNodesID); // verifier si n�c�ssaire
		setThreshold(255, 255, "none");
		rename(nom);
		masterTreeID=getImageID();
		run("Invert");
		run("Skeletonize");
		run("Invert");
		mapNodes (masterTreeID,destination2ID,1,0,"MasterNodes");
	} else {
		// method with the very small segments
		selectImage(nodesID);
		run("Select None");
		nom2="tempSmoothNodesMap";
		run("Duplicate...", "title=&nom2");
		smootMapID = getImageID();
		run("Options...", "iterations=1 count=3 edm=Overwrite do=Dilate");
		selectImage (mapSegment2Id);
		run("Invert");
		imageCalculator("Add create", mapSegment2Id,smootMapID);
		rename(nom);
		masterTreeID=getImageID();
		run("Invert");
		run("Skeletonize");
		run("Invert");
		mapNodes (masterTreeID,destination2ID,1,0,"MasterNodes");
	}
	// calulates and overlay the master segments
	imageCalculator("Multiply create", mapMasterJunctionID,masterTreeID);
	rename ("mapMasterTree");
	mapMasterTreeID = getImageID();
	run("Invert");
	getDimensions(width, height, channels, slices, frames);
	setForegroundColor (255, 255, 255);
	setLineWidth(1);drawRect(0, 0, width, height);
	setForegroundColor (0, 0, 0);
	run("Select None");
	run("Analyze Particles...", "size=0-Infinity circularity=0.00-1.00 show=Nothing clear include record");
	// quantify and plot master segments
	nbMasterSgment = nResults;
	if (nResults > 0) {
		for (i=0; i< nResults; i++) {
			selectImage (mapMasterTreeID);
			longueurSegment= floor (getResult ("Perim.",i)/2);
			doWand(getResult("XStart",i), getResult("YStart",i),0, "8-connected");
			getStatistics(area, mean, min, max, std, histogram);
			if (min != max) longueurSegment= floor (getResult ("Perim.",i));
			if (min == max) longueurSegment= floor (getResult ("Perim.",i)/2);
			// option to supress small master segments
			ereaseSegment=0;
			if (longueurSegment <= sizeExcludeSegments && supressSmallSegments ==1) {
				ereaseSegment=1;
			}
			xstartSelect=getResult("XStart",i);
			ystartSelect=getResult("YStart",i);
			if (ereaseSegment==0) {
				selectImage (mapSegment2Id);
				run("Restore Selection");
				if (showMasterSegments==1) run("Add Selection...", "stroke=&masterSegmentCoul width=2");
				run("Select None");
				selectImage (destination2ID);
				run("Restore Selection");
				if (showMasterSegments==1) run("Add Selection...", "stroke=&masterSegmentCoul width=2");
				run("Select None");
				selectImage (mapMasterTreeID);
				run("Select None");
				totalMasterSegmentLenght = totalMasterSegmentLenght + longueurSegment;
			} else {
				selectImage (mapSegment2Id);
				run("Restore Selection");
				run("Clear", "slice");
				run("Select None");
				selectImage (mapMasterJunctionID);
				run("Restore Selection");
				run("Enlarge...", "enlarge=2");
				run("Fill", "slice");
				run("Select None");
				nbMasterSgment=nbMasterSgment-1;
			}
		}
	}
	// plot the master junction (and quantify the merged nodes)
	selectImage (mapMasterJunctionID);
	run("Select None");
	run("Analyze Particles...", "size=0-Infinity circularity=0.00-1.00 show=Nothing clear include record");
	if (nResults > 0) {
		nbMasterJunctions=0;
		for (i=0; i< nResults; i++) {
			selectImage (mapMasterJunctionID);
			doWand(getResult("XStart",i), getResult("YStart",i),0, "8-connected");
			run("Select None");
			selectImage (destination);
			run("Restore Selection");
			run("Enlarge...", "enlarge=4");
			run("Add Selection...", "stroke=&masterJunctionCoul width=2 "); // show master junction as overlay
			run("Select None");
			nbMasterJunctions ++;
		}
	}
	if (showMesh==1 && analyseMasterTree==1) {
		selectImage (masterTreeID);
		run("Invert LUT");
		findLoop (masterTreeID,1);
	}
	nbMasterJunctionsFinal=nbMasterJunctions;
	if (totalMasterSegmentLenght > 0 && nbMasterJunctions > 0)   meshIndex=d2s((totalMasterSegmentLenght/(nbMasterJunctions)),3); else {
	meshIndex =0;}
}

// return the mean value of an histogram
function MeanStatHisto (histo,mini,maxi) {
	volhisto=0;meanhisto=0;nbPixel=0;
	if ((maxi-mini) > 0 ) {
		for (i=mini; i<= maxi; i++) {
			volhisto=volhisto +(i* histo[i]);
			nbPixel=nbPixel+(histo[i]);
		}
	}
	if (nbPixel >0) meanhisto=volhisto/nbPixel;
	return meanhisto;
}

function mapExtrema (binID,destination,overlay,keepPreviousMap) {
	//return mapExtremaID which contains the limits (as 1x1 pixels dots)
	nbExtrema=0;
	if (isOpen(mapExtremaID) && keepPreviousMap==0) {selectImage (mapExtremaID); close();} //else {exit ("no tree image available");}
	if (isOpen (binID)) {selectImage (binID);} else {exit ("no tree image available");}
	nomdArbre=getTitle;
	nom="mapExtrema";
	run("Duplicate...", "title=&nom");
	mapExtremaID = getImageID();
	roiManager("Show None");
	run("Remove Overlay");
	selectImage (mapExtremaID);
	getDimensions(width, height, channels, slices, frames);
	setBackgroundColor (255, 255, 255);
	setForegroundColor (0, 0, 0);
	run("Select All");
	run("Clear", "slice");run("Select None");
	showSelection=0;
	if (showOnlyFinal == 1 && pass == showPassChoice  && overlay==1) showSelection=1;
	if (showOnlyFinal == 0 && showExtrema==1 && overlay==1) showSelection=1;
	for (j=0; j<height; j++) {
		for (i=0; i<width; i++) {
			selectImage(binID);
			environs=0;
			snif = getPixel(i,j);
			if (snif > 0) {
				if (getPixel(i-1,j-1) == 255) environs ++;
				if (getPixel(i,j-1) == 255) environs ++;
				if (getPixel(i+1,j-1) == 255) environs ++;
				if (getPixel(i-1,j) == 255) environs ++;
				if (getPixel(i+1,j) == 255) environs ++;
				if (getPixel(i-1,j+1) == 255) environs ++;
				if (getPixel(i,j+1) == 255) environs ++;
				if (getPixel(i+1,j+1) == 255) environs ++;
				if (environs == 1)  {
					selectImage(mapExtremaID);
					setPixel(i, j, 255);
					showSelection=0;
					if (showOnlyFinal == 1 && pass == showPassChoice && showExtrema==1 && overlay==1) showSelection=1;
					if (showOnlyFinal == 0 && showExtrema==1 && overlay==1) showSelection=1;
					if (showSelection == 1 && isOpen (destination)) {
						selectImage (destination);
						makePoint(i,j);
						run("Add Selection...", "stroke=&extremityCoul width=1"); // show extremite as overlay
						makeRectangle((i-1), (j-1), 3, 3);
						run("Add Selection...", "stroke=&extremityEdgeCoul width=1"); // show extremite as overlay
					}
					nbExtrema ++;
				}
				run("Select None");
			}
		}
	}
	//mapExtremaID contains the limits
	return mapExtremaID;
}

function mapNodes (binID,destination,overlay,keepPreviousMap,calculate) {
	//returns mapNodesID which contains the map of the nodes (as 5x5 pixels circular dots)
	nbNodes=0; nbJunctions=0; nbMasterJunctions=0;
	if (isOpen(mapNodesID) && keepPreviousMap==0) {selectImage (mapNodesID); close();}
	if (isOpen (binID)) selectImage (binID); else exit ("no tree image available");
	if (isOpen (mapMasterJunctionID)) {selectImage (mapMasterJunctionID); close();}
	if (calculate=="MasterNodes") {
		getDimensions(width, height, channels, slices, frames);
		newImage("map Master Junctions", "8-bit White", width, height, 1);
		mapMasterJunctionID = getImageID();
	}
	selectImage (binID);
	nomdArbre=getTitle;
	nom="mapNodes";
	run("Duplicate...", "title=&nom");
	mapNodesID = getImageID();
	roiManager("Show None");
	run("Remove Overlay");
	nom="tempNode";
	selectImage (binID);
	run("Duplicate...", "title=&tempNode");
	tempNodeID=getImageID();
	selectImage (mapNodesID);
	getDimensions(width, height, channels, slices, frames);
	setBackgroundColor (255, 255, 255);
	setForegroundColor (0, 0, 0);
	run("Select All");
	run("Clear", "slice");run("Select None");
	j=0;
	for (j=0; j<height; j++) {
		for (i=0; i<width; i++) {
			selectImage(tempNodeID);
			environs=0;
			snif = getPixel(i,j);
			if (snif > 0) {
				if (getPixel(i-1,j-1) == 255) environs ++;
				if (getPixel(i,j-1) == 255) environs ++;
				if (getPixel(i+1,j-1) == 255) environs ++;
				if (getPixel(i-1,j) == 255) environs ++;
				if (getPixel(i+1,j) == 255) environs ++;
				if (getPixel(i-1,j+1) == 255) environs ++;
				if (getPixel(i,j+1) == 255) environs ++;
				if (getPixel(i+1,j+1) == 255) environs ++;
				if (environs > 2)  {
					selectImage(mapNodesID);
					//makeRectangle((i-1), (j-1), 3, 3);
					makeOval((i-3), (j-3), 7, 7);
					//if (pass == showPassChoice && supressVerySmallSegments == 1) makeOval((i-4), (j-4), 9, 9); //to avoid very small segments due to two near nodes (1 or 2 pixels distant)
					run("Fill", "slice");
					run("Select None");
					showSelection=0;
					if (showOnlyFinal == 1 && pass == showPassChoice && showNodes==1 && overlay==1) showSelection=1;
					if (showOnlyFinal == 0 && showNodes==1 && overlay==1) showSelection=1;
					if (showSelection == 1 && isOpen (destination)) {
						selectImage (destination);
						//makeRectangle((i-1), (j-1), 3, 3);
						makeOval((i-3), (j-3), 7, 7);
						if (calculate=="Nodes") run("Add Selection...", "stroke=&nodesCoul width=1 fill=&nodesCoul"); // show nodes as overlay
					}
					if (calculate=="Nodes") nbNodes ++;
				}
				run("Select None");
			}
		}
	}
	if (isOpen(tempNodeID)) { selectImage (tempNodeID); close ();}
	run("Select None");
	selectImage (mapNodesID);
	nom="mapJunction";
	run("Duplicate...", "title=&nom");
	mapJunctionID=getImageID();
	run("Set Measurements...", "area min perimeter redirect=None decimal=2");
	run("Analyze Particles...", "size=0-Infinity circularity=0.00-1.00 show=Nothing clear include record");
	for (i=0; i< nResults; i++) {
		selectImage (mapJunctionID);
		doWand(getResult("XStart",i), getResult("YStart",i),0, "8-connected");
		if (showSelection == 1 && isOpen (destination)) {
			selectImage (destination);
			run("Restore Selection");
			if (supressVerySmallSegments==1) {run("Enlarge...", "enlarge=0");} else {run("Enlarge...", "enlarge=1");}
			if (calculate=="Nodes") {
				run("Add Selection...", "stroke=&junctionCoul width=2 "); // show junctions as overlay
				run("Select None");
				nbJunctions ++;
			}
			if (calculate=="MasterNodes") {
				nbMasterJunctions ++;  // number to be recalculated in the final master node map
				selectImage (mapMasterJunctionID);
				run("Restore Selection");
				run("Fill", "slice");
				run("Select None");
			}
		}
		selectImage (mapJunctionID);
		run("Select None");
	}
	if (isOpen(mapJunctionID)) { selectImage (mapJunctionID); close ();}
	selectImage (mapNodesID);
	return mapNodesID;
	// mapNodes contains the nodes
}

// for segmentation of the tree
function histoAnalyse(Bright,method) {
	// get the auto threshold values
	partofHisto=0;histo1=0;
	method = method +" dark";
	setAutoThreshold(method);getThreshold (lower, upper);
	resetThreshold();
	// analyse of the histogram:
	getStatistics(area, mean, min, max, std, histogram);
	vol=newArray(3);
	vol[0]=volume (histogram,min,lower);
	vol[1]=volume (histogram,lower,upper);
	vol[2]=volume (histogram,upper,max);
	maxvol=maxOf(vol[1], vol[2]);
	for (a=1; a<3; a++) {if (maxvol == vol[a]) partofHisto = a;}
	if (partofHisto == 1) histo1=lower;
	if (partofHisto == 2) histo1=upper;
   	setThreshold(histo1, max,"black & white");
	run("Convert to Mask");
	if (Bright==0) run("Invert");
}

function volume (histo,mini,maxi) {
	volhisto=0;
	if ((maxi-mini) > 0) {
		for (i=mini; i<= maxi; i++) {
			volhisto=volhisto +(i* histo[i]);
		}
	}
	return volhisto;
}

function recordTheStepOfLimbing () {
	stepID=getImageID();
	recordStepOption=1;
	step=1;
	batchChoice=1;
	nomdimage = getTitle;
	if (lastIndexOf(nomdimage, sufficAn) > 1) {
		workingima = substring (nomdimage,0,lastIndexOf(nomdimage, sufficAn));
	} else {workingima=nomdimage;}
	nameInitImageorFolder=workingima;
	workingImagTemp=nameInitImageorFolder;
	windowTabName=angStatsTable; nameOfStatTab="["+windowTabName+"]";
	if ( !is("binary")) settings ();
	binStatus=0;
	if ( is("binary")) binStatus =1;
	imageFolder=getDirectory("image");
	imageFolderStep=imageFolder;
	for (stepByStep=-1; stepByStep<= nbIt; stepByStep++) {
		if (stepByStep < 0) {
			if ( !is("binary")) tree ();
			if (is("binary")) q3id=getImageID();
			workingImagTemp=workingima;
			totreatTemp=getImageID();
			if (isOpen(q3id)) {
			 	selectImage(q3id);
			 	titreStep=" pass="+(stepByStep)+"-step"+ (step++);
				recordSteps (q3id,workingImagTemp,titreStep);
				selectImage(q3id);
			}
			findLoop (totreatTemp,0);
			if (isOpen(q3id)) {
			 	selectImage(q3id);
			 	titreStep=" pass="+(stepByStep)+"-step"+ (step++);
				recordSteps (q3id,workingImagTemp,titreStep);
				selectImage(q3id);
			}
			batchChoice=1;
			loopIDTemp=getImageID();
			batchChoice=1;
			findNodesAndBranches ();
			batchChoice=1;
			imageFolder=imageFolderStep;
			if (isOpen(workid)) {
				selectImage(workid);
				run("Grays");
				titreStep=" pass="+(stepByStep)+"-step"+ (step++);
				recordSteps (workid,workingImagTemp,titreStep);
				selectImage(workid); close ();
			}
			if (isOpen(q3id)) {
			 	selectImage(q3id);
			 	run("Grays");
			 	titreStep=" pass="+(stepByStep)+"-step"+ (step++);
				recordSteps (q3id,workingImagTemp,titreStep);
				selectImage(q3id); close ();
			}
		} else {
			showPassChoice=stepByStep;
			if ( binStatus== 0) {
				selectImage (stepID);
				findNodes (1);
			} else {
				workingima =workingImagTemp;
			 	selectImage (stepID);
			 	rename(workingima);
				run("Duplicate...", "title=q3 bis");
				q3id=getImageID();
				findNodes (0);
			}
		}
	}
	recordStepOption=0; batchChoice=0;
}

function recordSteps (toRecordId,toRecordIma,curentStep) {

 	if (imageFolder  !="" && toRecordIma !="" )  {
 		recordFolderStep = imageFolder + toRecordIma + "=Steps" + File.separator;
 		recordFolderStepVisu = recordFolderStep + "Visus" + File.separator;
 		 File.makeDirectory(recordFolderStep);
		if (!File.exists(recordFolderStep)) exit ("Unable to create the \"" + recordFolderStep + "\"" +" directory.");
		File.makeDirectory(recordFolderStepVisu);
		if (!File.exists(recordFolderStepVisu)) exit ("Unable to create the \"" + recordFolderStepVisu + "\"" +" directory.");
		if (isOpen (toRecordId)) {
			selectImage (toRecordId);
			run("Flatten");
			pathVisu=recordFolderStepVisu + toRecordIma + curentStep + "-visu";
			saveAs("PNG",pathVisu);
			close();
			selectImage (toRecordId);
			saveAs("tiff",recordFolderStep + toRecordIma + curentStep);
		}
		if (isOpen(windowTabName)) {saveTab (recordFolderStep,windowTabName,toRecordIma);}
 	}
}

// tab for results part
// function building/managing a table window
function TabWindow (addLine) {
	// values:
	undoErease = "";
	if (! isOpen(windowTabName)) {
		run("New... ", "name="+nameOfStatTab+" type=Table");
		print(nameOfStatTab, "\\Headings:Image Name\tAnalysed area\tNb extrem.\tNb nodes\tNb Junctions\tNb master junction\tNb master segments\tTot. master segments lenght\tNb meshes\tTot.meshes area\tNb peaces\tNb segments\tNb branches\tNb isol. seg.\tTot. lenght\tTot. branching lenght\tTot. segments lenght\tTot. branches lenght\tTot. isol. branches lenght\tBranching interval\tMesh index\tMean Mesh Size\tPath");
	}
	if (addLine == 0) {print(nameOfStatTab, "\\Clear");}
	statsLine = workingima+ "\t" + analyzedArea + "\t" + nbExtrema + "\t" + nbNodes +  "\t" + nbJunctions + "\t" + nbMasterJunctionsFinal + "\t" + nbMasterSgment + "\t" + totalMasterSegmentLenght + "\t" + nbMesh + "\t" + totalMeshSize + "\t" + nbTotalPeaces +  "\t" + nbSegments+  "\t" + nbTwigs +  "\t" +nbIsolated +  "\t" + totalLenght +  "\t" + totalBranchingLenght +  "\t" + totalSegmentsLenght + "\t" + totalTwigsLenght  + "\t" + totalIsolatedLenght + "\t" + branchingIndex + "\t" + meshIndex + "\t" + meanMeshSize + "\t" + imageFolder;

	print(nameOfStatTab,  statsLine);
}

// function removing the last line of the tab
function rmLastLine () {
	if ( isOpen (windowTabName)) {
		selectWindow (windowTabName);
   		tabContent = getInfo();
   		linesInTab = split(tabContent, "\n");
		if (linesInTab[linesInTab.length-1] > 0) {
			print(nameOfStatTab, "\\Clear");
			resteLines="";
			for (i=1; i < (linesInTab.length -1); i++) {
				resteLines=resteLines+linesInTab[i] +"\n";
			}
			if (linesInTab.length > 2) print (nameOfStatTab,resteLines);
			if (linesInTab.length > 1) undoErease=linesInTab[linesInTab.length-1];
		}
	}
}

// function restoring the last ereased line in the table
function undormLastLine () {
	if (undoErease != "") print(nameOfStatTab,undoErease);
	undoErease="";
}

function openStatResultTable () {
    requires("1.39d");
  	path = File.openDialog("Select a File");
	name = File.getName(path);
	if (endsWith(name, ".xls")) {name=substring(name, 0, indexOf(name, ".xls"));} else {exit ("This file doesn't seam to be an Excel type file.");}
 	openTab (path,name);
 	imageFolder=""; path="";
}

function saveAndCloseAnalysis (closeWin) {
	workingImagTemp="";
	setBatchMode(true);
	if (isOpen (workid)) {
		selectImage (workid);
		workingImagTemp=getTitle();
	}
	if (isOpen(singleAnaMapExtID)) {
		titreStep="-Map of Extremities";
		recordSteps (singleAnaMapExtID,workingImagTemp,titreStep);
		if (closeWin==1) {selectImage (singleAnaMapExtID); close ();}
	}
	if (isOpen(singleAnaMapNodeID)) {
		titreStep="-Map of Nodes";
		recordSteps (singleAnaMapNodeID,workingImagTemp,titreStep);
		if (closeWin==1) {selectImage (singleAnaMapNodeID); close ();}
	}
	if (isOpen(masterTreeID)) {
		titreStep="-Master Tree";
		recordSteps (masterTreeID,workingImagTemp,titreStep);
		if (closeWin==1) {selectImage (masterTreeID); close ();}
	}
	if (isOpen(mapMasterTreeID)) {
		titreStep="-Map of Master Tree";
		recordSteps (mapMasterTreeID,workingImagTemp,titreStep);
		if (closeWin==1) {selectImage (mapMasterTreeID); close ();}
	}
	if (isOpen(initid)) {
		titreStep="-init";
		recordSteps (initid,workingImagTemp,titreStep);
		if (closeWin==1) {selectImage (initid); close ();}
	}
	if (isOpen(workid)) {
		titreStep="";
		recordSteps (workid,workingImagTemp,titreStep);
		if (closeWin==1) {selectImage (workid); close ();}
	}
	if (isOpen(q3id)) {
		titreStep="-Final Tree";
		selectImage (q3id);
		run("Remove Overlay");
		recordSteps (q3id,workingImagTemp,titreStep);
		if (closeWin==1) {selectImage (q3id); close ();}
	}
}

function settings () {
	Dialog.create("Settings for analysis");
	Dialog.addChoice("Kind of objects:", objectsChoices, objectsChoice);
	Dialog.show();
	objectsChoice = Dialog.getChoice();
	if (objectsChoice == objectsChoices[0]) objects = 1;
	if (objectsChoice == objectsChoices[1]) objects = 2;
}

function closeAll() {
	if (! is("Batch Mode")) {setBatchMode(true);}
	if (isOpen (mapExtremaID)) {selectImage (mapExtremaID); close ();}
	if (isOpen (mapNodesID)) {selectImage (mapNodesID); close ();}
	if (isOpen (q3id)) {selectImage (q3id); close ();}
	if (isOpen (workid)) {selectImage (workid); close ();}
	if (isOpen (initid)) {selectImage (initid); close ();}
	if (isOpen(mapMasterTreeID)) {selectImage (mapMasterTreeID); close();}
	if (isOpen(masterTreeID)) {selectImage (masterTreeID); close();}
	if (isOpen(mapMasterJunctionID)) {selectImage (mapMasterJunctionID); close();}
	if (isOpen(singleAnaMapExtID)) {selectImage (singleAnaMapExtID); close();}
	if (isOpen(singleAnaMapNodeID)) {selectImage (singleAnaMapNodeID); close();}
	if (isOpen(initSkeleton)) {selectImage (initSkeleton); close();}
}

// saving a tab as excel type file
function saveTab (path,WinTab,nameinit) {
	if (isOpen(windowTabName)) {
		if (path == "--" || path =="") {
			selectWindow (windowTabName);
			run("Input/Output...", "jpeg=75 gif=-1 file=.xls");
			saveAs("Text");
		}
		if (path != "--") {
			selectWindow(windowTabName);
			thepath = path+nameinit+"-"+WinTab+".xls";
			if (WinTab != angStatsTable) {
				thepath = path+ WinTab + ".xls";
			} else
				if (nameinit == "") {
					thepath = path+nameinit + WinTab +".xls";
				} else {
					thepath = path+nameinit+"-"+WinTab +".xls";
				}
			}
		saveAs("Text", thepath);
	} else {exit ("No Stat Results Table. Tab files have to be opened or generated by the Angiogenesis Analyzer");}
}

function openTab (path,name) {
	undoErease="";windowTabName=name; workingima="";
	lines=split(File.openAsString(path), "\n");
	if (lines.length < 2) { exit ("This file doesn't seam to contain data");}
  	headings = lines[0];
	titlesOfColumns = split(headings, ",\t");
	nameOfStatTab="["+windowTabName+"]";
	if (isOpen(windowTabName)) {selectWindow(windowTabName) ;run("Close");}
	firstLine="";
	for (i=0; i < (titlesOfColumns.length ); i++) {
		firstLine=firstLine+ titlesOfColumns [i];
		if ( i < (titlesOfColumns.length )-1) {firstLine=firstLine+ "\t";}
	}
	toPrint="";
	for (i=1; i < (lines.length ); i++) {
		toPrint=toPrint+lines[i]+"\n";
	}
	run("New... ", "name="+nameOfStatTab+" type=Table");
	print(nameOfStatTab, "\\Headings:"+firstLine+"");
	print(nameOfStatTab,toPrint);
}

//  batch processing adapter:
function batchProcessing () {
	selectedDir = getDirectory("Choose a Directory ");
	nameInitImageorFolder = File.getName(selectedDir);
	setBatchMode(true);
	settings ();
	workingImagTemp="";
	windowTabName =angStatsTable;
	nameOfStatTab="["+windowTabName+"]";	// ajuster le nom avec le nom de dossier master
	if (isOpen(windowTabName)) {print(nameOfStatTab, "\\Clear");}
	countBatch = 0; // number of file to analyse
	countFiles(selectedDir);
	incremcountBatch = 0;	fileProcessed=""; countBatchTreated=0;
	starttime = getTime();
	if (countBatch > 0) {
		showMessageWithCancel (countBatch + " Image(s) will be processed");
		if (isOpen(batchStatusWindow)) {print(batchStatus, "\\Clear");}
		processFiles(selectedDir);
		fileProcessed=fileProcessed + "Nb image(s) processed="+countBatchTreated+ ". Processing time= "+ HMS (getTime-starttime) ;
		saveTab (selectedDir,windowTabName,nameInitImageorFolder);
	} else {exit ("No .tiff file found");}
}

// count the number of file to process according to criteria, image suffix, folder name suffix, file type ?
function countFiles(selectedDir) {
	list = getFileList(selectedDir);
	for (i=0; i<list.length; i++) {
		if (endsWith(list[i], "/") && ! endsWith(list[i], TreatedFolderSuffix)) {
			countFiles(""+selectedDir+list[i]);
		} else {
			path = selectedDir+list[i];
			if (endsWith(path, ext1) || endsWith(path, ext2)) {
				countBatch++;
			}
		}
	}
}

function processFiles(selectedDir) {
	list = getFileList(selectedDir);
	for (i=0; i<list.length; i++) {
		if (endsWith(list[i], "/") && ! endsWith(list[i], TreatedFolderSuffix)){
			processFiles(""+selectedDir+list[i]);
		} else {
			showProgress(incremcountBatch++, countBatch);
			path = selectedDir+list[i];
			if (endsWith(path, ext1) || endsWith(path, ext2)) {
				countBatchTreated++;
				processFile(path,1);
				call("java.lang.System.gc"); // empty the garbage memomry
				print ("\\Clear");
				fileProcessed=fileProcessed + path+ "\n";
				// built the progess batch window
				if (!isOpen(batchStatusWindow)) {run("New... ", "name="+batchStatus+" type=Table");}
				print(batchStatus, "\\Clear");
     			print(batchStatus, fileProcessed);
				meanTime=(getTime-starttime)/countBatchTreated;
				timeRemaining=(countBatch-countBatchTreated)*meanTime;
				print(batchStatus, countBatchTreated + " image(s) processed on " + countBatch);
				print(batchStatus, "Performed processing time: " + HMS (getTime-starttime));
				print(batchStatus, "Mean processing time per image= " + HMS (meanTime));
				if ((countBatch-countBatchTreated) != 0) print(batchStatus, "Estimated processing time remaining: " + HMS (timeRemaining));
				if (isOpen(batchStatusWindow)) {
					selectWindow (batchStatusWindow);
					setLocation(0, (ecranY/2));
				}
			}
		}
	}
}

function HMS (milliSec) {
	seconds = milliSec/1000;
	Shours= (floor (seconds/3600));
	Sminutes= floor ((seconds-(3600*Shours)) /60);
	Sseconds= d2s ((seconds-(3600*Shours)-(60*Sminutes)),2);
	theTime=toString(Shours) + " h " + toString(Sminutes) + " min "+ toString(Sseconds) + " sec";
	return theTime;
}

function processFile(path,saveOption) {
	open(path);
	folderpath=getDirectory("image");
	returnedImageID=Treatment();
	if (saveOption ==1) {
		format="tiff";
		if (onFoldeByImage==1) {
			folderTreatedpath=folderpath + imageName + TreatedFolderSuffix +File.separator; // on folder for each image
		}
		if (onFoldeByImage==0) {
			foldername=File.getName(folderpath);
			folderTreatedpath=folderpath + foldername + TreatedFolderSuffix +File.separator; // on folder for all images
		}
		File.makeDirectory(folderTreatedpath);
		if (!File.exists(folderTreatedpath)) exit ("Unable to create the \"" + folderTreatedpath + "\"" +" directory.");
		selectImage(returnedImageID);nom = getTitle;
		if (folderTreatedpath != "") {
			if (isOpen (workid)) {
				selectImage (workid);
				saveAs (format, folderTreatedpath+nom);
				pathVisu=folderTreatedpath+nom+"-visu";
				run("Flatten");
				saveAs("jpeg",pathVisu);
				close();
			}
		}
	}
	closeAll (); // for residual analysis images
}

// Image treatment function
function Treatment() {
	analysedID=findNodes (1);
	return analysedID;
}

macro "Hide overlay [h]" {
	run("Hide Overlay");
}

macro "Show overlay [s]" {
	run("Show Overlay");
}

// for blinking to evaluate the pertinancy
macro "Blink overlay [b]" {
	blinkOverlay () ;
}

function blinkOverlay () {
	source=getImageID();
 	blink (source);
}

function blink (image) {
	selectImage (image);
	start = getTime;
	while (click()==0 && (getTime-start)<4000) {
		wait(500);click();run("Hide Overlay");
		wait(500);run("Show Overlay");
	}
}

// adapted from the GetCursorLocDemo macro available at the
// http://rsb.info.nih.gov/ij/macros/GetCursorLocDemo.txt
function click() {
	rightButton=4;leftButton=16;
	getCursorLoc(x, y, z, flags);
	if (flags&leftButton!=0) exit;
    if (flags&rightButton!=0) exit;
    return 0;
}

// additional utilities
function aboutTheTools () {
	requires("1.46i");
	about="--------------------------- \"Angiogenesis Analyzer\" ------------------------------\n\n";
	about= about+"This set of tools allows analysis of cellular networks. Typically it can detect and analyze the\n";
	about= about+"pseudo vascular organization of endothelial cells cultured in gel medium.\n";
	about= about+"Version 1.0.a for numerical support of 2012 edition of the  ImageJ Conferences, Luxembourg.\n";
	about= about+"\n";
	about= about+"------------------------------------------------------------------------------\n";
	about= about+"Installation: the tools file has to be stored in the \"ImageJ/macros/toolset\" repertory \n";
	about= about+"------------------------------------------------------------------------------\n";
	about= about+"Short documentation:\n\n";
	about= about+"- \" Network Analysis Menu Tool\" provides functions to analyze phase contrast and fluorescent\n";
	about= about+"   images of endothelial cell networks.\n";
	about= about+"- \" Blurred Mask Tool\" removes the gradient in a selected area, useful to remove aggregates of\n";
	about= about+"   dirt in phase contrast images.\n";
	about= about+"- \" Tuning Functions Menu Tool\" regroups a series of functions to test each step of analysis \n";
	about= about+"   the best settings for a new image origin. It also provides a tool to get \n";
	about= about+"   customized overlays of detected elements.\n";
	about= about+"- \" Batch Image Treatment Tool\" allows analysis of batch of images in several level of directories.\n";
	about= about+"   The function records individual results summarized in a table as an Excel like file and\n";
	about= about+"   analyzed images containing the detected elements as overlay.\n";
	about= about+"- \" Measurement Documents Menu Tool\" regroups functions to open, edit and save the result\n";
	about= about+"- \" Online Documentation and Demo\" tool bar menu gives some internet resources including\n";
	about= about+"  online documentation and downloadable images samples.\n";
	about= about+"- \" Version and Update Information\" tool bar menu provides version and update information.\n";
	about=about + "\n------------------------------------------------------------------------------";
	about=about +"\nAuthor : Gilles Carpentier"+"\nFaculte des Sciences et Technologie"+"\nUniversite Paris Est Creteil Val de Marne, France.";
	about=about + "\n------------------------------------------------------------------------------\n";
	// from PrintToTextWindow macro available at the http://rsbweb.nih.gov/ij/macros/PrintToTextWindow.txt
	// author: Wayne Rasband
	title1 = "Infos for the \"Angiogenesis Analyzer\"";
	title2 = "["+title1+"]";
	f = title2;
	if (isOpen(title1)) {
		print(f, "\\Update:"); // clears the window
  		print(f, about);
		selectWindow (title1);
	} else {
		run("New... ", "name="+title2+" type=[Text File] width=80 height=16");
  		print(f, about);
	}
}

function netTest () {
	if (indexOf (File.openUrlAsString(urllist), errorNetMessage) >0) exit("You need an internet access to run this function.");
}

function doc () {
	netTest ();
	showMessageWithCancel  ("A notice is avaible on line. Open it with your default web browser?");
	run("URL...", "url=["+onlinedoclink +"]");
}

function OpenImageLink(link,name,question) {
	// Check if already downloaded.
	demoimalocation = getDirectory("startup");
	fildestination = demoimalocation+ "Downloaded Demo Images/" + name;
	if (File.exists(fildestination)) {
		if (question ==1 ) showMessageWithCancel ("The \"" + name + "\" has already been downloaded. Open it?");
		open(fildestination);
	}
	else {
		netTest ();
		showMessageWithCancel ("ImageJ will download a demo image. Continue?");
		run("URL...", "url=["+link+"]");
		imageid = getImageID();
		nomdimage = getTitle;
		// Create a <Downloaded Demo Images> repertory in ImageJ folder.
		ImaDemo = demoimalocation+"Downloaded Demo Images"+File.separator;
		File.makeDirectory(ImaDemo);
		if (!File.exists(ImaDemo)) exit("Unable to create directory, something wrong in the ImageJ folder");
		selectWindow(nomdimage);
		save(""+ImaDemo+""+ nomdimage +"");
	}
}

////// end addition utilities
// -------------------*** Additionnal code for on line update resources ***-----------------------------

//Developer info
//Kind:Toolset
//Title:"Angiogenesis Analyzer"
//Version:1.0.c
//Date: 03 Dec 2013
//Origin:NIH
//NotUpdateThisFile
//End info

function VersionInfos () {
	// variables for on line update resources
	beginsign="//Developer info";endsign="//End info";
	kind="toolsets/";
	urlrep="http://image.bio.methods.free.fr/ij/ijmacro/Angiogenesis/";
	name="Angiogenesis Analyzer.txt";
	namedev="Angiogenesis Analyzer.txt";
	favoritefoldername= "Macros";
	version=versionMessage();
	if (indexOf(version, "install it?" ) > 0 ) {
		macrotext=getdistantmacro (namedev,urlrep);macrolocal="";
		macropath=getDirectory("macros")+kind+namedev;
		if (File.exists(macropath)) {macrolocal=File.openAsString(macropath);}
		if (macrotext != macrolocal) {
			//perfom the installation
			Dialog.create("New version installation option");
			Dialog.addMessage(version);
			Dialog.addCheckbox("Install a Plugin Shortcut?", 0);
			Dialog.addMessage("(This option provides a shortcut in the plugins menu of ImageJ, making easier\nthe next use of the new installed version).");
			Dialog.show();
			plugin= Dialog.getCheckbox();
			f= File.open(macropath);
			print (f,macrotext);
			File.close(f);
			if (plugin ==1) {InstallPluginsStarter(namedev);}
			message="The installation of the "+giveDevInfo (macrotext,1)+ " "+ giveDevInfo (macrotext,2)+ "is completed.";
			message=message+ " Do you want to run it?";
			showMessageWithCancel(message);
			run("Install...", "install=["+macropath+"]");
		}
	} else {showMessage (version);}// comment without installation available
}

function versionMessage() {
	version="";
	if (getDirectory("startup") == 0) exit ("Unable to find the startup directory, something wrong in the ImageJ folder");
	if (getDirectory("macros") == 0) exit ("Unable to find the macros directory, something wrong in the ImageJ folder");
	MacroPath=getDirectory("macros");thismacropath=MacroPath+kind+name;
	if (! File.exists(thismacropath)) exit ("This macro has to be recorded under the name of \"" +name+"\"\ninto the \"macros/"+kind+"\" folder of ImageJ.");
	macrotext=File.openAsString(thismacropath);
	macrotextdistant=getdistantmacro (namedev,urlrep);
	version="";macrolocal="";
	version=version + "\n \nThis version of the " + giveDevInfo (macrotext,1) + " " + giveDevInfo (macrotext,2);
	version=version + "is provided by the " + giveDevInfo (macrotext,5)+ " web site.";
	version=version + "\nVersion number: " + giveDevInfo (macrotext,3)+ " - " + giveDevInfo (macrotext,4) +".";
	if (macrotextdistant !="" ) {
		new=giveDevInfo (macrotextdistant,3);old=giveDevInfo (macrotext,3);
		if (new > old) {
			macropath=getDirectory("macros")+kind+namedev;
			if (File.exists(macropath)) {macrolocal=File.openAsString(macropath);}
			if (macrotextdistant != macrolocal) {
				update="\n \nA new version "+new+ " is available on the "  +giveDevInfo (macrotextdistant,5)+ " web site: ";
				update=update+ "\n \nDo you want to install it?";
			} else {
				update ="\n \nThe latest "+new+" version called \"" +namedev+ "\" provided by \nthe "+giveDevInfo (macrotextdistant,5) +" web site has already be installed";
				update = update+ " in the \"" +kind+ "\" repertory \nof ImageJ.";
			}
		} else {
			update="No new version available.";
		}
		version=version +"\n" + update ;
	}
	return version;
}

function giveDevInfo (text,n) {
	lines=split(text,"\n");
	if ( (indexOf(text, beginsign)<0) || (indexOf(text, endsign)<0) ) exit ("Not upgradable macro code.");
	for (i=0; lines[i] != endsign; i ++) {}
	for (j=i; lines[j] != beginsign; j --) {}
	infotext=newArray(i-j-1);
	for (i=0; i < infotext.length; i ++) {infotext[i]=lines[i+j+1];}
	info=infotext[n-1]; signature=":";
	cut = indexOf(info, signature);
	info = substring(info,(cut+lengthOf(signature)),lengthOf(info));
	return info;
}

// Function giving the content of a distant macro (name) located at the distant repertory (urlrep).
function getdistantmacro (name,urlrep) {
	macrotextnih="";
	erNetMessage ="Error: ";
	testlink = "http://rsb.info.nih.gov/ij/macros/Arrays.txt";
	if (indexOf (File.openUrlAsString(testlink), erNetMessage) < 0) {
		distantmacrolink = urlrep + name;
		if (indexOf(distantmacrolink, " ") > -1) {
			while (indexOf(distantmacrolink, " ") > -1) {
				distantmacrolink=substring(distantmacrolink, 0, (indexOf(distantmacrolink, " ")))+"%20"+substring(distantmacrolink, (indexOf(distantmacrolink, " ")+1),lengthOf(distantmacrolink) );
			}
		}
		showStatus("Internet link...");
		macrotextnih =File.openUrlAsString(distantmacrolink);
		showStatus("");
	} else { showMessage ("No internet connection to looks for new version.");}
	return macrotextnih;
}

function InstallPluginsStarter(macroname) {
	// from MacroPluginShortcutsTool.txt
	codestarter = "run\(\"Install...\", \"install=[\"+getDirectory(\"macros\")+\""+kind+ macroname + "\]\"\);";
	if (getDirectory("plugins") == "") exit ("Unable to find the Plugins directory; something wrong in the ImageJ folder.");
	if (endsWith(macroname, ".txt") || endsWith(macroname, ".ijm")) pluginname = substring(macroname, 0, (lengthOf(macroname)-4));
	StarterDir = getDirectory("plugins")+favoritefoldername+File.separator;
	File.makeDirectory(StarterDir);
	if (!File.exists(StarterDir)) exit ("Unable to create "+favoritefoldername+" Macros directory, something wrong in the ImageJ folder.");
	starterplugin = StarterDir + pluginname +"_ .ijm";
	f= File.open(StarterDir + pluginname +"_ .ijm");
	print (f,codestarter);
	File.close(f);
	showMessage ("The plugin shortcut \"" +pluginname+ "\" will be available after\nImageJ restarting, in the \"Plugins->" + favoritefoldername + "\" menu.");
}

// *** End of additionnal code for on line update ressources ***
//sectionEnd>
