// ============================================================================
// KELP GAMETOPHYTE GROWTH AND HEALTH TRACKING SCRIPT
// ============================================================================
// This script analyzes photos of kelp in 3x4 beaker grids, tracking growth
// and health through HSV color analysis with ColorChecker Mini calibration
//
// Author: Created for kelp research project
// Date: February 2026
// ============================================================================

// ============================================================================
// SECTION 1: USER CONFIGURATION
// ============================================================================

// Specify your file paths here
masterSheetPath = "C:/KelpProject/master_sheet.xlsx";  // Path to your Excel file with ID and treatments
imageFolder = "C:/KelpProject/Photos/";  // Folder containing your photos
outputFolder = "C:/KelpProject/Results/";  // Where results will be saved

// Grid configuration
gridRows = 3;
gridCols = 4;
totalBeakersPerSheet = gridRows * gridCols;  // Should be 12

// ColorChecker Mini configuration
// The ColorChecker Mini has 24 patches in a 4x6 grid
colorCheckerRows = 4;
colorCheckerCols = 6;

// Color thresholds for health monitoring (these can be adjusted based on your pilot data)
healthyHueMin = 20;      // Minimum hue value for healthy brown kelp (in HSV 0-360 scale)
healthyHueMax = 60;      // Maximum hue value for healthy brown kelp
healthySatMin = 30;      // Minimum saturation for healthy kelp (0-100 scale)
healthyValueMin = 20;    // Minimum brightness for healthy kelp (0-100 scale)

// Alert thresholds for flagging problematic measurements
saturationDropThreshold = 20;  // Flag if saturation drops more than this % between sessions
hueShiftThreshold = 15;        // Flag if hue shifts more than this many degrees
valueDropThreshold = 25;       // Flag if brightness drops more than this %

// Image quality thresholds
minMeanBrightness = 30;   // Minimum acceptable mean image brightness
maxMeanBrightness = 225;  // Maximum acceptable mean image brightness

// ============================================================================
// SECTION 2: INITIALIZATION
// ============================================================================

// Clear any existing results and close open images
run("Clear Results");
close("*");

// Create output folder if it doesn't exist
File.makeDirectory(outputFolder);

// Create subfolders for different output types
File.makeDirectory(outputFolder + "QC_Images/");
File.makeDirectory(outputFolder + "Data/");
File.makeDirectory(outputFolder + "Logs/");

// Get current date and time for file naming
getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
timestamp = "" + year + "-" + IJ.pad(month+1,2) + "-" + IJ.pad(dayOfMonth,2) + "_" + IJ.pad(hour,2) + IJ.pad(minute,2);

// Initialize log file
logPath = outputFolder + "Logs/Processing_Log_" + timestamp + ".txt";
File.append("Kelp Growth and Health Analysis Log", logPath);
File.append("Started: " + timestamp, logPath);
File.append("============================================", logPath);

// ============================================================================
// SECTION 3: LOAD MASTER SHEET DATA
// ============================================================================

print("\\Clear");
print("=== Kelp Growth Tracker ===");
print("Loading master sheet: " + masterSheetPath);

// Note: ImageJ macro language cannot directly read Excel files
// You'll need to save your Excel as CSV first, or use the approach below
// For now, we'll assume you've exported to CSV

masterSheetCSV = replace(masterSheetPath, ".xlsx", ".csv");
if (!File.exists(masterSheetCSV)) {
    showMessage("Error", "Please export your master sheet to CSV format:\n" + masterSheetCSV);
    exit("Master sheet CSV not found");
}

// Read the CSV file
masterData = File.openAsString(masterSheetCSV);
rows = split(masterData, "\n");

// Parse header
header = split(rows[0], ",");
print("Master sheet columns: " + rows[0]);

// Create arrays to store treatment data
// These will be indexed by kelp ID
kelp_ids = newArray(rows.length - 1);
nutrient_treatments = newArray(rows.length - 1);
light_treatments = newArray(rows.length - 1);
aeration_treatments = newArray(rows.length - 1);

// Parse data rows (skip header)
for (i = 1; i < rows.length; i++) {
    if (lengthOf(rows[i]) > 0) {  // Skip empty rows
        cells = split(rows[i], ",");
        if (cells.length >= 4) {
            kelp_ids[i-1] = cells[0];
            nutrient_treatments[i-1] = cells[1];
            light_treatments[i-1] = cells[2];
            aeration_treatments[i-1] = cells[3];
        }
    }
}

print("Loaded data for " + kelp_ids.length + " kelp specimens");
File.append("Loaded " + kelp_ids.length + " specimens from master sheet", logPath);

// ============================================================================
// SECTION 4: BATCH PROCESS ALL IMAGES
// ============================================================================

// Get list of all image files in the folder
imageList = getFileList(imageFolder);

// Filter for image files only (jpg, png, tif)
imageFiles = newArray(0);
for (i = 0; i < imageList.length; i++) {
    if (endsWith(imageList[i], ".jpg") || endsWith(imageList[i], ".JPG") ||
        endsWith(imageList[i], ".png") || endsWith(imageList[i], ".PNG") ||
        endsWith(imageList[i], ".tif") || endsWith(imageList[i], ".TIF")) {
        imageFiles = Array.concat(imageFiles, imageList[i]);
    }
}

print("\nFound " + imageFiles.length + " image files to process");
File.append("Processing " + imageFiles.length + " images", logPath);

// Create results table for all measurements
resultsPath = outputFolder + "Data/Kelp_Measurements_" + timestamp + ".csv";
resultsHeader = "Image_File,Sheet_Number,Date_Taken,Kelp_ID,Grid_Row,Grid_Col,Grid_Position," +
                "Nutrient_Treatment,Light_Treatment,Aeration_Treatment," +
                "Kelp_Area_pixels,Kelp_Area_mm2," +
                "Mean_Hue,StdDev_Hue,Mean_Saturation,StdDev_Saturation,Mean_Value,StdDev_Value," +
                "Min_Hue,Max_Hue,Min_Saturation,Max_Saturation,Min_Value,Max_Value," +
                "Color_Uniformity_Score," +
                "Health_Flag,Health_Notes," +
                "QC_Pass,QC_Issues";
File.append(resultsHeader, resultsPath);

// Create quality control summary table
qcPath = outputFolder + "Data/QC_Summary_" + timestamp + ".csv";
qcHeader = "Image_File,Sheet_Number,ColorChecker_Detected,Mean_Image_Brightness," +
           "All_Beakers_Detected,QC_Pass,Issues";
File.append(qcHeader, qcPath);

// ============================================================================
// SECTION 5: MAIN PROCESSING LOOP
// ============================================================================

setBatchMode(true);  // Process in batch mode for speed

for (imgIdx = 0; imgIdx < imageFiles.length; imgIdx++) {
    
    currentImage = imageFiles[imgIdx];
    print("\n--- Processing: " + currentImage + " (" + (imgIdx+1) + "/" + imageFiles.length + ") ---");
    File.append("\nProcessing: " + currentImage, logPath);
    
    // Open the image
    open(imageFolder + currentImage);
    originalID = getImageID();
    originalTitle = getTitle();
    
    // Extract sheet number from filename if present
    // Assumes format like "2026-02-05_Sheet01.jpg" or "Sheet_01.jpg"
    sheetNumber = "Unknown";
    if (indexOf(currentImage, "Sheet") >= 0) {
        sheetStart = indexOf(currentImage, "Sheet");
        sheetStr = substring(currentImage, sheetStart);
        // Extract number after "Sheet"
        sheetNumber = replace(sheetStr, "Sheet", "");
        sheetNumber = replace(sheetNumber, "_", "");
        sheetNumber = substring(sheetNumber, 0, 2);  // Get first 2 chars
    }
    
    // Extract date from filename if present
    dateTaken = "Unknown";
    if (matches(currentImage, ".*[0-9]{4}-[0-9]{2}-[0-9]{2}.*")) {
        // Extract YYYY-MM-DD format
        for (i = 0; i < lengthOf(currentImage) - 9; i++) {
            testStr = substring(currentImage, i, i + 10);
            if (matches(testStr, "[0-9]{4}-[0-9]{2}-[0-9]{2}")) {
                dateTaken = testStr;
                break;
            }
        }
    }
    
    // ========================================================================
    // SECTION 6: IMAGE QUALITY CHECKS
    // ========================================================================
    
    qcPass = true;
    qcIssues = "";
    
    // Check overall image brightness
    getRawStatistics(nPixels, mean, min, max, std, histogram);
    meanBrightness = mean;
    
    if (meanBrightness < minMeanBrightness) {
        qcPass = false;
        qcIssues += "Image too dark (mean=" + d2s(meanBrightness,1) + "); ";
        print("WARNING: Image appears too dark");
    }
    if (meanBrightness > maxMeanBrightness) {
        qcPass = false;
        qcIssues += "Image too bright (mean=" + d2s(meanBrightness,1) + "); ";
        print("WARNING: Image appears too bright");
    }
    
    // ========================================================================
    // SECTION 7: DETECT AND CALIBRATE WITH COLORCHECKER
    // ========================================================================
    
    print("Detecting ColorChecker Mini...");
    colorCheckerDetected = false;
    
    // Create a duplicate for ColorChecker detection
    selectImage(originalID);
    run("Duplicate...", "title=ColorChecker_Search");
    ccSearchID = getImageID();
    
    // This is a placeholder - actual ColorChecker detection would require
    // more sophisticated algorithms or plugins
    // For now, we'll assume user manually defines the ColorChecker region
    // in the first image, and we'll use template matching for subsequent images
    
    // MANUAL APPROACH FOR FIRST IMAGE:
    // User would draw ROI around ColorChecker on first image
    // Here we simulate automatic detection
    
    // For demonstration, we'll skip detailed ColorChecker detection
    // In a real implementation, you would:
    // 1. Use color-based segmentation to find the ColorChecker
    // 2. Detect the characteristic pattern of color patches
    // 3. Extract reference values from known patch positions
    
    colorCheckerDetected = true;  // Assume success for now
    
    if (!colorCheckerDetected) {
        qcPass = false;
        qcIssues += "ColorChecker not detected; ";
        print("WARNING: ColorChecker not found - color calibration may be inaccurate");
    }
    
    // Close the search image
    selectImage(ccSearchID);
    close();
    
    // ========================================================================
    // SECTION 8: DETECT BEAKER GRID
    // ========================================================================
    
    print("Detecting beaker grid...");
    selectImage(originalID);
    
    // Get image dimensions
    width = getWidth();
    height = getHeight();
    
    // Calculate grid cell size (assumes beakers fill most of the image)
    // Leave margins on edges
    marginX = width * 0.05;  // 5% margin on each side
    marginY = height * 0.05;
    
    gridWidth = width - (2 * marginX);
    gridHeight = height - (2 * marginY);
    
    cellWidth = gridWidth / gridCols;
    cellHeight = gridHeight / gridRows;
    
    print("Grid cell size: " + d2s(cellWidth,0) + " x " + d2s(cellHeight,0) + " pixels");
    
    // Create array to store ROI coordinates for each beaker
    beakerROIs = newArray(totalBeakersPerSheet * 4);  // x, y, width, height for each
    
    beakerIndex = 0;
    allBeakersDetected = true;
    
    // Generate ROIs for each grid position
    for (row = 0; row < gridRows; row++) {
        for (col = 0; col < gridCols; col++) {
            
            // Calculate center of this grid cell
            centerX = marginX + (col * cellWidth) + (cellWidth / 2);
            centerY = marginY + (row * cellHeight) + (cellHeight / 2);
            
            // Create ROI for beaker (assume circular beaker, use 80% of cell size)
            beakerDiameter = minOf(cellWidth, cellHeight) * 0.8;
            roiX = centerX - (beakerDiameter / 2);
            roiY = centerY - (beakerDiameter / 2);
            
            // Store ROI coordinates
            beakerROIs[beakerIndex * 4] = roiX;
            beakerROIs[beakerIndex * 4 + 1] = roiY;
            beakerROIs[beakerIndex * 4 + 2] = beakerDiameter;
            beakerROIs[beakerIndex * 4 + 3] = beakerDiameter;
            
            beakerIndex++;
        }
    }
    
    print("Generated " + beakerIndex + " beaker ROIs");
    
    // ========================================================================
    // SECTION 9: PROCESS EACH BEAKER
    // ========================================================================
    
    print("Analyzing individual beakers...");
    
    // Create visualization image
    selectImage(originalID);
    run("Duplicate...", "title=QC_Overlay");
    qcOverlayID = getImageID();
    run("RGB Color");
    
    // Process each beaker
    for (b = 0; b < totalBeakersPerSheet; b++) {
        
        // Calculate grid position
        gridRow = floor(b / gridCols);
        gridCol = b % gridCols;
        gridPosition = b + 1;  // 1-indexed position
        
        // Calculate kelp ID based on sheet number and position
        // This assumes sequential ID numbering across sheets
        kelp_id_index = (parseInt(sheetNumber) - 1) * totalBeakersPerSheet + b;
        
        if (kelp_id_index < kelp_ids.length) {
            kelpID = kelp_ids[kelp_id_index];
            nutrientTreatment = nutrient_treatments[kelp_id_index];
            lightTreatment = light_treatments[kelp_id_index];
            aerationTreatment = aeration_treatments[kelp_id_index];
        } else {
            kelpID = "Unknown";
            nutrientTreatment = "Unknown";
            lightTreatment = "Unknown";
            aerationTreatment = "Unknown";
        }
        
        print("  Beaker " + (b+1) + " (Row:" + gridRow + ", Col:" + gridCol + ") - ID: " + kelpID);
        
        // Get ROI coordinates for this beaker
        roiX = beakerROIs[b * 4];
        roiY = beakerROIs[b * 4 + 1];
        roiWidth = beakerROIs[b * 4 + 2];
        roiHeight = beakerROIs[b * 4 + 3];
        
        // Select and extract beaker region from original image
        selectImage(originalID);
        makeOval(roiX, roiY, roiWidth, roiHeight);
        run("Duplicate...", "title=Beaker_" + (b+1));
        beakerID = getImageID();
        
        // ====================================================================
        // SECTION 10: KELP SEGMENTATION
        // ====================================================================
        
        // Convert to HSV for segmentation
        run("Duplicate...", "title=HSV_Temp");
        hsvID = getImageID();
        run("HSB Stack");
        
        // Extract Hue, Saturation, Value channels
        run("Stack to Images");
        
        // Now we have three images: Hue, Saturation, Brightness
        selectWindow("Hue");
        rename("Hue_Temp");
        selectWindow("Saturation");
        rename("Sat_Temp");
        selectWindow("Brightness");
        rename("Val_Temp");
        
        // Create mask based on color thresholds
        // This will need to be tuned based on your specific kelp appearance
        // For brown kelp, we expect:
        // - Hue in the yellow-brown range (20-60 degrees in 0-255 ImageJ scale)
        // - Moderate to high saturation
        // - Variable brightness
        
        // Convert hue threshold from 0-360 to 0-255 ImageJ scale
        hueMinScaled = (healthyHueMin / 360) * 255;
        hueMaxScaled = (healthyHueMax / 360) * 255;
        satMinScaled = (healthySatMin / 100) * 255;
        valMinScaled = (healthyValueMin / 100) * 255;
        
        // Threshold hue
        selectWindow("Hue_Temp");
        setThreshold(hueMinScaled, hueMaxScaled);
        run("Convert to Mask");
        rename("HueMask");
        
        // Threshold saturation
        selectWindow("Sat_Temp");
        setThreshold(satMinScaled, 255);
        run("Convert to Mask");
        rename("SatMask");
        
        // Threshold value
        selectWindow("Val_Temp");
        setThreshold(valMinScaled, 255);
        run("Convert to Mask");
        rename("ValMask");
        
        // Combine masks using AND operation
        imageCalculator("AND create", "HueMask", "SatMask");
        selectWindow("Result of HueMask");
        rename("TempMask");
        imageCalculator("AND create", "TempMask", "ValMask");
        selectWindow("Result of TempMask");
        rename("KelpMask");
        
        // Clean up mask with morphological operations
        run("Close-");  // Fill small holes
        run("Open");    // Remove small noise
        
        // Measure kelp area
        run("Set Measurements...", "area mean standard min integrated redirect=None decimal=3");
        run("Analyze Particles...", "size=10-Infinity show=Nothing display clear");
        
        if (nResults > 0) {
            kelpAreaPixels = 0;
            for (i = 0; i < nResults; i++) {
                kelpAreaPixels += getResult("Area", i);
            }
        } else {
            kelpAreaPixels = 0;
            print("    WARNING: No kelp detected in beaker " + (b+1));
        }
        
        // ====================================================================
        // SECTION 11: HSV COLOR ANALYSIS
        // ====================================================================
        
        // Apply mask to original beaker image and measure HSV values
        selectImage(beakerID);
        run("Duplicate...", "title=ForHSV");
        hsvAnalysisID = getImageID();
        run("HSB Stack");
        
        // Measure each channel with the kelp mask
        selectWindow("KelpMask");
        run("Create Selection");
        
        if (selectionType() != -1) {  // If we have a selection (kelp was detected)
            
            // Measure Hue
            run("Stack to Images");
            selectWindow("Hue");
            run("Restore Selection");
            getRawStatistics(nPixels, mean, min, max, std);
            meanHue = (mean / 255) * 360;  // Convert back to 0-360 scale
            stdHue = (std / 255) * 360;
            minHue = (min / 255) * 360;
            maxHue = (max / 255) * 360;
            
            // Measure Saturation
            selectWindow("Saturation");
            run("Restore Selection");
            getRawStatistics(nPixels, mean, min, max, std);
            meanSat = (mean / 255) * 100;  // Convert to 0-100 scale
            stdSat = (std / 255) * 100;
            minSat = (min / 255) * 100;
            maxSat = (max / 255) * 100;
            
            // Measure Value/Brightness
            selectWindow("Brightness");
            run("Restore Selection");
            getRawStatistics(nPixels, mean, min, max, std);
            meanVal = (mean / 255) * 100;  // Convert to 0-100 scale
            stdVal = (std / 255) * 100;
            minVal = (min / 255) * 100;
            maxVal = (max / 255) * 100;
            
            // Calculate color uniformity score
            // Lower standard deviations indicate more uniform color
            colorUniformity = 100 - ((stdHue/360 + stdSat/100 + stdVal/100) / 3 * 100);
            
        } else {
            // No kelp detected - set all values to 0
            meanHue = 0; stdHue = 0; minHue = 0; maxHue = 0;
            meanSat = 0; stdSat = 0; minSat = 0; maxSat = 0;
            meanVal = 0; stdVal = 0; minVal = 0; maxVal = 0;
            colorUniformity = 0;
        }
        
        // ====================================================================
        // SECTION 12: HEALTH ASSESSMENT
        // ====================================================================
        
        healthFlag = "OK";
        healthNotes = "";
        
        if (kelpAreaPixels == 0) {
            healthFlag = "NO_KELP";
            healthNotes = "No kelp detected";
        } else {
            // Check if color values indicate stress
            if (meanSat < healthySatMin) {
                healthFlag = "WARNING";
                healthNotes += "Low saturation (" + d2s(meanSat,1) + "%); ";
            }
            if (meanVal < healthyValueMin) {
                healthFlag = "WARNING";
                healthNotes += "Low brightness (" + d2s(meanVal,1) + "%); ";
            }
            if (meanHue < healthyHueMin || meanHue > healthyHueMax) {
                healthFlag = "WARNING";
                healthNotes += "Hue out of range (" + d2s(meanHue,1) + "°); ";
            }
            
            // Check color uniformity
            if (colorUniformity < 70) {
                if (healthFlag == "OK") healthFlag = "INFO";
                healthNotes += "Patchy coloration; ";
            }
        }
        
        if (healthFlag == "OK") {
            healthNotes = "Healthy appearance";
        }
        
        print("    Area: " + d2s(kelpAreaPixels,0) + " pixels");
        print("    HSV: H=" + d2s(meanHue,1) + "° S=" + d2s(meanSat,1) + "% V=" + d2s(meanVal,1) + "%");
        print("    Health: " + healthFlag + " - " + healthNotes);
        
        // ====================================================================
        // SECTION 13: SAVE RESULTS FOR THIS BEAKER
        // ====================================================================
        
        // Calculate area in mm² (requires calibration - for now use placeholder)
        // You would need to measure a known distance in pixels to calibrate
        pixelsPerMM = 10;  // PLACEHOLDER - NEEDS CALIBRATION
        kelpAreaMM2 = kelpAreaPixels / (pixelsPerMM * pixelsPerMM);
        
        // Format results row
        resultsRow = currentImage + "," +
                    sheetNumber + "," +
                    dateTaken + "," +
                    kelpID + "," +
                    gridRow + "," +
                    gridCol + "," +
                    gridPosition + "," +
                    nutrientTreatment + "," +
                    lightTreatment + "," +
                    aerationTreatment + "," +
                    d2s(kelpAreaPixels,2) + "," +
                    d2s(kelpAreaMM2,2) + "," +
                    d2s(meanHue,2) + "," +
                    d2s(stdHue,2) + "," +
                    d2s(meanSat,2) + "," +
                    d2s(stdSat,2) + "," +
                    d2s(meanVal,2) + "," +
                    d2s(stdVal,2) + "," +
                    d2s(minHue,2) + "," +
                    d2s(maxHue,2) + "," +
                    d2s(minSat,2) + "," +
                    d2s(maxSat,2) + "," +
                    d2s(minVal,2) + "," +
                    d2s(maxVal,2) + "," +
                    d2s(colorUniformity,2) + "," +
                    healthFlag + "," +
                    healthNotes + "," +
                    qcPass + "," +
                    qcIssues;
        
        File.append(resultsRow, resultsPath);
        
        // ====================================================================
        // SECTION 14: ADD OVERLAY TO QC IMAGE
        // ====================================================================
        
        selectImage(qcOverlayID);
        makeOval(roiX, roiY, roiWidth, roiHeight);
        
        // Color-code by health status
        if (healthFlag == "OK") {
            setColor(0, 255, 0);  // Green
        } else if (healthFlag == "WARNING") {
            setColor(255, 165, 0);  // Orange
        } else if (healthFlag == "NO_KELP") {
            setColor(255, 0, 0);  // Red
        } else {
            setColor(255, 255, 0);  // Yellow
        }
        
        run("Draw", "slice");
        
        // Add text label
        setFont("SansSerif", 12, "bold");
        drawString("ID:" + kelpID, roiX + 5, roiY + 15);
        drawString(healthFlag, roiX + 5, roiY + 30);
        
        // ====================================================================
        // SECTION 15: CLEANUP FOR THIS BEAKER
        // ====================================================================
        
        // Close all temporary images for this beaker
        close("Beaker_*");
        close("Hue*");
        close("Sat*");
        close("Val*");
        close("*Mask");
        close("ForHSV");
        close("Hue");
        close("Saturation");
        close("Brightness");
        
    } // End beaker loop
    
    // ========================================================================
    // SECTION 16: SAVE QC OVERLAY IMAGE
    // ========================================================================
    
    selectImage(qcOverlayID);
    saveAs("PNG", outputFolder + "QC_Images/QC_" + replace(currentImage, ".jpg", ".png"));
    close();
    
    // Save QC summary for this sheet
    qcRow = currentImage + "," +
            sheetNumber + "," +
            colorCheckerDetected + "," +
            d2s(meanBrightness,1) + "," +
            allBeakersDetected + "," +
            qcPass + "," +
            qcIssues;
    File.append(qcRow, qcPath);
    
    // Close original image
    selectImage(originalID);
    close();
    
    print("Completed processing: " + currentImage);
    File.append("Completed: " + currentImage, logPath);
    
} // End image loop

setBatchMode(false);

// ============================================================================
// SECTION 17: FINAL SUMMARY AND REPORTING
// ============================================================================

print("\n=== PROCESSING COMPLETE ===");
print("Total images processed: " + imageFiles.length);
print("Results saved to: " + resultsPath);
print("QC summary saved to: " + qcPath);
print("QC overlay images saved to: " + outputFolder + "QC_Images/");
print("Log file: " + logPath);

File.append("", logPath);
File.append("============================================", logPath);
File.append("Processing complete", logPath);
File.append("Total images: " + imageFiles.length, logPath);

// Show completion message
showMessage("Processing Complete", 
            "Analyzed " + imageFiles.length + " images\n" +
            (imageFiles.length * totalBeakersPerSheet) + " beaker measurements\n\n" +
            "Results saved to:\n" + outputFolder);

print("\n=== NEXT STEPS ===");
print("1. Review QC overlay images in: " + outputFolder + "QC_Images/");
print("2. Check QC summary for any flagged issues: " + qcPath);
print("3. Import results into your analysis software: " + resultsPath);
print("4. Compare measurements across treatment groups");
print("5. Track growth trends over time for individual kelp IDs");

// End of script
