# Kelp Gametophyte Growth and Health Tracking Script - User Guide

## Overview

This ImageJ macro script automatically analyzes photographs of kelp specimens in a 3x4 beaker grid, tracking growth metrics and health indicators through HSV color analysis. The script is designed for batch processing of multiple photo sheets across time series experiments.

---

## Table of Contents

1. [Prerequisites](#prerequisites)
2. [Initial Setup](#initial-setup)
3. [File Organization](#file-organization)
4. [Running the Script](#running-the-script)
5. [Understanding the Output](#understanding-the-output)
6. [Troubleshooting](#troubleshooting)
7. [Customization Guide](#customization-guide)
8. [Script Components Explained](#script-components-explained)

---

## Prerequisites

### Software Requirements
- **ImageJ** or **Fiji** (recommended - Fiji is ImageJ with useful plugins pre-installed)
- Download from: https://fiji.sc/ or https://imagej.net/

### Data Requirements
- Digital photographs of your 3x4 beaker grids (JPG, PNG, or TIF format)
- Each photo must include a ColorChecker Mini for color calibration
- Master spreadsheet (Excel) containing kelp ID and treatment information
- Consistent photography setup (same camera, lighting, distance, angle)

---

## Initial Setup

### Step 1: Export Your Master Sheet to CSV

The ImageJ macro language cannot directly read Excel files, so you need to export your master spreadsheet:

1. Open your master Excel file
2. Go to **File → Save As**
3. Choose **CSV (Comma delimited) (*.csv)** format
4. Save it with the same name in the same location
   - Example: If your Excel file is `master_sheet.xlsx`, save as `master_sheet.csv`

**Required columns in your CSV:**
- Column 1: `id` (kelp specimen ID, e.g., 1, 2, 3...)
- Column 2: `nutrient treatment` (e.g., "High", "Low")
- Column 3: `light treatment` (e.g., "High", "Low")
- Column 4: `aeration` (e.g., "Yes", "No")

### Step 2: Organize Your Photo Files

Create a folder structure like this:
```
KelpProject/
├── master_sheet.xlsx
├── master_sheet.csv
├── Photos/
│   ├── 2026-02-05_Sheet01.jpg
│   ├── 2026-02-05_Sheet02.jpg
│   ├── 2026-02-12_Sheet01.jpg
│   └── ... (more photos)
└── Results/
    └── (will be created automatically)
```

**Photo naming convention (recommended):**
- Include date in YYYY-MM-DD format: `2026-02-05_Sheet01.jpg`
- Include sheet number: `Sheet01`, `Sheet02`, etc.
- This helps the script automatically extract metadata

### Step 3: Configure the Script

Open the script file `Kelp_Growth_Health_Tracker.ijm` in a text editor and modify the file paths in **SECTION 1: USER CONFIGURATION**:

```javascript
// Update these paths to match your setup:
masterSheetPath = "C:/KelpProject/master_sheet.xlsx";
imageFolder = "C:/KelpProject/Photos/";
outputFolder = "C:/KelpProject/Results/";
```

**Important:** Use forward slashes `/` even on Windows, or double backslashes `\\`

---

## File Organization

### How Kelp IDs Map to Grid Positions

The script assumes your beakers are organized sequentially across sheets:

**Sheet 1:**
- Position 1 (Row 0, Col 0) = Kelp ID 1
- Position 2 (Row 0, Col 1) = Kelp ID 2
- Position 3 (Row 0, Col 2) = Kelp ID 3
- Position 4 (Row 0, Col 3) = Kelp ID 4
- Position 5 (Row 1, Col 0) = Kelp ID 5
- ... continuing through position 12

**Sheet 2:**
- Position 1 = Kelp ID 13
- Position 2 = Kelp ID 14
- ... continuing through position 12 = Kelp ID 24

And so on for all 11 sheets (124 beakers total).

If your organization is different, you'll need to modify the ID mapping logic in Section 9 of the script.

---

## Running the Script

### Method 1: Run from ImageJ Interface

1. Open ImageJ/Fiji
2. Go to **Plugins → Macros → Run...**
3. Navigate to and select `Kelp_Growth_Health_Tracker.ijm`
4. Click **Open**
5. The script will run automatically - watch the Log window for progress

### Method 2: Install as Plugin (for repeated use)

1. Copy `Kelp_Growth_Health_Tracker.ijm` to your ImageJ plugins folder:
   - On Windows: `C:/Users/YourName/ImageJ/plugins/`
   - On Mac: `/Applications/Fiji.app/plugins/`
2. Restart ImageJ
3. The script will appear in **Plugins → Macros → Kelp Growth Health Tracker**

### What Happens During Execution

The script will:
1. Load your master sheet data (CSV file)
2. Find all image files in your Photos folder
3. Process each image:
   - Detect the ColorChecker for calibration
   - Identify the 3x4 beaker grid
   - Analyze each beaker individually
   - Segment kelp from background
   - Measure area and HSV color properties
   - Assess health based on color thresholds
   - Generate QC overlay images
4. Save all results to CSV files
5. Display a summary when complete

**Processing time:** Approximately 30-60 seconds per image, depending on your computer

---

## Understanding the Output

The script creates several output files in your Results folder:

### 1. Main Results File
**Location:** `Results/Data/Kelp_Measurements_YYYY-MM-DD_HHMM.csv`

This is your primary data file with one row per beaker per image. Columns include:

**Image Information:**
- `Image_File`: Name of the photo file
- `Sheet_Number`: Which sheet (extracted from filename)
- `Date_Taken`: Date photo was taken (extracted from filename)

**Beaker Information:**
- `Kelp_ID`: Unique identifier for this kelp specimen
- `Grid_Row`, `Grid_Col`: Position in the 3x4 grid
- `Grid_Position`: 1-12 position number

**Treatment Information:**
- `Nutrient_Treatment`: High or Low
- `Light_Treatment`: High or Low
- `Aeration_Treatment`: Yes or No

**Growth Measurements:**
- `Kelp_Area_pixels`: Total kelp area in pixels
- `Kelp_Area_mm2`: Total kelp area in mm² (requires calibration)

**Color Analysis (HSV values):**
- `Mean_Hue`: Average hue (0-360°) - indicates color (brown vs green vs yellow)
- `StdDev_Hue`: Variation in hue - higher = more color variation
- `Mean_Saturation`: Average saturation (0-100%) - how vivid the color is
- `StdDev_Saturation`: Variation in saturation
- `Mean_Value`: Average brightness (0-100%)
- `StdDev_Value`: Variation in brightness
- `Min_Hue`, `Max_Hue`: Range of hue values
- `Min_Saturation`, `Max_Saturation`: Range of saturation
- `Min_Value`, `Max_Value`: Range of brightness
- `Color_Uniformity_Score`: 0-100 score (higher = more uniform color)

**Health Assessment:**
- `Health_Flag`: OK, WARNING, NO_KELP, or INFO
- `Health_Notes`: Explanation of health status

**Quality Control:**
- `QC_Pass`: TRUE or FALSE
- `QC_Issues`: Description of any problems detected

### 2. Quality Control Summary
**Location:** `Results/Data/QC_Summary_YYYY-MM-DD_HHMM.csv`

One row per image processed, showing:
- Whether ColorChecker was detected
- Mean image brightness
- Whether all beakers were detected
- Overall QC pass/fail status
- List of any issues

### 3. QC Overlay Images
**Location:** `Results/QC_Images/QC_[original_filename].png`

Visual verification images showing:
- Circles drawn around each detected beaker
- Color-coded by health status:
  - **Green** = Healthy (OK)
  - **Orange** = Warning (potential issues)
  - **Red** = No kelp detected
  - **Yellow** = Info (minor notes)
- Kelp ID and health flag labeled on each beaker

**Always review these images** to verify the script correctly identified beaker positions!

### 4. Processing Log
**Location:** `Results/Logs/Processing_Log_YYYY-MM-DD_HHMM.txt`

Detailed log of script execution with timestamps and any errors or warnings.

---

## Troubleshooting

### Common Issues and Solutions

#### Problem: "Master sheet CSV not found"
**Solution:** Make sure you've exported your Excel file to CSV format with the same filename in the same location.

#### Problem: All health flags show "NO_KELP"
**Solutions:**
- Your kelp might not match the default color thresholds
- Adjust `healthyHueMin`, `healthyHueMax`, `healthySatMin`, and `healthyValueMin` in Section 1
- Start by expanding the ranges, then narrow them based on your actual kelp appearance

#### Problem: Beakers detected in wrong positions
**Solutions:**
- Verify your photos show the full 3x4 grid
- Check that beakers aren't cut off at edges
- Adjust `marginX` and `marginY` values in Section 8 if needed
- Ensure consistent photo framing across all images

#### Problem: ColorChecker not detected
**Note:** The current script version has placeholder ColorChecker detection. You may need to:
- Manually define the ColorChecker region in your first image
- Use ImageJ's ROI Manager to save ColorChecker positions
- Or implement more sophisticated ColorChecker detection (this would require additional code)

#### Problem: Incorrect Kelp IDs assigned to beakers
**Solutions:**
- Verify your sheet numbering in photo filenames matches your expectations
- Check that your master CSV has 124 rows (one per beaker)
- Review the ID mapping logic in Section 9 of the script

#### Problem: Script runs very slowly
**Solutions:**
- The script uses batch mode for speed - make sure it's enabled
- Process fewer images at once if memory is limited
- Close other programs to free up RAM
- Reduce image resolution if photos are extremely large (>4000x3000 pixels)

---

## Customization Guide

### Adjusting Color Thresholds for Your Kelp

The default thresholds may not perfectly match your kelp species or photography conditions. To calibrate:

1. **Measure your healthy kelp manually:**
   - Open one of your photos in ImageJ
   - Draw a selection around healthy-looking kelp
   - Go to **Image → Type → HSB Stack**
   - Go to **Image → Stacks → Stack to Images**
   - Measure each channel: **Analyze → Measure**
   - Note the mean values for Hue, Saturation, and Brightness

2. **Convert ImageJ values to script values:**
   - Hue: Multiply ImageJ value by (360/255)
   - Saturation: Multiply ImageJ value by (100/255)
   - Value: Multiply ImageJ value by (100/255)

3. **Update script thresholds in Section 1:**
   ```javascript
   healthyHueMin = [your minimum hue];
   healthyHueMax = [your maximum hue];
   healthySatMin = [your minimum saturation];
   healthyValueMin = [your minimum brightness];
   ```

### Setting Up Spatial Calibration (Pixels to mm)

To convert pixel measurements to real-world dimensions:

1. **Include a scale bar in your photos:**
   - Place a ruler or object of known size in the same plane as the beakers
   - Keep it in the same position for all photos

2. **Measure pixels per mm:**
   - Open a photo in ImageJ
   - Draw a line along your known distance (e.g., 10mm on a ruler)
   - Go to **Analyze → Measure** (or press M)
   - Note the "Length" value in pixels
   - Calculate: `pixelsPerMM = Length_in_pixels / Known_distance_in_mm`

3. **Update script in Section 13:**
   ```javascript
   pixelsPerMM = [your calculated value];  // Replace placeholder
   ```

### Modifying Grid Detection

If your beakers have unusual spacing or you need different margins:

**In Section 8**, adjust these values:
```javascript
marginX = width * 0.05;  // Increase for more left/right margin
marginY = height * 0.05;  // Increase for more top/bottom margin

// Adjust beaker size within cells:
beakerDiameter = minOf(cellWidth, cellHeight) * 0.8;  // Change 0.8 to 0.7 for smaller, 0.9 for larger
```

---

## Script Components Explained

### Section-by-Section Breakdown

**Sections 1-2: Setup**
- User configuration variables
- File system initialization
- Creates output folders

**Section 3: Data Loading**
- Reads your master CSV file
- Parses treatment information
- Creates arrays indexed by kelp ID

**Section 4: Image Discovery**
- Finds all JPG, PNG, TIF files in your Photos folder
- Prepares batch processing

**Section 5: Main Loop**
- Processes each image one at a time
- Extracts metadata from filenames

**Section 6: Quality Checks**
- Verifies image isn't too dark or too bright
- Flags potential photography problems

**Section 7: ColorChecker Detection**
- Currently a placeholder
- In full implementation, would detect and extract reference colors

**Section 8: Grid Detection**
- Calculates positions of all 12 beakers
- Creates circular ROIs for each beaker position

**Section 9: Beaker Loop**
- Processes each of the 12 beakers individually
- Maps grid position to kelp ID using master sheet

**Section 10: Segmentation**
- Converts to HSV color space
- Creates masks based on color thresholds
- Identifies kelp pixels vs background

**Section 11: Color Measurement**
- Measures HSV values only within kelp regions
- Calculates means, standard deviations, min/max
- Computes color uniformity score

**Section 12: Health Assessment**
- Compares color values to healthy ranges
- Sets health flags (OK, WARNING, NO_KELP, INFO)
- Generates health notes

**Section 13: Data Export**
- Formats all measurements for this beaker
- Appends row to results CSV

**Section 14: Visualization**
- Draws colored circles around beakers on QC image
- Adds text labels

**Section 15: Cleanup**
- Closes temporary images
- Frees memory

**Section 16: Image Completion**
- Saves QC overlay image
- Records QC summary

**Section 17: Final Summary**
- Reports total processing statistics
- Shows completion message

---

## Advanced Usage Tips

### Tracking Growth Over Time

To analyze growth trends for individual kelp specimens:

1. Process all your time series photos
2. Import the results CSV into R, Python, or Excel
3. Filter by `Kelp_ID` to get all measurements for one specimen
4. Plot metrics (Area, HSV values) vs `Date_Taken`
5. Calculate growth rates between time points

### Comparing Treatment Effects

To test your experimental hypotheses:

1. Group data by treatment combinations
2. Calculate mean and standard error for each treatment
3. Perform statistical tests (ANOVA, t-tests) on:
   - Growth rates (change in area over time)
   - Final size at end of experiment
   - Color metrics (saturation, hue as health proxies)

### Quality Control Workflow

**Recommended QC process:**
1. Run the script on all images
2. Review QC Summary CSV - identify any failed images
3. Look at QC overlay images for failed/flagged sheets
4. Re-photograph problematic sheets if needed
5. Manually verify a random sample of measurements
6. Remove or flag poor-quality data points in your analysis

---

## Citation and Support

**Script created:** February 2026  
**For:** Kelp gametophyte and sporophyte growth tracking research

### Acknowledgments

This script uses ImageJ/Fiji for image processing. If you publish research using this tool, please cite:

Schindelin, J., Arganda-Carreras, I., Frise, E., et al. (2012). "Fiji: an open-source platform for biological-image analysis." Nature Methods 9(7): 676-682.

### Getting Help

If you encounter issues not covered in this guide:

1. Check the ImageJ/Fiji documentation: https://imagej.net/imaging
2. Review the ImageJ macro language reference: https://imagej.net/ij/developer/macro/functions.html
3. Post questions on the ImageJ forum: https://forum.image.sc/

---

## Appendix: Quick Reference


### Color Thresholds to Adjust
```javascript
healthyHueMin = 20;       // Yellow-brown minimum
healthyHueMax = 60;       // Yellow-brown maximum  
healthySatMin = 30;       // Minimum color vividness
healthyValueMin = 20;     // Minimum brightness
```

### Calibration to Set
```javascript
pixelsPerMM = 10;  // MUST BE CALIBRATED for your setup
```

### Grid Settings to Modify
```javascript
gridRows = 3;
gridCols = 4;
marginX = width * 0.05;
marginY = height * 0.05;
beakerDiameter = minOf(cellWidth, cellHeight) * 0.8;
```

---

## Version History

**Version 1.0** (February 2026)
- Initial release
- HSV color analysis
- Batch processing support
- Health assessment based on color metrics
- QC overlay image generation
- Treatment integration from master sheet

---

*End of User Guide*
