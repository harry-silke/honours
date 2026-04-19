# Kelp Tracking Script - Quick Start Checklist

## Before First Run

### ☐ Software Installation
- [ ] Download and install Fiji (ImageJ) from https://fiji.sc/
- [ ] Test that ImageJ opens successfully

### ☐ Data Preparation
- [ ] Create project folder structure:
  ```
  KelpProject/
  ├── master_sheet.xlsx
  ├── Photos/
  └── (Results/ will be created automatically)
  ```

- [ ] Prepare master sheet Excel file with columns:
  - [ ] Column 1: id (numeric: 1, 2, 3... 124)
  - [ ] Column 2: nutrient treatment (text: "High" or "Low")
  - [ ] Column 3: light treatment (text: "High" or "Low")
  - [ ] Column 4: aeration (text: "Yes" or "No")

- [ ] Export master sheet to CSV:
  - [ ] Open master_sheet.xlsx
  - [ ] File → Save As → CSV format
  - [ ] Save as master_sheet.csv in same folder

- [ ] Organize photo files:
  - [ ] All photos in Photos/ folder
  - [ ] Use consistent naming (recommended: YYYY-MM-DD_SheetXX.jpg)
  - [ ] Verify all photos are clear, well-lit, in focus
  - [ ] Check that ColorChecker Mini is visible in each photo

### ☐ Script Configuration
- [ ] Open Kelp_Growth_Health_Tracker.ijm in text editor
- [ ] Update line ~27: `masterSheetPath = "YOUR_PATH/master_sheet.xlsx";`
- [ ] Update line ~28: `imageFolder = "YOUR_PATH/Photos/";`
- [ ] Update line ~29: `outputFolder = "YOUR_PATH/Results/";`
- [ ] Save the file

### ☐ Initial Calibration (Optional but Recommended)

**Color Threshold Calibration:**
- [ ] Open one representative photo in ImageJ
- [ ] Manually measure HSV values of healthy kelp
- [ ] Update color thresholds in script (lines 40-43)

**Spatial Calibration:**
- [ ] Measure a known distance in pixels (e.g., ruler in photo)
- [ ] Calculate pixels per mm
- [ ] Update `pixelsPerMM` in script (around line 663)

## Running the Script

### ☐ First Test Run
- [ ] Open ImageJ/Fiji
- [ ] Plugins → Macros → Run...
- [ ] Select Kelp_Growth_Health_Tracker.ijm
- [ ] Click Open
- [ ] Watch Log window for progress
- [ ] Wait for completion message

### ☐ Verify Results
- [ ] Check that Results/Data/ folder was created
- [ ] Check that Results/QC_Images/ folder was created
- [ ] Check that Results/Logs/ folder was created
- [ ] Open Kelp_Measurements_[timestamp].csv - should have 12 rows per image
- [ ] Open QC_Summary_[timestamp].csv - should have 1 row per image
- [ ] Look at QC overlay images:
  - [ ] Are circles drawn around all 12 beakers?
  - [ ] Are beaker positions correct?
  - [ ] Are kelp IDs labeled correctly?
  - [ ] Are health flags reasonable?

## Common First-Run Issues

### If NO_KELP appears for all beakers:
- [ ] Color thresholds need adjustment - kelp doesn't match defaults
- [ ] Open one photo manually in ImageJ
- [ ] Convert to HSV (Image → Type → HSB Stack)
- [ ] Measure your actual kelp colors
- [ ] Update thresholds in script

### If beakers aren't detected correctly:
- [ ] Check photo framing - is full grid visible?
- [ ] Adjust margin settings in script (Section 8)
- [ ] Verify beakers aren't cut off at edges

### If wrong IDs assigned to beakers:
- [ ] Verify sheet numbering in filenames
- [ ] Check master CSV has exactly 124 rows
- [ ] Verify ID sequence matches your physical arrangement

### If script gives error about CSV:
- [ ] Make sure you exported Excel to CSV
- [ ] Check CSV is in same folder as .xlsx
- [ ] Verify CSV has comma separators (not semicolons)

## After Successful Run

### ☐ Data Analysis Setup
- [ ] Import results CSV into your analysis software (R, Python, Excel)
- [ ] Create growth curves for individual kelp IDs
- [ ] Group by treatment combinations
- [ ] Calculate statistics

### ☐ Ongoing Monitoring
- [ ] Establish routine photo schedule
- [ ] Use consistent naming convention
- [ ] Review QC images after each batch
- [ ] Flag any anomalies or problems
- [ ] Re-photograph problematic sheets

### ☐ Backup Strategy
- [ ] Regularly backup raw photos
- [ ] Save results CSVs in multiple locations
- [ ] Keep processing logs for troubleshooting

## Calibration Validation

### ☐ Verify Accuracy
- [ ] Manually measure 5-10 beakers in ImageJ
- [ ] Compare to script results
- [ ] Ensure <5% difference in area measurements
- [ ] Check that color values make biological sense

### ☐ Document Your Settings
- [ ] Record final threshold values used
- [ ] Note pixels per mm calibration
- [ ] Save example photos showing what "healthy" looks like
- [ ] Keep notes on any adjustments made

## Support Resources

### Documentation
- [ ] Read full User Guide (Kelp_Tracking_User_Guide.md)
- [ ] Review script comments for detailed explanations
- [ ] Check ImageJ documentation for macro functions

### Troubleshooting
- [ ] Check processing log for error messages
- [ ] Review QC summary for systematic issues
- [ ] Compare results across similar beakers
- [ ] Test on subset of images first

## Notes Section

Use this space to record your specific settings and observations:

**My Color Thresholds:**
- Hue range: _______ to _______
- Saturation minimum: _______
- Value minimum: _______

**My Calibration:**
- Pixels per mm: _______
- Measured from: _______

**My Observations:**
- Typical healthy kelp appearance: _______________________
- Common issues to watch for: _______________________
- Best photography conditions: _______________________

**Custom Modifications:**
- Changes made to script: _______________________
- Reason for changes: _______________________

---

## Ready to Go!

Once all items above are checked, you're ready to process your kelp photos!

**Remember:**
✓ Always review QC images after processing
✓ Validate results with manual checks periodically  
✓ Keep consistent photography conditions
✓ Document any unusual observations
✓ Backup your data regularly

Good luck with your kelp research! 🌿
