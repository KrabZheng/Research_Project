// ==== USER SETTINGS ====
inputDir  = "I:/File_Backup/SFA_Exp_18/Separated/H1 H1cse H2 H2cse 4x tcl 4x cse Merged/";   // folder containing series_*.ome.tif
outputDir = "I:/File_Backup/SFA_Exp_18/Processed/H1H2/";          // will be created if missing
rollingBall = 300;                           // ~2–3× radius (62 px) → start at ~200
useCLAHE = false;                            // set true if needed (see block below)
// ========================

File.makeDirectory(outputDir);
list = getFileList(inputDir);
for (i=0; i<list.length; i++) {
    if (endsWith(list[i], ".ome.tif") || endsWith(list[i], ".ome.tiff")) {
        open(inputDir + list[i]);
        // Don’t apply B&C; work in raw intensities
        run("Subtract Background...", "rolling="+rollingBall+" sliding");
        // Optional: gentle bandpass to suppress tiny noise & very large trends
        // run("Bandpass Filter...", "filter_large=250 filter_small=2 suppress=None tolerance=5 autoscale");
        if (useCLAHE) {
            // Block size ~2–3× diameter (124 px) → try 256–384; slope 3.0 is conservative
            run("CLAHE ", "blocksize=320 histogram=256 maximum=3 mask=*None*");
        }
        // Convert to 8-bit to cut size & speed up CP (safe after background correction)
        run("8-bit");
        saveAs("Tiff", outputDir + replace(list[i], ".ome.tif", "_proc.tif"));
        close();
    }
}
