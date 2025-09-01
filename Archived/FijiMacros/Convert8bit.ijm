// Batch 16-bit → 8-bit converter with optional inversion, background subtraction, CLAHE, recursion.
// Saves as OME-TIFF (Bio-Formats Exporter) or plain TIFF.
// Deterministic scaling: 0..65535 → 0..255 (no per-image autoscale).

macro "Batch 16→8-bit (OME-TIFF)" {
    // --- UI ---
    inputDir  = getDirectory("Choose INPUT folder");
    if (inputDir=="") exit("No input folder chosen.");
    outputDir = getDirectory("Choose OUTPUT folder");
    if (outputDir=="") exit("No output folder chosen.");

    Dialog.create("16→8-bit Batch Conversion");
        Dialog.addCheckbox("Recurse into subfolders", false);
        Dialog.addCheckbox("Invert (for dark spheres on light background)", false);
        Dialog.addCheckbox("Apply background subtraction (rolling-ball)", false);
        Dialog.addNumber("Rolling-ball radius (px)", 200);
        Dialog.addCheckbox("Apply CLAHE (contrast limit)", false);
        Dialog.addNumber("CLAHE block size (px)", 320);
        Dialog.addChoice("Save format", newArray("OME-TIFF (Bio-Formats)", "Plain TIFF"), "OME-TIFF (Bio-Formats)");
        Dialog.addChoice("Compression (OME-TIFF only)", newArray("Uncompressed", "LZW"), "Uncompressed");
    Dialog.show();

    doRecurse     = Dialog.getCheckbox();
    doInvert      = Dialog.getCheckbox();
    doBG          = Dialog.getCheckbox();
    rbRadius      = Dialog.getNumber();
    doCLAHE       = Dialog.getCheckbox();
    claheBlock    = Dialog.getNumber();
    saveChoice    = Dialog.getChoice();
    compChoice    = Dialog.getChoice();

    saveAsOME = (saveChoice == "OME-TIFF (Bio-Formats)");
    if (compChoice == "LZW")
        compFlag = "compression=LZW";
    else
        compFlag = "compression=Uncompressed";

    // --- Run ---
    setBatchMode(true);
    nDone = 0; nSkipped = 0; nErrors = 0;

    if (doRecurse)
        processDir(inputDir, true);
    else
        processDir(inputDir, false);

    setBatchMode(false);
    showMessage("Done",
        "Converted: " + nDone + "\n" +
        "Skipped:   " + nSkipped + "\n" +
        "Errors:    " + nErrors + "\n\n" +
        "Output → " + outputDir);
}

// ---- Helpers ----

function processDir(dir, recurse) {
    list = getFileList(dir);
    for (i=0; i<list.length; i++) {
        name = list[i];
        path = dir + name;

        // Skip system/hidden entries
        if (startsWith(name, ".")) continue;

        if (File.isDirectory(path)) {
            if (recurse) processDir(path, true);
            continue;
        }

        // Accept TIFF/OME-TIFF (case-insensitive)
        lname = toLowerCase(name);
        if (endsWith(lname, ".tif") || endsWith(lname, ".tiff")) {
            // Avoid reprocessing outputs (e.g., *_8bit.*)
            if (indexOf(lname, "_8bit.") >= 0) { nSkipped++; continue; }

            showStatus("Processing: " + name);
            runOne(path);
        }
    }
}

function runOne(inPath) {
    // Open with Bio-Formats (safer for OME)
    run("Bio-Formats Importer", "open=[" + inPath + "] autoscale color_mode=Default view=Hyperstack stack_order=Default");

    // Convert to 8-bit deterministically
    depth = bitDepth();
    if (depth==16) {
        setOption("ScaleConversions", true);
        run("8-bit");
    } else if (depth==32) {
        // Map floats into 0..255 via current display range
        setOption("ScaleConversions", true);
        run("8-bit");
    } else if (depth!=8) {
        print("Skipping (unsupported bit depth): " + inPath);
        close();
        nSkipped++;
        return;
    }

    // Optional pre-processing (after 8-bit for speed)
    if (doBG) {
        run("Subtract Background...", "rolling="+rbRadius+" sliding");
    }
    if (doCLAHE) {
        run("CLAHE ", "blocksize="+claheBlock+" histogram=256 maximum=3 mask=*None*");
    }
    if (doInvert) {
        run("Invert");
    }

    // Build output name
    base = File.getName(inPath);
    base = replace(base, ".ome.tiff", "");
    base = replace(base, ".ome.tif", "");
    base = replace(base, ".tiff", "");
    base = replace(base, ".tif", "");
    outName = base + "_8bit";

    if (saveAsOME)
        outPath = outputDir + outName + ".ome.tif";
    else
        outPath = outputDir + outName + ".tif";

    // Save
    if (saveAsOME) {
        run("Bio-Formats Exporter", "save=[" + outPath + "] " + compFlag);
    } else {
        saveAs("Tiff", outPath);
    }

    close();
    nDone++;
}