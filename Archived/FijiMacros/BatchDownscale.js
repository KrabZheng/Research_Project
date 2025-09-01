// Batch BIN 2× (Average), 16-bit in/out, add _DS2 suffix
setBatchMode(true);
inputDir  = getDirectory("Choose input");
out = getDirectory("Choose output");
list = getFileList(inputDir);
for (i=0; i<list.length; i++) {
  if (endsWith(list[i], ".tif") || endsWith(list[i], ".tiff")) {
    open(inputDir+list[i]);
    // ensure 16-bit
    run("16-bit");
    // Bin 2× with averaging
    run("Bin...", "x=2 y=2 bin=Average");
    // (Optional) double the pixel size to keep µm/pixel correct
    unit=getUnit(); pw=getPixelWidth(); ph=getPixelHeight();
    run("Set Scale...", "distance=1 known="+(pw*2)+" pixel=1 unit="+unit+" global");
    // save
    base = File.nameWithoutExtension(list[i]);
    saveAs("Tiff", out+base+"_DS2.tif");
    close();
  }
}
setBatchMode(false);