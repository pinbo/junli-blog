document.getElementById("reference").addEventListener("change", loadRef, false);
document.getElementById("fastq").addEventListener("change", loadFq, false);

let align = new Aioli("subread-align/2.0.1");
let buildindex = new Aioli("subread-buildindex/2.0.1");
let samtools = new Aioli("samtools/latest"); // for indexing only
// Initialize
buildindex.init();

align.init()
.then(() => align.exec("-v"))
.then(d => console.log("STDOUT", d.stdout, "STDERR", d.stderr));

samtools.init()
.then(() => samtools.exec("index"))
.then(d => console.log("STDOUT", d.stdout, "STDERR", d.stderr));

// make all bams
async function makeAll(){
    let suffix1 =  document.getElementById("suffix1").value; // R1 suffix
    let filenames = document.getElementById("demoFq").innerHTML.split("\t");
    document.getElementById("download-btn").style.visibility = "visible"; // in case promise.all did not finish
    let promises = [];
    for (i = 0; i < filenames.length; i++) {
        let ff = filenames[i];
        if (ff.includes(suffix1)) {
            let prefix = ff.replace(suffix1, "");
            promises.push(makeBam(prefix));
        }
    }
    await Promise.all(promises);
    await transferBam();
    await indexAll(); // wait indexing
    document.getElementById("bam").innerHTML = "All the files have been processed!";
}

// make single bam
async function makeBam (prefix, index="my_index") {
    let wd = "/data/";
    align.setwd(wd);
    let suffix1 =  document.getElementById("suffix1").value; // R1 suffix
    let suffix2 =  document.getElementById("suffix2").value; // R2 suffix
    let R1 = prefix + suffix1;
    let R2 = prefix + suffix2;
    let out = prefix + ".bam";
    let cmd = ["-i", index, "-r", R1, "-R", R2, "-o", out, "-t 1 -I 16 --sv --sortReadsByCoordinates"].join(' ');
    console.log(cmd);
    let std = await align.exec(cmd);
    document.getElementById("bam").innerHTML = "Finished making bam file " + out;
    console.log(std.stderr);
    console.log("Finished writing ", out);
    document.getElementById("stderr").value += std.stderr + "Created " + out + "\n";
    scrollLogToBottom("stderr");
    return 0;
}

function loadFq(event)
{
    document.getElementById("demoFq").innerHTML = "";
    var files = event.target.files;
    for (var i = 0, f; f = files[i]; i++) {
        Aioli.mount(f, null, null, align);// only to worker align
        document.getElementById("demoFq").innerHTML += f.name + "\t";
    }
}

async function loadRef(event)
{
    let files = event.target.files;
    let f = files[0];
    document.getElementById("demoRef").innerHTML = f.name;
    await Aioli.mount(f, null, null, buildindex); // only to worker buildindex
    await delay(500);
    buildindex.ls("/data").then(console.log);
    // index
    buildindex.exec("-M 1000 -o /data/my_index /data/" + f.name)
    .then(std => document.getElementById("stderr").value += std.stderr + "\n");
    let ff = await buildindex.ls("/data");
    if (ff.length < 8) {
        document.getElementById("indexErr").innerHTML = "Indexing the reference FAILED. Please refresh the page to retry!";
    } else {
        scrollLogToBottom("stderr");
        await transferIndex();
        buildindex.worker.terminate();
    }
}

// delay before
const delay = ms => new Promise(res => setTimeout(res, ms));
// transfer sam files to samtools /data
async function transferIndex(){
    let files = await buildindex.ls("/data"); // an array of files
    // let promises = [];
    for (var i = 0, f; f = files[i]; i++) {
        if (f.includes("my_index")) {
            Aioli.transfer("/data/" + f, "/data/" + f, buildindex, align);
        }
    }
    // await Promise.all(promises);
    await delay(1000);
    console.log("Finished transfering index files!");
}

// transfer bam files to samtools worker for making index
// because I found a lot of indexes are broken from subread align
async function transferBam(){
    let files = await align.ls("/data"); // an array of files
    // let promises = [];
    for (var i = 0, f; f = files[i]; i++) {
        if (f.endsWith(".bam")) {
            Aioli.transfer("/data/" + f, "/data/" + f, align, samtools);
        }
    }
    // await Promise.all(promises);
    await delay(1000);
    console.log("Finished transfering bam files!");
}

// index single bams
async function indexBam (filename) {
    let cmd = ["index", filename].join(' ');
    console.log(cmd);
    let std = await samtools.exec(cmd);
    document.getElementById("bam").innerHTML = "Finished indexing " + filename;
    document.getElementById("bamErr").innerHTML = std.stderr;
    console.log(std.stderr);
}

// index all bams
// make all bams
async function indexAll(){
    let wd = "/data/";
    samtools.setwd(wd);
    let files = await samtools.ls("/data"); // an array of files
    let promises = [];
    for (i = 0; i < files.length; i++) {
        let ff = files[i];
        if (ff.endsWith(".bam")) {
            promises.push(indexBam(ff));
        }
    }
    return Promise.all(promises);
    // document.getElementById("bam").innerHTML = "All the files have been processed!";
    // document.getElementById("download-btn").style.visibility = "visible";
}


// download all the output files as a zip file
// samtools.downloadBinary("/samtools/examples/out2.bam.bai").then(d => saveAs(d, "download.bam"));
async function downloadBam(){
    let files = await samtools.ls("/data"); // an array of files
    let fn = files.length; // if no files were transfered promise.all(making bam) was blocked, download files from align /data
    if (fn == 2) files = await align.ls("/data");
    let zip = new JSZip();
    let promises = [];
    for (let i = 0, f; f = files[i]; i++) {
        if (f.includes(".bam")) {
            console.log("Prepare downloading ", f);
            if (fn > 2){
                let aa = samtools.downloadBinary("/data/" + f).then(d => d.arrayBuffer()).then(d => zip.file("bams_by_subread/" + f, d));
                promises.push(aa);
            } else {
                let aa = align.downloadBinary("/data/" + f).then(d => d.arrayBuffer()).then(d => zip.file("bams_by_subread/" + f, d));
                promises.push(aa);
            }
        }
    }
    let bb = download_stderr().arrayBuffer().then(d => zip.file("Subread_running_summary.txt", d));
    promises.push(bb);
    let cc = merge_indels().then(d => d.arrayBuffer()).then(d => zip.file("Summary_of_indels_less_or_equal_16_bp.csv", d));
    promises.push(cc);
    let dd = merge_sv().then(d => d.arrayBuffer()).then(d => zip.file("Summary_of_indels_more_than_16_bp.csv", d));
    promises.push(dd);
    
    await Promise.all(promises);
    console.log("Finished preparing downloanding bams!");
	zip.generateAsync({type:"blob"})
		.then(function(content) {
			saveAs(content, "bams_and_indel_summary_Subread.zip");
		});
}

// download the indexed reference files
async function downloadIndex(){
    let files = await buildindex.ls("/data"); // an array of files
    let zip = new JSZip();
    let promises = [];
    for (let i = 0, f; f = files[i]; i++) {
        if (f.includes("my_index")) {
            console.log("Prepare downloading ", f);
            let aa = await buildindex.downloadBinary("/data/" + f).then(d => d.arrayBuffer()).then(d => zip.file(f, d));
            promises.push(aa);
        }
    }
    const d = await Promise.all(promises);
    console.log("Finished preparing downloanding indexed references!");
	zip.generateAsync({type:"blob"})
		.then(function(content) {
			saveAs(content, "indexed_reference_files.zip");
		});
}

// merge all the indel.vcf files
async function merge_indels(){
    let files = await align.ls("/data"); // an array of files
    let promises = [];
    let indelSummary = "Sample,Gene,POS,REF,ALT,TotalCoverage,indelCoverage,indelPercent,indelSize\n";
    for (let i = 0, f; f = files[i]; i++) {
        if (f.includes(".bam.indel.vcf")) {
            let aa = await process_indel_vcf(f);
            promises.push(aa);
        }
    }
    const d = await Promise.all(promises);
    indelSummary += promises.join("");
    let blob = new Blob([indelSummary], { type: "text/plain;charset=utf-8" });
    // saveAs(blob, "Summary_of_indels_less_or_equal_16_bp.csv");
    return blob;
}

// process indel vcf content for only 1 file
async function process_indel_vcf(f){//filename
    let filename = f.replace(".bam.indel.vcf", "");
    let vcf = await align.cat("/data/" + f);
    let lines = vcf.split(/\r?\n/);
    let summary = "";
    for (let line of lines){
        if (line && !line.includes("#")){
            let ss = line.split(/\t/);
            let DP = ss[7].replace("INDEL;DP=", "").split(";SR="); // DP and SR
            if (DP[1] - DP[0] > 0) DP[0] = DP[1];
            let pct = (parseInt(DP[1]) / parseInt(DP[0]) * 100).toFixed(1); // percent of indels
            let size = String(ss[4].length - ss[3].length);
            summary += [filename, ss[0], ss[1], ss[3], ss[4], DP[0], DP[1], pct, size].join(',') + "\n";
        }
    }
    return summary;
}

// merge all the .bam.breakpoints.vcf files
async function merge_sv(){ // structure variations and big indels > 15bp
    let files = await align.ls("/data"); // an array of files
    let promises = [];
    let indelSummary = "Sample,Gene,Start,End,indelCoverage,indelSize\n";
    for (let i = 0, f; f = files[i]; i++) {
        if (f.includes(".bam.breakpoints.vcf")) {
            let aa = await process_sv_vcf(f);
            promises.push(aa);
        }
    }
    const d = await Promise.all(promises);
    indelSummary += promises.join("");
    let blob = new Blob([indelSummary], { type: "text/plain;charset=utf-8" });
    // saveAs(blob, "Summary_of_indels_more_than_16_bp.csv");
    return blob;
}

// process indel vcf content for only 1 file
async function process_sv_vcf(f){//filename
    let filename = f.replace(".bam.breakpoints.vcf", "");
    let vcf = await align.cat("/data/" + f);
    let lines = vcf.split(/\r?\n/);
    let summary = "";
    let n = 0;
    let start = "";
    for (let line of lines){
        if (line && !line.includes("#")){
            n += 1;
            let ss = line.split(/\t/);
            if(n % 2 == 1){
                summary += [filename, ss[0], ss[1]].join(',') + ",";
                start = ss[1];
            } else {
                let SR = ss[7].split(";SR=")[1]; // SR
                let size = String(parseInt(ss[1]) - parseInt(start) - 1); // indel size
                summary += [ss[1], SR, size].join(',') + "\n";
            }
        }
    }
    return summary;
}

// function to download running summary
function download_stderr(){
    runningSummary = document.getElementById("stderr").value;
    let blob = new Blob([runningSummary], { type: "text/plain;charset=utf-8" });
    // saveAs(blob, "Software_running_summary.txt");
    return blob;
}
// download all files at once
function downloadAll(){
    document.getElementById("download").innerHTML = "... Preparing files for downloading. It may take some time ...";
    downloadBam();
    // merge_indels();
    // merge_sv();
    // download_stderr();
}

// function to show the bottom of the textarea
function scrollLogToBottom(id) {
    logTa = document.getElementById(id)
    logTa.scrollTop = logTa.scrollHeight;
}