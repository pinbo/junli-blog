document.getElementById("reference").addEventListener("change", loadRef, false);
document.getElementById("fastq").addEventListener("change", loadFq, false);

let align = new Aioli("hisat2-align-s/2.2.1");
let buildindex = new Aioli("hisat2-build-s/2.2.1");
let samtools = new Aioli("samtools/latest"); // for indexing only
let gzip = new Aioli("gzip/1.12");
// Initialize
buildindex.init();

align.init();
// .then(() => align.exec("-v"))
// .then(d => console.log("STDOUT", d.stdout, "STDERR", d.stderr));

samtools.init()
.then(() => samtools.exec("index"))
.then(d => console.log("STDOUT", d.stdout, "STDERR", d.stderr));

gzip.init();
// .then(() => gzip.exec("-h"))
// .then(d => console.log("STDOUT", d.stdout, "STDERR", d.stderr));

// to call variants
let exactSNP = new Aioli("exactSNP/2.0.1");
// exactSNP will be Aioli.workers[1]
exactSNP
.init()
.then(() => exactSNP.exec("-v"))
.then(d => console.log("STDOUT\n", d.stdout, "STDERR\n", d.stderr));

// make all bams
async function makeAll(){
    let suffix1 =  document.getElementById("suffix1").value.replace(".gz",""); // R1 suffix
    let filenames = document.getElementById("demoFq").innerHTML.split("\t");
    document.getElementById("download-btn").style.visibility = "visible"; // in case promise.all did not finish
    let promises = [];
    for (i = 0; i < filenames.length; i++) {
        let ff = filenames[i];
        if (ff.includes(suffix1)) {
            let prefix = ff.replace(suffix1, "");
            promises.push(makeSam(prefix));
        }
    }
    await Promise.all(promises);
    await transferSam();
    await indexAll(); // wait indexing
    await callAll(); // call snps
    document.getElementById("bam").innerHTML = "All the files have been processed!";
}

// make single bam
async function makeSam (prefix, index="my_index") {
    let wd = "/data/";
    align.setwd(wd);
    let suffix1 =  document.getElementById("suffix1").value.replace(".gz",""); // R1 suffix
    let suffix2 =  document.getElementById("suffix2").value.replace(".gz",""); // R2 suffix
    let R1 = prefix + suffix1;
    let R2 = prefix + suffix2;
    let out = prefix + ".sam";
    let addPara = document.getElementById("hisat2").value; // additional parameters
    //-x index/refhisat -1 filtered_R02F01_R1_001.fastq.gz -2 filtered_R02F01_R2_001.fastq.gz -S testhisat.sam
    let cmd = ["-x", index, "-1", R1, "-2", R2, "-S", out, "--pen-noncansplice 0", addPara].join(' ').trim().replace(/  +/g, ' ');
    console.log(cmd);
    let std = await align.exec(cmd);
    document.getElementById("bam").innerHTML = "Finished making sam file " + out;
    console.log(std.stderr);
    console.log("Finished writing ", out);
    document.getElementById("stderr").value += std.stderr + "Created " + out + "\n";
    scrollLogToBottom("stderr");
    return 0;
}

async function loadFq(event)
{
    document.getElementById("demoFq").innerHTML = "";
    var files = event.target.files;
    for (var i = 0, f; f = files[i]; i++) {
        Aioli.mount(f, null, null, gzip);// only to worker align
        // document.getElementById("demoFq").innerHTML += f.name + "\t";
    }
    await delay(1000);
    await unzip();
}

async function unzip(){
    gzip.setwd("/data");
    let cmd = "-d";
    let files = await gzip.ls("/data"); // an array of files
    for (var i = 0, f; f = files[i]; i++) {
        if (f.endsWith("gz"))
          cmd =  cmd + " " + f;
    }
    console.log("gzip", cmd);
    let dd = await gzip.exec(cmd);
    // document.getElementById("bamErr").innerHTML = dd.stderr;
    console.log(dd.stderr);

    files = await gzip.ls("/data"); // an array of files
    // let promises = [];
    for (var i = 2, f; f = files[i]; i++) {
        Aioli.transfer("/data/" + f, "/data/" + f, gzip, align);
        document.getElementById("demoFq").innerHTML += f + "\t";
    }
    // await Promise.all(promises);
    await delay(1000);
    console.log("Finished transfering unzipped fq files!");
    gzip.worker.terminate();
}

async function loadRef(event)
{
    let files = event.target.files;
    let f = files[0];
    document.getElementById("demoRef").innerHTML = "formatted_references.fa";
    await Aioli.mount(f, "ref.fa", null, buildindex); // only to worker buildindex
    await delay(500);
    buildindex.ls("/data").then(console.log);
    //reformat to "formatted_references.fa"
    let refContent = await buildindex.cat("/data/ref.fa");
    let newfile = {};
    newfile.name = "/data/formatted_references.fa";
    newfile.content = formatFasta(refContent);
    buildindex.write(newfile);
    // index
    buildindex.exec("/data/formatted_references.fa /data/my_index")
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
    //transfer fasta file to exactSNP
    Aioli.transfer("/data/formatted_references.fa", "/data/formatted_references.fa", buildindex, exactSNP);
    // await Promise.all(promises);
    await delay(1000);
    console.log("Finished transfering index files!");
}

// transfer bam files to samtools worker for making index
// because I found a lot of indexes are broken from hisat2 align
async function transferSam(){
    let files = await align.ls("/data"); // an array of files
    // let promises = [];
    for (var i = 0, f; f = files[i]; i++) {
        if (f.endsWith(".sam")) {
            Aioli.transfer("/data/" + f, "/data/" + f, align, samtools);
            Aioli.transfer("/data/" + f, "/data/" + f, align, exactSNP);
        }
    }
    // await Promise.all(promises);
    await delay(1000);
    console.log("Finished transfering bam files!");
}

// index single bams
async function indexBam (filename) {
    let bamfile = filename.replace("sam", "bam");
    let cmd = ["sort", filename, "-o", bamfile].join(' ');
    console.log(cmd);
    let std = await samtools.exec(cmd);
    cmd = ["index", bamfile].join(' ');
    console.log(cmd);
    std = await samtools.exec(cmd);
    document.getElementById("bam").innerHTML = "Finished indexing " + bamfile;
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
        if (ff.endsWith(".sam")) {
            promises.push(indexBam(ff));
        }
    }
    return Promise.all(promises);
    // document.getElementById("bam").innerHTML = "All the files have been processed!";
    // document.getElementById("download-btn").style.visibility = "visible";
}


// download all the bams as a zip file
// samtools.downloadBinary("/samtools/examples/out2.bam.bai").then(d => saveAs(d, "download.bam"));
async function downloadBam(){
    document.getElementById("download").innerHTML = "... Preparing files for downloading. It may take some time ...";
    let files = await samtools.ls("/data"); // an array of files
    let zip = new JSZip();
    let promises = [];
    for (let i = 0, f; f = files[i]; i++) {
        if (f.includes(".bam")) {
            console.log("Prepare downloading ", f);
            let aa = samtools.downloadBinary("/data/" + f).then(d => d.arrayBuffer()).then(d => zip.file("bams/"+f, d));
            promises.push(aa);
        }
    }

    let bb = exactSNP.downloadBinary("/data/formatted_references.fa").then(d => d.arrayBuffer()).then(d => zip.file("formatted_references.fa", d));
    promises.push(bb);
    let cc = merge_indels().then(d => d.arrayBuffer()).then(d => zip.file("Summary_of_SNPs_and_indels.txt", d));
    promises.push(cc);
    let dd = download_stderr().arrayBuffer().then(d => zip.file("hisat2_running_summary.txt", d));
    promises.push(dd);
    await Promise.all(promises);
    console.log("Finished preparing downloanding bams!");
    zip.generateAsync({type:"blob"})
        .then(function(content) {
            saveAs(content, "bams_and_variant_summary_hisat2.zip");
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

// call SNPs for all bams
async function callAll(){
    let filenames = await exactSNP.ls("/data");
    document.getElementById("download-btn").style.visibility = "visible"; // in case promise.all did not finish
    let promises = [];
    for (i = 0; i < filenames.length; i++) {
        let ff = filenames[i];
        if (ff.includes(".sam")) {
            console.log("Processing: ", ff);
            promises.push(callVar(ff));
        }
    }
    await Promise.all(promises);
    document.getElementById("bam").innerHTML = "All the files have been processed!";
}

// call single sam
async function callVar (sam) {
    let ref = document.getElementById("demoRef").innerHTML;
    let wd = "/data/";
    exactSNP.setwd(wd);
    let out = sam + ".vcf";
    let addPara = document.getElementById("exactSNP").value; // additional parameters
    let cmd = ["-i", sam, "-g", ref, "-o", out, addPara].join(' ').trim().replace(/  +/g, ' '); // 2022-05-03: I updated aioli.worker, so no need to trim and replace multiple spaces here
    console.log(cmd);
    let std = await exactSNP.exec(cmd);
    document.getElementById("bam").innerHTML = "Finished calling SNPs for " + sam;
    console.log(std.stderr);
    console.log("Finished writing ", out);
    // document.getElementById("stderr").value += std.stderr + "\n";
}
// merge all the indel.vcf files
async function merge_indels(){
    let files = await exactSNP.ls("/data"); // an array of files
    let promises = [];
    let indelSummary = "Sample\tGene\tPOS\tREF\tALT\tTotalCoverage\tmutCoverage\tindelPercent\tindelSize\n";
    for (let i = 0, f; f = files[i]; i++) {
        if (f.includes(".vcf")) {
            let aa = await process_indel_vcf(f);
            promises.push(aa);
        }
    }
    await Promise.all(promises);
    indelSummary += promises.join("");
    let blob = new Blob([indelSummary], { type: "text/plain;charset=utf-8" });
    // saveAs(blob, "Summary_of_SNPs_and_indels_less_or_equal_16_bp.csv");
    return blob;
}

// process indel vcf content for only 1 file
async function process_indel_vcf(f){//filename
    let filename = f.replace(".sam.vcf", "");
    let vcf = await exactSNP.cat("/data/" + f);
    let lines = vcf.split(/\r?\n/);
    let summary = "";
    for (let line of lines){
        if (line && !line.includes("#")){
            let ss = line.split(/\t/);
            if (ss[7].includes("MM")){// SNPs
                let ee = ss[7].split(/;/);
                let DP = ee[0].replace("DP=", ""); // WT counts
                let SR = ee[1].replace("MMsum=", ""); // all mut alleles counts
                let SRsingle = ee[2].replace("MM=", ""); // "3,5" for 2 alt alleles
                let pct = (parseInt(SR) / parseInt(DP) * 100).toFixed(1); // percent of mut
                let size = "0"; // all SNPs
                summary += [filename, ss[0], ss[1], ss[3], ss[4], DP, SRsingle, pct, size].join('\t') + "\n";
            } else { // indels
                let DP = ss[7].replace("INDEL;DP=", "").split(";SR="); // DP and SR
                let pct = (parseInt(DP[1]) / (parseInt(DP[0])) * 100).toFixed(1); // percent of indels
                let size = String(ss[4].length - ss[3].length);
                summary += [filename, ss[0], ss[1], ss[3], ss[4], parseInt(DP[0]), DP[1], pct, size].join('\t') + "\n";
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


// function to show the bottom of the textarea
function scrollLogToBottom(id) {
    logTa = document.getElementById(id)
    logTa.scrollTop = logTa.scrollHeight;
}

// function to format fasta to fixed length
function formatFasta (fileContent, width = 60) {// 60 bp per line
    let lines = fileContent.split(/\r?\n/);
    let dna = "";   // all the bases of 1 seq (removed new lines)
    let title = ""; // title of seq
    let line = "";  // content of each line
    let i = 0;      // line number
    let newContent = "";
    for (;;){
        if(i==lines.length || (line=lines[i].trim())[0]=='>'){// trim to avoid space at the beginning
            if(dna.length != 0){
                newContent += title + "\n";
                while(dna.length != 0){
                    let n = Math.min(dna.length, width);
                    newContent += dna.substring(0,n) + "\n";
                    dna = dna.substring(n);
                }
            }
            if(i===lines.length) break;
            title = line.replace(/> +/, '>');// in case something like "> seq1"
            dna = "";
        } else dna += line.trim().replace(/ +/g, ""); // in case space inside
        ++i;
    }
    return newContent;
}