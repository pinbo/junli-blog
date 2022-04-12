let bwa = new Aioli("strobealign/0.7");
let samtools = new Aioli("samtools/latest"); // null before init
// Initialize bwa and output the version
bwa
.init()
.then(() => bwa.exec("-h"))
.then(d => console.log(d.stdout, "ERRRRR", d.stderr));
// init samtools
samtools.init()
.then(() => samtools.exec("--version"))
.then(d => console.log("Samtools: ", d.stdout, "\nSTDERR\n", d.stderr));
// to call variants
let exactSNP = new Aioli("exactSNP/2.0.1");
// exactSNP will be Aioli.workers[1]
exactSNP
.init()
.then(() => exactSNP.exec("-v"))
.then(d => console.log("STDOUT\n", d.stdout, "STDERR\n", d.stderr));
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
            // let bb = samtools.downloadBinary("/data/" + bamfile).then(d => d.arrayBuffer()).then(d => zip.file(bamfile, d));
            // let cc = await Promise.all([aa, bb]);
            promises.push(aa);
        }
    }
    // bams in exactSNP worer
    // files = await exactSNP.ls("/data"); // an array of files
    // for (let i = 0, f; f = files[i]; i++) {
    //     if (f.includes(".bam")) {
    //         console.log("Prepare downloading ", f);
    //         let aa = exactSNP.downloadBinary("/data/" + f).then(d => d.arrayBuffer()).then(d => zip.file("bams/"+f, d));
    //         promises.push(aa);
    //     }
    // }
    let cc = merge_indels().then(d => d.arrayBuffer()).then(d => zip.file("Summary_of_SNPs_and_small_indels.txt", d));
    promises.push(cc);
    await Promise.all(promises);
    console.log("Finished preparing downloanding bams!");
    zip.generateAsync({type:"blob"})
        .then(function(content) {
            saveAs(content, "bams_and_variant_summary_bwamem.zip");
        });
}

// all in one
async function analyzeBam(){
    await makeSam();
    await makeBam();
    // await transferBam();
    await callAll();
    // document.getElementById("download-btn").style.visibility = "visible";
}

// make bams
async function makeBam(){
    await transferSam(); // first transfer all the sam files to samtools worker
    let files = await samtools.ls("/data"); // an array of files
    console.log(files);
    let promises = [];
    for (let i = 0, f; f = files[i]; i++) {
        if (f.includes(".sam")) {
            promises.push(samtoBam("/data/" + f));
        }
    }
    const d = await Promise.all(promises);
    console.log("Finished make bams!")
    document.getElementById("sort").innerHTML = "BAM files have been created. Variants will be called";
}

// sam to bam for single files
// samtools.ls("/data").then(d => console.log(d))
// samtoBam("/data/out_1.sam")
// samtools.exec("index /data/out_1.bam")
async function samtoBam(samfile){ // samfile full path
    console.log("sam file is: ", samfile);
    let bamfile = samfile.replace(".sam", ".bam");
    console.log("bam file is: ", bamfile);
    document.getElementById("sort").innerHTML = "... Making bam file " + bamfile;
    // await samtools.ls("/data").then(console.log)
    let cmd = ["sort", samfile, bamfile].join(' '); // I edited the samtools c file to make the 2nd argument bam output
    console.log(cmd);
    let std = await samtools.exec(cmd);//.then(() => samtools.exec("index " + bamfile));
    console.log("STDERR\n", std.stderr);
    samtools.rm(samfile); // remove sam files after getting the bam files to save space
    // let indexfile = bamfile + ".bai";
    let cmd2 = ["index", bamfile].join(' ');
    console.log(cmd2);
    let std2 = await samtools.exec("index " + bamfile);
    console.log("STDERR\n", std2.stderr);
    return bamfile;
}

// delay before
const delay = ms => new Promise(res => setTimeout(res, ms));

// transfer sam files to samtools /data
async function transferSam(){
    // Initialize samtools and output the version
    // await samtools.init()
    // .then(() => samtools.exec("--version"))
    // .then(d => console.log("Samtools: ", d.stdout, "\nSTDERR\n", d.stderr));
    let files = await bwa.ls("/data"); // an array of files
    for (var i = 0, f; f = files[i]; i++) {
        if (f.includes(".sam")) {
            Aioli.transfer("/data/" + f, "/data/" + f, bwa, samtools);
            Aioli.transfer("/data/" + f, "/data/" + f, bwa, exactSNP);
        }
    }
    let ref = document.getElementById("demoRef").innerHTML;
    Aioli.transfer("/data/" + ref, "/data/" + ref, bwa, exactSNP);
    await delay(1000);
    console.log("Finished transfering files!");
    return 0;
}

document.getElementById("reference").addEventListener("change", loadRef, false);
document.getElementById("fastq").addEventListener("change", loadFq, false);

function loadFq(event)
{
    document.getElementById("demoFq").innerHTML = "";
    // var filenames = [];
    var files = event.target.files;
    for (var i = 0, f; f = files[i]; i++) {
        Aioli.mount(f, null, null, bwa);
        document.getElementById("demoFq").innerHTML += f.name + "\t";
    }
}

// load fasta file
async function loadRef(event)
{
    let files = event.target.files;
    let f = files[0];
    document.getElementById("demoRef").innerHTML = f.name;
    await Aioli.mount(f, null, null, bwa) // return new file path
    .then(() => bwa.ls("/data"))
    .then(d => console.log(d));
    // index
    // await delay(1000); // mount did not really await
    // await bwa.exec("index /data/" + f.name)
    // .then(d => console.log(d.stdout, "End of stdout\n", d.stderr, "End of stderr"));
    // let ff = await bwa.ls("/data");
    // if (ff.length < 8) {
    //     document.getElementById("indexErr").innerHTML = "Indexing the reference FAILED. Please refresh the page to retry!";
    // }
}

// run bwa mem on all fastq files
// loop through all fastq.gz files
async function makeSam(){
    // let reference = document.getElementById("demoRef").innerHTML.split("\t")[0];
    let filenames = document.getElementById("demoFq").innerHTML.split("\t");
    let reference = document.getElementById("demoRef").innerHTML;
    let suffix1 =  document.getElementById("suffix1").value; // R1 suffix
    console.log("FASTQ files\n", filenames);
    console.log(reference);
    let promises = [];
    for (i = 0; i < filenames.length; i++) {
        let ff = filenames[i];
        if (ff.includes(suffix1)) {
            console.log("Processing: ", ff);
            let prefix = ff.replace(suffix1, "");
            promises.push(bwamem(prefix, reference));
        }
    }
    await Promise.all(promises);
    document.getElementById("bwa").innerHTML = "Finished mapping reads.";
}

// make sam file with bwa mem
// bwamem("2", "references.fa").then(d => console.log(d));
async function bwamem (prefix, reference) {
    let suffix1 =  document.getElementById("suffix1").value; // R1 suffix
    let suffix2 = document.getElementById("suffix2").value;
    bwa.setwd("/data/"); // set working directory
    let R1 = prefix + suffix1;
    let R2 = prefix + suffix2;
    let out = "out_" + prefix + ".sam";
    // let rg = "-R \@RG\\tID:" + prefix + "\\tSM:" + prefix; // read group tag
    reference = reference;
    // bwa mem
    // let cmd = ["mem", reference, R1, R2, out].join(' '); // I modifed fastmap.c to use the 4th arguments as output
    let cmd = ["-A 4 -o", out, reference, R1, R2].join(' '); // I modifed fastmap.c to use the 4th arguments as output
    console.log(cmd);
    let std = await bwa.exec(cmd);
    console.log(std.stderr);
    console.log("Finished writing ", out);
    document.getElementById("bwa").innerHTML = "... Mapping " + prefix;
    // return out; 
}

// transfer sam files to samtools /data
async function transferBam(){//transfer refernces and bam files to exactSNP worker
    let ref = document.getElementById("demoRef").innerHTML;
    Aioli.transfer("/data/" + ref, "/data/" + ref, bwa, exactSNP);
    let files = await samtools.ls("/data"); // an array of files
    for (var i = 0, f; f = files[i]; i++) {
        if (f.endsWith(".bam")) {
            Aioli.transfer("/data/" + f, "/data/" + f, samtools, exactSNP);
        }
    }
    await delay(1000);
    console.log("Finished transfering bam files!");
    return 0;
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
    document.getElementById("sort").innerHTML = "All the files have been processed!";
}

// call single bam
async function callVar (bam) {
    let ref = document.getElementById("demoRef").innerHTML;
    let wd = "/data/";
    exactSNP.setwd(wd);
    let out = bam + ".vcf";
    let cmd = ["-i", bam, "-g", ref, "-o", out].join(' ');
    console.log(cmd);
    let std = await exactSNP.exec(cmd);
    document.getElementById("sort").innerHTML = "Finished calling SNPs for " + bam;
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
    let filename = f.replace(".bam.vcf", "");
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