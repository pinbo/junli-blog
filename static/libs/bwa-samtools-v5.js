// v2: 2022-05-06: add editcall to call potential indels, big deletions and inversions
// v2: 2022-05-30: add editcall c version and deleted exactSNP (due to wrong calls for SNPs)
// v4.1: 2022-12-27: can process single end reads too (like merged reads)
// v4.2: 2023-01-03: add -q option for bwa mem to keep SA mapping quality
// v4.3: 2024-03-03: add -q option for bwa editcall (0.7.17JZv3.1) to filter lower mapping quality reads when calling variants
let bwa = new Aioli("bwa2/0.7.17JZv4.0");
let samtools = new Aioli("samtools/latest"); // null before init
// Initialize bwa and output the version
bwa
.init()
.then(() => bwa.exec("index"))
.then(d => console.log(d.stdout, "ERRRRR", d.stderr));
// init samtools
samtools.init()
.then(() => samtools.exec("--version"))
.then(d => console.log("Samtools: ", d.stdout, "\nSTDERR\n", d.stderr));

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
    let bb = bwa.downloadBinary("/data/formattedReferences.fa").then(d => d.arrayBuffer()).then(d => zip.file("formattedReferences.fa", d));
    promises.push(bb);
    let cc = merge_indels().then(d => d.arrayBuffer()).then(d => zip.file("Summary_of_SNPs_and_small_indels_editcall.txt", d));
    promises.push(cc);
    let dd = bwa.downloadBinary("/data/read_count_table.txt").then(d => d.arrayBuffer()).then(d => zip.file("read_count_table.txt", d));
    promises.push(dd);
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
            // Aioli.transfer("/data/" + f, "/data/" + f, bwa, exactSNP);
            // Aioli.transfer("/data/" + f, "/data/" + f, bwa, editcall);
        }
    }
    // let ref = document.getElementById("demoRef").innerHTML;
    // Aioli.transfer("/data/" + ref, "/data/" + ref, bwa, exactSNP);
    // Aioli.transfer("/data/" + ref, "/data/" + ref, bwa, editcall);
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
    let infile = "formattedReferences.fa";
    document.getElementById("demoRef").innerHTML = infile; //f.name;
    await Aioli.mount(f, "ref.fa", null, bwa) // return new file path
    .then(() => bwa.ls("/data"))
    .then(d => console.log(d));
    // index
    await delay(1000); // mount did not really await
    let refContent = await bwa.cat("/data/ref.fa");
    let newfile = {};
    newfile.name = "/data/formattedReferences.fa";
    newfile.content = formatFasta(refContent);
    bwa.write(newfile);
    await delay(1000);
    await bwa.exec("index /data/" + infile)
    .then(d => console.log(d.stdout, "End of stdout\n", d.stderr, "End of stderr"));
    let ff = await bwa.ls("/data");
    if (ff.length < 8) {
        document.getElementById("indexErr").innerHTML = "Indexing the reference FAILED. Please refresh the page to retry!";
    }
}

// run bwa mem on all fastq files
// loop through all fastq.gz files
async function makeSam(){
    // let reference = document.getElementById("demoRef").innerHTML.split("\t")[0];
    let filenames = document.getElementById("demoFq").innerHTML.split("\t");
    let reference = document.getElementById("demoRef").innerHTML;
    let suffix1 =  document.getElementById("suffix1").value; // R1 suffix
    let suffix2 = document.getElementById("suffix2").value;
    console.log("FASTQ files\n", filenames);
    console.log(reference);
    let promises = [];
    for (i = 0; i < filenames.length; i++) {
        let ff = filenames[i];
        if (ff.includes(suffix1)) {
            console.log("Processing: ", ff);
            let prefix = ff.replace(suffix1, "");
            let R2 = ff.replace(suffix1, suffix2);
            if (!(filenames.includes(R2))) R2="";
            promises.push(bwamem(prefix, ff, R2, reference));
        }
    }
    await Promise.all(promises);
    document.getElementById("bwa").innerHTML = "Finished mapping reads.";
}

// make sam file with bwa mem
// bwamem("2", "references.fa").then(d => console.log(d));
async function bwamem (prefix, R1, R2, reference) {
    // let suffix1 =  document.getElementById("suffix1").value; // R1 suffix
    // let suffix2 = document.getElementById("suffix2").value;
    bwa.setwd("/data/"); // set working directory
    // let R1 = prefix + suffix1;
    // let R2 = prefix + suffix2;
    let out = "out_" + prefix + ".sam";
    // let rg = "-R \@RG\\tID:" + prefix + "\\tSM:" + prefix; // read group tag
    reference = reference;
    // bwa mem
    // let cmd = ["mem", reference, R1, R2, out].join(' '); // I modifed fastmap.c to use the 4th arguments as output
    let cmd = ["mem -q -o", out, reference, R1, R2].join(' ').trim(); // no modify of supp alignments
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
    bwa.setwd("/data/");
    let filenames = await bwa.ls("/data");
    let samFileList = []; // sam files
    document.getElementById("download-btn").style.visibility = "visible"; // in case promise.all did not finish
    let promises = [];
    for (i = 0; i < filenames.length; i++) {
        let ff = filenames[i];
        if (ff.endsWith(".sam")) {
            console.log("Processing: ", ff);
            samFileList.push(ff);
            promises.push(callVar(ff));
        }
    }
    await Promise.all(promises);
    let cmd = "count -l 30 -o read_count_table.txt " + samFileList.join(' ');
    console.log(cmd);
    let std = await bwa.exec(cmd); // get read counts by chromsome for all the sam files
    console.log("Finished writing read_count_table.txt");
    console.log(std.stderr);
    document.getElementById("sort").innerHTML = "All the files have been processed!";
}

// call single sam
async function callVar (sam) {
    let ref = document.getElementById("demoRef").innerHTML;
    let wd = "/data/";
    bwa.setwd(wd);
    let out = sam + ".out.txt";
    let mq = "-q " + document.getElementById("mq").value;
    let cmd = ["editcall", "-b -c 2 -f", ref, mq, "-o", out, sam].join(' ');
    console.log(cmd);
    let std = await bwa.exec(cmd);
    document.getElementById("sort").innerHTML = "Finished calling SNPs for " + sam + "with exactSNP";
    console.log(std.stderr);
    console.log("Finished writing ", out);
    // document.getElementById("stderr").value += std.stderr + "\n";
}

// merge all the indel.vcf files
async function merge_indels(){
    let files = await bwa.ls("/data"); // an array of files
    let promises = [];
    let indelSummary = "Sample\tGene\trefStartPos\trefEndPos\tREF\tALT\tindelSize\ttype\tmutCoverage\tTotalCoverage\tindelPercent\n";
    for (let i = 0, f; f = files[i]; i++) {
        if (f.includes(".out.txt")) {
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
    let filename = f.replace(".sam.out.txt", "");
    let vcf = await bwa.cat("/data/" + f);
    let lines = vcf.trim().split(/\r?\n/);
    let summary = "";
    for (let line of lines){
        summary += [filename, line].join('\t') + "\n";
    }
    return summary;
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
