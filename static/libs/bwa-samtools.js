let bwa = new Aioli("bwa2/latest");
let samtools = new Aioli("samtools/latest"); // null before init
// Initialize bwa and output the version
bwa
.init()
.then(() => bwa.exec("index"))
.then(d => console.log(d.stdout, "ERRRRR", d.stderr));

// download all the bams as a zip file
// samtools.downloadBinary("/samtools/examples/out2.bam.bai").then(d => saveAs(d, "download.bam"));
async function downloadBam(){
    document.getElementById("download").innerHTML = "... Preparing files for downloading. It may take some time ...";
    let files = await samtools.ls("/data"); // an array of files
    let zip = new JSZip();
    let promises = [];
    for (let i = 0, f; f = files[i]; i++) {
        if (f.includes(".bam.bai")) {
            // let indexfile = "/data/" + f;
            let bamfile = f.replace(".bai", "");
            console.log("Prepare downloading ", bamfile, " and ", f);
            // let indexblob = samtools.downloadBinary("/data/" + f);
            // let bamblob = samtools.downloadBinary("/data/" + bamfile);
            // zip.file(f, indexblob);
            // zip.file(bamfile, bamblob);
            let aa = samtools.downloadBinary("/data/" + f).then(d => d.arrayBuffer()).then(d => zip.file(f, d));
            let bb = samtools.downloadBinary("/data/" + bamfile).then(d => d.arrayBuffer()).then(d => zip.file(bamfile, d));
            let cc = await Promise.all([aa, bb]);
            promises.push(cc);
        }
    }
    const d = await Promise.all(promises);
    console.log("Finished preparing downloanding bams!");
    zip.generateAsync({type:"blob"})
        .then(function(content) {
            saveAs(content, "indexed_bams.zip");
        });
}

// all in one
async function analyzeBam(){
    await makeSam();
    await makeBam();
    document.getElementById("download-btn").style.visibility = "visible";
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
    document.getElementById("sort").innerHTML = "BAM files have been created. You can click the DOWNLOAD button below.";
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
    let indexfile = bamfile + ".bai";
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
    await samtools.init()
    .then(() => samtools.exec("--version"))
    .then(d => console.log("Samtools: ", d.stdout, "\nSTDERR\n", d.stderr));
    let files = await bwa.ls("/data"); // an array of files
    for (var i = 0, f; f = files[i]; i++) {
        if (f.includes(".sam")) {
            Aioli.transfer("/data/" + f, "/data/" + f, Aioli.workers[0], Aioli.workers[1]);
        }
    }
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
        loadSingleFile(f);
        document.getElementById("demoFq").innerHTML += f.name + "\t";
    }
}
// load fasta file
async function loadRef(event)
{
    let files = event.target.files;
    let f = files[0];
    document.getElementById("demoRef").innerHTML = f.name;
    await Aioli.mount(f) // return new file path
    .then(() => bwa.ls("/data"))
    .then(d => console.log(d));
    // index
    await bwa.exec("index /data/" + f.name)
    .then(d => console.log(d.stdout, "End of stdout\n", d.stderr, "End of stderr"));
    let files = await bwa.ls("/data");
    if (files.length < 8) {
        document.getElementById("indexErr").innerHTML = "Indexing the reference FAILED. Please refresh the page to retry!";
    }
}

// run bwa mem on all fastq files
// loop through all fastq.gz files
async function makeSam(){
    // let reference = document.getElementById("demoRef").innerHTML.split("\t")[0];
    let filenames = document.getElementById("demoFq").innerHTML.split("\t");
    let reference = document.getElementById("demoRef").innerHTML;
    let suffix =  document.getElementById("suffix").value; // R1 suffix
    // bwa index reference
    // let wd = "/bwa2/examples/";
    let wd = "/data/";
    console.log("FASTQ files\n", filenames);
    console.log(reference);
    let promises = [];
    for (i = 0; i < filenames.length; i++) {
        let ff = filenames[i];
        if (ff.includes(suffix)) {
            console.log("Processing: ", ff);
            let prefix = ff.replace(suffix, "");
            promises.push(bwamem(prefix, reference));
        }
    }
    await Promise.all(promises);
    document.getElementById("bwa").innerHTML = "Finished mapping reads.";
}

// make sam file with bwa mem
// bwamem("2", "references.fa").then(d => console.log(d));
async function bwamem (prefix, reference) {
    // let wd = "/bwa2/examples/";
    let suffix =  document.getElementById("suffix").value; // R1 suffix
    let R2suffix = suffix.replace("R1", "R2");
    let wd = "/data/";
    let R1 = wd + prefix + suffix;
    let R2 = wd + prefix + R2suffix;
    let out = wd + "out_" + prefix + ".sam";
    reference = wd + reference;
    // bwa mem
    let cmd = ["mem", reference, R1, R2, out].join(' '); // I modifed fastmap.c to use the 4th arguments as output
    console.log(cmd);
    let std = await bwa.exec(cmd);
    console.log(std.stderr);
    console.log("Finished writing ", out);
    document.getElementById("bwa").innerHTML = "... Mapping " + prefix;
    // return out; 
}
