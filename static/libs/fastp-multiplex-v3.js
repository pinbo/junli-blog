// fastp worker to load data and store results
let fastp = new Aioli("fastp/0.20.1");
// Initialize fastp and output the version
fastp.init();
// .then(() => fastp.exec("--help"))
// .then(d => console.log(d.stdout, "ERRRRR", d.stderr));

// fastp2 for processing files and terminate after each file to release memory
let fastp2 = new Aioli("fastp/0.20.1");

// print help
function printHelp(){
    fastp.exec("--help")
    .then(d => document.getElementById("help").innerHTML = d.stderr.replace(/(?:\r\n|\r|\n)/g, "<br>"));
}

// clear help content
function clearHelp(){
    document.getElementById("help").innerHTML = "";
}

// delay before
const delay = ms => new Promise(res => setTimeout(res, ms));

// load fastq files
document.getElementById("fastq").addEventListener("change", loadFq, false);

function loadFq(event)
{
    document.getElementById("demoFq").innerHTML = "";
    var files = event.target.files;
    for (var i = 0, f; f = files[i]; i++) {
        Aioli.mount(f, null, null, fastp);// only to worker fastp
        document.getElementById("demoFq").innerHTML += f.name + "\t";
    }
}
// transfer file from
async function transferFile(ff, from, to){//ff is file name
    Aioli.transfer("/data/"+ff, "/data/"+ff, from, to); 
    await delay(5);
}
// transfer filtered files back
async function transferFiltered(){
    let files = await fastp2.ls("/data"); // an array of files
    for (var i = 0, f; f = files[i]; i++) {
        if (f.includes("filtered_")){
            Aioli.transfer("/data/" + f, "/data/" + f, fastp2, fastp);
        }
    }
    await delay(10);
}
// make all bams
async function makeAll(){
    document.getElementById("stderr").value = "";
    let suffix1 =  document.getElementById("suffix1").value; // R1 suffix
    let suffix2 =  document.getElementById("suffix2").value; // R2 suffix
    let filenames = document.getElementById("demoFq").innerHTML.split("\t");
    // document.getElementById("download-btn").style.visibility = "visible"; // in case promise.all did not finish
    document.getElementById("download-btn").style.display = "block";
    // reiniatialize
    // fastp2.init();
    // fastp2.setwd("/data"); // set working directory
    // let promises = [];
    for (i = 0; i < filenames.length; i++) {
        let ff = filenames[i];
        if (ff){// if not empty string
            if (document.getElementById("interleaved").checked){ // in case blank filenames
                await fastp2.init();
                fastp2.setwd("/data");
                await transferFile(ff, fastp, fastp2);
                // promises.push(filter(ff, "")); // only R1
                await filter(ff, "");
                await transferFiltered();
                fastp2.worker.terminate();
            } else {
                if (ff.includes(suffix1)) {
                    await fastp2.init();
                    fastp2.setwd("/data");
                    await transferFile(ff, fastp, fastp2);
                    let R2 = ff.replace(suffix1, suffix2);
                    if (filenames.includes(R2)) {
                        await transferFile(R2, fastp, fastp2);
                    } else R2="";
                    // promises.push(filter(ff, R2));
                    await filter(ff, R2);
                    await transferFiltered();
                    fastp2.worker.terminate();
                }
            }
        }
    }
    // await Promise.all(promises);
    document.getElementById("stdout").innerHTML = "All the files have been processed!";
}

// run fastp on one pair of fastq(.gz) files
async function filter(read1, read2=""){
    let input = "";
    let output = "";
    let reportfile = "-j " + "fastp_" + read1 + ".json " + " -h " + "fastp_" + read1 + ".html";
    if (read2) {// paired end
        input = "-i " + read1 + " -I " + read2;
        output = "-o filtered_" + read1 + " -O filtered_" + read2;
    } else { // single end: there is a bug: cannot trim adapters for .gz file (okay for non-gz file), so disable it
        input = "-i " + read1;
        // a bug: gz file cannot be processed with adapter trimming, so disable it
        if (read1.endsWith(".gz")) output = "-A -o filtered_" + read1; // disable adapter trimming
        else output = "-o filtered_" + read1; // okay for non-gz files
    }
    // let adapterTim = "-A";
    // if (document.getElementById("trimAdapter").checked){adapterTim = ""}
    let interleaved = "";
    if (document.getElementById("interleaved").checked){
        interleaved = "--interleaved_in";
        output = "-o filtered_R1_" + read1 + " -O filtered_R2_" + read1;
    }
    if (document.getElementById("merge").checked){
        if ((read2 === "" && document.getElementById("interleaved").checked) || read2){
            output = output.replace(/filtered/g, "filtered_unmerged") + " -m --merged_out filtered_merged_" + read1;
        } else {
            alert("ERROR: Your input file cannot be merged! Make your input is paired end: read1 and read2 fastq files OR 1 interleaved fastq file!");
            return 1;
        }
    }
    if (document.getElementById("interleaved_out").checked){
        if (document.getElementById("merge").checked){
            document.getElementById("error").innerHTML = "Warning: merged reads cannot be interleaved. Option Ignored.";
            document.getElementById("interleaved_out").checked = false;
        } else if (read2 === "" && !document.getElementById("interleaved").checked) {
            document.getElementById("error").innerHTML = "Warning: Single end reads cannot be interleaved. If your input is an interleaved fastq, please check the 'Interleaved PE'. Option Ignored.";
            document.getElementById("interleaved_out").checked = false;
        } else {
            output = "--interleaved_out -o filtered_interleaved_" + read1;
        }
    }
    let baseQuality = "-q " + document.getElementById("basequality").value;
    let addopt = document.getElementById("addopt").value.replace(/(?:\r\n|\r|\n)/g, " ");
    // let cmd = [input, interleaved, adapterTim, baseQuality, addopt, output].join(" ").replace(/  +/g, ' ');
    let cmd = [input, interleaved, baseQuality, addopt, output, reportfile].join(" ").replace(/  +/g, ' ');

    // fastp.setwd("/data") // set working directory
    console.log("input is", cmd);
    document.getElementById("stdout").innerHTML = "Filtering " + read1;
    let dd = await fastp2.exec(cmd);
    // await delay(100); // mount did not really await
    // document.getElementById("stdout").innerHTML += dd.stderr;
    document.getElementById("stderr").value += "================= Filtering " + read1 + " =================\n" + dd.stderr + "\n\n";
    // document.getElementById("download-btn").style.display = "block";
    // return 0;
}

// download all the bams as a zip file
// fastp.downloadBinary("/fastp/examples/out2.bam.bai").then(d => saveAs(d, "download.bam"));
async function download(){
    document.getElementById("stdout").innerHTML = "Preparing downloading file ... It might take a while"
    let files = await fastp.ls("/data"); // an array of files
    let zip = new JSZip();
    // download running log
    runningSummary = document.getElementById("stderr").value;
    let blob = new Blob([runningSummary], { type: "text/plain;charset=utf-8" });
    zip.file("fastp-running-log.txt", blob)
    let promises = [];
    for (let i = 0, f; f = files[i]; i++) {
        if (f.includes("filtered_")) {
            console.log("Prepare downloading ", f);
            let aa = fastp.downloadBinary("/data/" + f).then(d => d.arrayBuffer()).then(d => zip.file(f, d));
            promises.push(aa);
        }
    }
    const d = await Promise.all(promises);
    console.log("Finished preparing downloanding!");
    zip.generateAsync({type:"blob"})
        .then(function(content) {
            saveAs(content, "filtered_fastqs.zip");
        });
}
