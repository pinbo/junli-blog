// v2 uses https://cdn.biowasm.com/v2/aioli/latest/aioli.js
// aioli is quite different from v1, so several changes are made in this script

// let CLI = await new Aioli("fastp/0.20.1");

// let CLI = await new Aioli("fastp/0.20.1", {
//     urlAioli: "/tools/aioli/v2.4.0/aioli.worker.js",  // Optional: custom path to aioli.js and aioli.worker.js; for local Aioli development (default=biowasm CDN)
//     printInterleaved: true,                 // Optional: whether `exec()` returns interleaved stdout/stderr; if false, returns object with stdout/stderr keys (default=true)
//     debug: false,                           // Optional: set to true to see console log messages for debugging (default=false)
// });
// let temp = await CLI.exec("fastp -i /fastp/testdata/R1.fq -o filtered-R1.fq.gz");
// console.log(temp);

let CLI = await new Aioli({
    tool: "fastp",
    version: "0.20.1-v2",
    // urlPrefix: "/tools/fastp/0.20.1-v2",
    // urlPrefix: `${window.location.origin}/tools/fastp/0.20.1-v2`,
    // reinit: true, // even true did not release the worker memory
},{
    urlAioli: "/tools/aioli/v2.4.0/aioli.worker.js",  // Optional: custom path to aioli.js and aioli.worker.js; for local Aioli 
});



// load fastq files
document.getElementById("fastq").addEventListener("change", loadFq, false);
async function loadFq(event)
{
    document.getElementById("demoFq").innerHTML = "";
    // First, mount the file(s) to a virtual file system
    const files = event.target.files;
    // The function `.mount()` returns the absolute paths of each file mounted
    const paths = await CLI.mount(files);
    console.log(paths);
    // List files in the current folder
    console.log("ls:", await CLI.ls("."));
    for (var i = 0, f; f = files[i]; i++) {
        document.getElementById("demoFq").innerHTML += f.name + "\t";
    }
}

// print help
function printHelp(){
    CLI.exec("fastp --help")
    .then(d => document.getElementById("help").innerHTML = d.stderr.replace(/(?:\r\n|\r|\n)/g, "<br>"));
}

// clear help content
function clearHelp(){
    document.getElementById("help").innerHTML = "";
}

// delay before
const delay = ms => new Promise(res => setTimeout(res, ms));

// make all bams
async function makeAll(){
    document.getElementById("stderr").value = "";
    let suffix1 =  document.getElementById("suffix1").value; // R1 suffix
    let suffix2 =  document.getElementById("suffix2").value; // R2 suffix
    let filenames = document.getElementById("demoFq").innerHTML.split("\t");
    // document.getElementById("download-btn").style.visibility = "visible"; // in case promise.all did not finish
    document.getElementById("download-btn").style.display = "block";
    // fastp.setwd("/data"); // set working directory
    let promises = [];
    for (let i = 0; i < filenames.length; i++) {
        let ff = filenames[i];
        if (ff){// if not empty string
            if (document.getElementById("interleaved").checked){ // in case blank filenames
                promises.push(filter(ff, "")); // only R1
            } else {
                if (ff.includes(suffix1)) {
                    let R2 = ff.replace(suffix1, suffix2);
                    if (!(filenames.includes(R2))) R2="";
                    promises.push(filter(ff, R2));
                }
            }
        }
    }
    await Promise.all(promises);
    document.getElementById("stdout").innerHTML = "All the files have been processed!";
}
document.getElementById("start").addEventListener("click", makeAll, false);
// run fastp on one pair of fastq(.gz) files
async function filter(read1, read2=""){
    let input = "";
    let output = "";
    let reportfile = "-j " + "fastp_" + read1 + ".json " + " -h " + "fastp_" + read1 + ".html";
    if (read2) {// paired end
        input = "-i " + read1 + " -I " + read2;
        output = "-o filtered_" + read1 + " -O filtered_" + read2;
    } else {
        input = "-i " + read1;
        output = "-o filtered_" + read1;
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
    let cmd = ["fastp", input, interleaved, baseQuality, addopt, output, reportfile].join(" ").replace(/  +/g, ' ');

    // fastp.setwd("/data") // set working directory
    console.log("input is", cmd);
    document.getElementById("stdout").innerHTML = "Filtering " + read1;
    let dd = await CLI.exec(cmd);
    // await delay(100); // mount did not really await
    // document.getElementById("stdout").innerHTML += dd.stderr;
    document.getElementById("stderr").value += "================= Filtering " + read1 + " =================\n" + dd + "\n\n";
    // document.getElementById("download-btn").style.display = "block";
    // return 0;
}

// download all the bams as a zip file
// fastp.downloadBinary("/fastp/examples/out2.bam.bai").then(d => saveAs(d, "download.bam"));
async function download(){
    document.getElementById("stdout").innerHTML = "Preparing downloading file ... It might take a while"
    let files = await CLI.ls("."); // an array of files
    let zip = new JSZip();
    // download running log
    let runningSummary = document.getElementById("stderr").value;
    let blob = new Blob([runningSummary], { type: "text/plain;charset=utf-8" });
    zip.file("fastp-running-log.txt", blob);
    // let promises = [];
    for (let i = 0, f; f = files[i]; i++) {
        if (f.includes("filtered_") || f.includes("fastp_")) {
            console.log("Prepare downloading ", f);
            // let blob = new Blob([ await CLI.cat(f) ]);
            // let blob = await CLI.downloadBinary(f);
            let blob = await CLI.download(f).then(d => d.arrayBuffer());
            zip.file(f, blob);
            // promises.push(aa);
        }
    }
    // const d = await Promise.all(promises);
    console.log("Finished preparing downloanding!");
    zip.generateAsync({type:"blob"})
        .then(function(content) {
            saveAs(content, "filtered_fastqs.zip");
        });
}
document.getElementById("download").addEventListener("click", download, false);