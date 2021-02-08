let fastp = new Aioli("fastp/0.20.1");
// Initialize fastp and output the version
fastp
.init();
// .then(() => fastp.exec("--help"))
// .then(d => console.log(d.stdout, "ERRRRR", d.stderr));

// print help
function printHelp(){
    fastp.exec("--help")
    .then(d => document.getElementById("help").innerHTML = d.stderr.replace(/(?:\r\n|\r|\n)/g, "<br>"));
}

// clear help content
function clearHelp(){
    document.getElementById("help").innerHTML = "";
}

// load fastq files
document.getElementById("fastq").addEventListener("change", loadFq, false);

function loadFq(event)
{
    document.getElementById("demoFq").innerHTML = "";
    var files = event.target.files;
    for (var i = 0, f; f = files[i]; i++) {
        loadSingleFile(f);
        document.getElementById("demoFq").innerHTML += f.name + "\t";
    }
    // fastp.ls("/fastp2/examples").then(d => console.log(d));
}

// test wethere 2 init will change things
function loadSingleFile(file)
{
    return Aioli
    .mount(file); // First mount the file
}

// run fastp mem on all fastq files
// loop through all fastq.gz files
async function filter(){
    // let reference = document.getElementById("demoRef").innerHTML.split("\t")[0];
    let filenames = document.getElementById("demoFq").innerHTML.trim().split("\t");
    console.log("Files are: ", filenames);
    let wd = "/data/";
    if (filenames.length < 1 || filenames.length > 2) {
        document.getElementById("error").innerHTML = "ERROR: Incorrect number of fastq files.<br> Please give 1 file for single end (SE) sequencing or 2 files for paired end (PE) sequencing or 1 file for interleaved PE!"
    }
    let input = "";
    let output = "";
    if (filenames.length == 1) {
        input = "-i " + filenames[0];
        output = "-o filtered_" + filenames[0];
    }
    else {
        input = "-i " + filenames[0] + " -I " + filenames[1];
        output = "-o filtered_" + filenames[0] + " -O filtered_" + filenames[1];
    }
    let adapterTim = "-A";
    if (document.getElementById("trimAdapter").checked){adapterTim = ""}
    let interleaved = "";
    if (document.getElementById("interleaved").checked){
        interleaved = "--interleaved_in";
        output = "-o filtered_R1_" + filenames[0] + " -O filtered_R2_" + filenames[1];
    }
    if (document.getElementById("interleaved_out").checked){output = "--stdout"}
    let baseQuality = "-q " + document.getElementById("basequality").value;
    let addopt = document.getElementById("addopt").innerHTML.replace(/(?:\r\n|\r|\n)/g, " ");
    let cmd = [input, interleaved, adapterTim, baseQuality, addopt, output].join(" ");

    fastp.setwd("/data") // set working directory
    console.log("input is", cmd);
    document.getElementById("stdout").innerHTML = "Filtering ...";
    let dd = await fastp.exec(cmd);
    document.getElementById("stdout").innerHTML = dd.stderr;

    if (document.getElementById("interleaved_out").checked){
        let outfilename = "";
        if (filenames[0].includes(".gz")){
            outfilename = "filtered_" + filenames[0];
        } else {
            outfilename = "filtered_" + filenames[0] + ".gz";
        }
        let out2 = new stdInfo();
        out2.name = outfilename;
        out2.content = await pako.gzip(dd.stdout);
        fastp.write(out2);
    }
    document.getElementById("download-btn").style.visibility = "visible";
}

// for writing stdout to file
class stdInfo {
    constructor(name, content) {
        this.name = name;
        this.content = content;
    }
}

// download all the bams as a zip file
// fastp.downloadBinary("/fastp/examples/out2.bam.bai").then(d => saveAs(d, "download.bam"));
async function download(){
    let files = await fastp.ls("/data"); // an array of files
    let zip = new JSZip();
    let promises = [];
    for (let i = 0, f; f = files[i]; i++) {
        if (f.includes("filtered_") || f.includes("fastp.")) {
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