let velveth = new Aioli("velveth/1.2.10"); // null before init
let velvetg = new Aioli("velvetg/1.2.10"); 
// Initialize 
velveth.init();
velvetg.init();

// print help
function printHelp(){
    velveth.exec("-h") // hashing
    .then(d => document.getElementById("help").innerHTML = d.stdout.replace(/(?:\r\n|\r|\n)/g, "<br>"));
    velvetg.exec("-h") // assembling
    .then(d => document.getElementById("help").innerHTML += "<br><br>" + d.stdout.replace(/(?:\r\n|\r|\n)/g, "<br>"));
}

// clear help content
function clearHelp(){
    document.getElementById("help").innerHTML = "";
}

// load fastq files
function loadFq(event)
{
    document.getElementById("demoFq").innerHTML = "";
    var files = event.target.files;
    for (var i = 0, f; f = files[i]; i++) {
        Aioli.mount(f, null, null, velveth);
        document.getElementById("demoFq").innerHTML += f.name + "\t";
    }
}

document.getElementById("fastq").addEventListener("change", loadFq, false);
// delay before
const delay = ms => new Promise(res => setTimeout(res, ms));
// Function to launch velvetg
async function run() {
	// CLI.mount returns the absolute path of each file mounted
	// const files = document.getElementById("fastq").files;
	// const paths = await CLI.mount(files);
    // let nfiles = paths.length; // number of loaded files
 
    let filenames = document.getElementById("demoFq").innerHTML.trim().split("\t");
    console.log("Files are: ", filenames);
    // await velveth.mkdir("/data/out");
    let nfiles = filenames.length;
    if (nfiles < 1 || nfiles > 2) {
        document.getElementById("error").innerHTML = "ERROR: Incorrect number of fastq files.<br> Please give 1 file for single end (SE) sequencing or 2 files for paired end (PE) sequencing or 1 file for interleaved PE!"
    }
    // parameters for velveth
    let format = "-fastq.gz"; document.getElementById('fileFormat').value;
    console.log("Format is ", format);
    let hashLength = document.getElementById("hashLength").value;
    let readType = "-shortPaired"; // default for paired-end short reads
    let fileLayout = "-separate"; // default for separate read1 and read2 (2 fastq files)
    if (document.getElementById("interleaved").checked)  fileLayout = "-interleaved"; 
    if (nfiles == 1 && fileLayout == "-separate")        readType = "-short";
    let inputFiles = ""; //paths[0] + " " + paths[1]
    if (readType == "-short" || fileLayout == "-interleaved") inputFiles = filenames[0];
    else inputFiles = filenames[0] + " " + filenames[1];

    if ((nfiles == 2 && fileLayout == "-interleaved") )
    {
        alert("Warning: You selected interleaved input fastq but give 2 files, use ONLY the first file");
    }
    let addopt1 = document.getElementById("addopt1").value.replace(/(?:\r\n|\r|\n)/g, " "); // additional options for velveth
    let cmd1 = ["/data", hashLength, format, readType, fileLayout, inputFiles, addopt1].join(" ").replace(/  +/g, ' ').trim();

    // parameters for velvetg
    // velvetg output -cov_cutoff 3 -ins_length 250 -exp_cov 25
    let minCov = document.getElementById("minCov").value; // minimum coverage for assembly output
    let insertSize = document.getElementById("insertSize").value; // average insert size
    let expCov = document.getElementById("expCov").value; // expected coverage
    let addopt2 = document.getElementById("addopt2").value.replace(/(?:\r\n|\r|\n)/g, " "); // additional options for velvetg
    let cmd2 = ["/data -cov_cutoff", minCov, "-ins_length", insertSize, "-exp_cov", expCov, addopt2].join(" ").replace(/  +/g, ' ').trim();

    console.log(cmd1);
    console.log(cmd2);
    document.getElementById("stdout").innerHTML = "Hashing ...";
    // await CLI.ls(".");
    velveth.setwd("/data");
    let r1 = await velveth.exec(cmd1);
    console.log("stdout:\n",r1.stdout, "\n\nstdErr\n",r1.stderr);
    let lsfiles = await velveth.ls("/data");
    console.log(lsfiles);
    // transfer output to velvetg workspace
    // await velvetg.mkdir("/data/out");
    let files = await velveth.ls("/data"); // an array of files
    for (var i = 2, f; f = files[i]; i++) {// 0 and 1 are "." and ".."
        Aioli.transfer("/data/" + f, "/data/" + f, velveth, velvetg);
    }
    await delay(1000);
    console.log("Finished transfering files!");
    // start running velvetg
    velvetg.setwd("/data");
    let r2 = await velvetg.exec(cmd2);
    console.log("stdout:\n",r2.stdout, "\n\nstdErr\n",r2.stderr);
    lsfiles = await velvetg.ls("/data");
    console.log(lsfiles);
    document.getElementById("stdout").innerHTML = r2.stderr;
    document.getElementById("download-btn").style.display = "block";

    // add donwload link
    // let url = await velvetg.download("out/contigs.fa");
    // document.getElementById("downloadLink").href = url;
}

// download the assembly
async function download(){
    await velvetg.downloadBinary("/data/contigs.fa").then(d => saveAs(d, "contigs.fa"));
}

//
document.getElementById("run").addEventListener("click", run);
document.getElementById("printHelp").addEventListener("click", printHelp);
document.getElementById("clearHelp").addEventListener("click", clearHelp);
document.getElementById("download").addEventListener("click", download);
