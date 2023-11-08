//const CLI = await new Aioli(["samtools/1.10"]);
// const baseURL = "http://localhost:4321/tools";
const baseURL = document.location.protocol + '//' + document.location.host + "/" + "tools";
// const CLI = await new Aioli([
// {
//     tool: "velvet",
//     version: "1.2.10",
//     program: "velvetg",
//     urlPrefix: baseURL + "/velvet/1.2.10"
// },
// {
//     tool: "velvet",
//     version: "1.2.10",
//     program: "velveth",
//     urlPrefix: baseURL + "/velvet/1.2.10"
// }
// ], {
//     printInterleaved: false,  // Optional: whether to return interleaved stdout/stderr; if false, returns object with stdout/stderr keys (default: true)
//     debug: true,            // Optional: set to true to see console log messages for debugging (default: false)
// });

// load velveth first, because it uses remove() function which only works in base system.
const CLI = await new Aioli([
    {
        tool: "velveth",
        urlPrefix: baseURL + "/velvet/1.2.10"
    },
    {
        tool: "velvetg",
        urlPrefix: baseURL + "/velvet/1.2.10",
        loading: "lazy"
    },
    ], {
        printInterleaved: false,  // Optional: whether to return interleaved stdout/stderr; if false, returns object with stdout/stderr keys (default: true)
        debug: true,            // Optional: set to true to see console log messages for debugging (default: false)
    });

// print help
function printHelp(){
    CLI.exec("velveth") // hashing
    .then(d => document.getElementById("help").innerHTML = d.stdout.replace(/(?:\r\n|\r|\n)/g, "<br>"));
    CLI.exec("velvetg") // assembling
    .then(d => document.getElementById("help").innerHTML += "<br><br>" + "=".repeat(20) + "<br>" + d.stdout.replace(/(?:\r\n|\r|\n)/g, "<br>"));
}

// clear help content
function clearHelp(){
    document.getElementById("help").innerHTML = "";
}


// Function to launch samtools
async function run() {
    // CLI.mount returns the absolute path of each file mounted
    const files = document.getElementById("fastq").files;
    const paths = await CLI.mount(files);

    let nfiles = paths.length; // number of loaded files
    // ./velveth output 31 -fastq.gz -separate -shortPaired  ../sample_R1_001.fastq.gz  ../sample_R2_001.fastq.gz
    if (nfiles < 1 || nfiles > 2) {
        document.getElementById("error").innerHTML = "ERROR: Incorrect number of fastq files.<br> Please give 1 file for single end (SE) sequencing or 2 files for paired end (PE) sequencing or 1 file for interleaved PE!"
        return 1;
    }
    await CLI.mkdir("outFolder");
    console.log("Current folder is ", await CLI.pwd());
    console.log(await CLI.ls("."));
    let format = "-fastq";
    if (paths[0].endsWith(".gz")) format = "-fastq.gz";
    let hashLength = document.getElementById("hashLength").value;
    let readType = "-shortPaired"; // default for paired-end short reads
    let fileLayout = "-separate"; // default for separate read1 and read2 (2 fastq files)
    if (document.getElementById("interleaved").checked)  fileLayout = "-interleaved"; 
    if (nfiles == 1 && fileLayout == "-separate")        readType = "-short";
    let inputFiles = ""; //paths[0] + " " + paths[1]
    if (readType == "-short" || fileLayout == "-interleaved") inputFiles = paths[0];
    else inputFiles = paths[0] + " " + paths[1];

    if ((nfiles == 2 && fileLayout == "-interleaved") )
    {
        alert("Warning: You selected interleaved input fastq but give 2 files, use ONLY the first file");
    }
    let addopt1 = document.getElementById("addopt1").value.replace(/(?:\r\n|\r|\n)/g, " "); // additional options for velveth
    let cmd1 = ["velveth /shared/data/outFolder", hashLength, format, readType, fileLayout, inputFiles, addopt1].join(" ").replace(/  +/g, ' ').trim();
    // velvetg output -cov_cutoff 3 -ins_length 250 -exp_cov 25
    let minCov = document.getElementById("minCov").value; // minimum coverage for assembly output
    let insertSize = document.getElementById("insertSize").value; // average insert size
    let expCov = document.getElementById("expCov").value; // expected coverage
    let addopt2 = document.getElementById("addopt2").value.replace(/(?:\r\n|\r|\n)/g, " "); // additional options for velvetg
    let cmd2 = ["velvetg /shared/data/outFolder -cov_cutoff", minCov, "-ins_length", insertSize, "-exp_cov", expCov, addopt2].join(" ").replace(/  +/g, ' ').trim();

    console.log(cmd1);
    console.log(cmd2);
    document.getElementById("stdout").innerHTML = "Assembling ...";
    // await CLI.ls(".");
    let r1 = await CLI.exec(cmd1);
    console.log("stdout:\n",r1.stdout, "\n\nstdErr\n",r1.stderr);
    let lsfiles = await CLI.ls("outFolder");
    console.log(lsfiles);
    let r2 = await CLI.exec(cmd2);
    console.log("stdout:\n",r2.stdout, "\n\nstdErr\n",r2.stderr);
    lsfiles = await CLI.ls("outFolder");
    console.log(lsfiles);
    document.getElementById("stdout").innerHTML = r2.stderr;
    document.getElementById("download-btn").style.display = "block";

    // add donwload link
    let url = await CLI.download("outFolder/contigs.fa");
    document.getElementById("downloadLink").href = url;
    let url2 = await CLI.download("outFolder/Log");
    document.getElementById("downloadlog").href = url2;
    let url3 = await CLI.download("outFolder/stats.txt");
    document.getElementById("downloadstats").href = url3;
}

// download the assembly
// async function download(){
//     let url = await CLI.download("outFolder/contigs.fa");
//     document.getElementById("download-btn").href = url;
// }

//
document.getElementById("run").addEventListener("click", run);
// document.getElementById("printHelp").addEventListener("click", printHelp);
document.getElementById("printHelp1").addEventListener("click", function() {
    CLI.exec("velveth") // hashing
    .then(d => document.getElementById("help").innerHTML = d.stdout.replace(/(?:\r\n|\r|\n)/g, "<br>"));
});
document.getElementById("printHelp2").addEventListener("click", function() {
    CLI.exec("velvetg") // hashing
    .then(d => document.getElementById("help").innerHTML = d.stdout.replace(/(?:\r\n|\r|\n)/g, "<br>"));
});
document.getElementById("clearHelp").addEventListener("click", clearHelp);
// export { run, printHelp, clearHelp };