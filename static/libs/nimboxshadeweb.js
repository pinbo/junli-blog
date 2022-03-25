let boxshade = new Aioli("nimBoxshadeMini/1.0");
// Initialize boxshade and output the version
boxshade.init()
.then(() => boxshade.exec("-h"))
.then(d => console.log(d.stdout, "ERRRRR", d.stderr));


let infile = "input.fa"
let outfile = "/data/output.rtf";

async function prepareInput(evt) {
	var files = evt.target.files; // FileList object
	if(files.length==0) return;
    // infile = files[0].name.replace(/ |-/g, "_");
    Aioli.mount(files[0], infile, null, boxshade);// only to worker boxshade
    // inputFileContent = await readTextFileAsync(files[0]);
    return 0;
}

// load the SNP file
document.getElementById("snpfile").addEventListener("change", prepareInput, false);

// loop through all fastq.gz files
async function process(){
    // if (document.getElementById("output").value) outfile = "/data/" + document.getElementById("output").value;
    let ruler = "";
    let consensus = "";
    let seqnum = "";
    let fraction = "-t=" + document.getElementById("fraction").value;
    let outlen = "-w=" + document.getElementById("outlen").value;
    if (document.getElementById("ruler").checked) ruler = "-r";
    if (document.getElementById("consensus").checked) consensus = "-c";
    if (document.getElementById("seqnum").checked) seqnum = "-n";
    let cmd = ["/data/" + infile, "-o=" + outfile, ruler, consensus, fraction, outlen, seqnum].join(" ").replace(/  +/g, ' ');

    // fastp.setwd("/data") // set working directory
    console.log("input is", cmd);
    document.getElementById("stdout").innerHTML = "Filtering ...";
    let dd = await boxshade.exec(cmd);
    console.log("stdout :", dd.stdout)
    document.getElementById("stdout").innerHTML = dd.stderr;
    document.getElementById("download-btn").style.display = "block";
}

// download
async function download(){
    // document.getElementById("stdout").innerHTML = "Preparing downloading file ... It might take a while";
    // let result = document.getElementById("stdout").innerHTML;
    // let blob = new Blob([result], { type: "text/plain;charset=utf-8" });
    // saveAs(blob, "KASP-primers.csv");
    boxshade.downloadBinary(outfile).then(d => saveAs(d, outfile.split("/")[2]));
}

