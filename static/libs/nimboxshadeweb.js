let boxshade = new Aioli("nimBoxshadeMini/1.2");
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
    let dna = "";
    let fraction = "-t=" + document.getElementById("fraction").value;
    let outlen = "-w=" + document.getElementById("outlen").value;
    let ifgc = "--ifg=" + hexToRgb(document.getElementById("ifgc").value);
    let ibgc = "--ibg=" + hexToRgb(document.getElementById("ibgc").value);
    let sfgc = "--sfg=" + hexToRgb(document.getElementById("sfgc").value);
    let sbgc = "--sbg=" + hexToRgb(document.getElementById("sbgc").value);
    if (document.getElementById("ruler").checked) ruler = "-r";
    if (document.getElementById("consensus").checked) consensus = "-c";
    if (document.getElementById("seqnum").checked) seqnum = "-n";
    if (document.getElementById("dna").checked) dna = "-d";
    // prepare input from text box
    let fasta = document.getElementById("paste").value.trim();
    if (fasta.startsWith('>')){
        let newfile = {};
        newfile.name = "/data/" + infile;
        newfile.content = fasta;
        boxshade.write(newfile); // cannot use await for write function. never finishs
    }
    let cmd = ["/data/" + infile, ruler, consensus, fraction, outlen, seqnum, dna, "-o=" + outfile, ifgc, ibgc, sfgc, sbgc].join(" ").replace(/  +/g, ' ');

    // fastp.setwd("/data") // set working directory
    console.log("input is", cmd);
    document.getElementById("stdout").innerHTML = "Filtering ...";
    let dd = await boxshade.exec(cmd);
    console.log("stdout :", dd.stdout)
    document.getElementById("stdout").innerHTML = dd.stderr;
    document.getElementById("download-btn").style.display = "block";
}

// clear input sequences
function clearseq() {
    document.getElementById("paste").value = "";
};

// input an example
function putExample(){
    document.getElementById("paste").value = ">seq1\nAYTATLVTPT\n>seq2\nTYKVKFITPE"
}

// download
async function download(){
    // document.getElementById("stdout").innerHTML = "Preparing downloading file ... It might take a while";
    // let result = document.getElementById("stdout").innerHTML;
    // let blob = new Blob([result], { type: "text/plain;charset=utf-8" });
    // saveAs(blob, "KASP-primers.csv");
    boxshade.downloadBinary(outfile).then(d => saveAs(d, outfile.split("/")[2]));
}

// hex color to rgb values
function hexToRgb(hex) {
    var result = /^#?([a-f\d]{2})([a-f\d]{2})([a-f\d]{2})$/i.exec(hex);
//     return result ? {
//       r: parseInt(result[1], 16),
//       g: parseInt(result[2], 16),
//       b: parseInt(result[3], 16),
//     } : null;
    return parseInt(result[1], 16) + "," + parseInt(result[2], 16) + "," + parseInt(result[3], 16);
}
