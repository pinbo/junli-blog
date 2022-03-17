let boxshade = new Aioli("boxshade/3.31");
// Initialize boxshade and output the version
boxshade.init()
.then(() => boxshade.exec("-help"))
.then(d => console.log(d.stdout, "ERRRRR", d.stderr));

// function to read a text file async, for loadRef
function readTextFileAsync(file) {
    return new Promise((resolve, reject) => {
        let reader = new FileReader();
        reader.onload = () => {
            resolve(new TextDecoder('utf-8').decode(reader.result, {stream: true}));
        };
        reader.onerror = reject;
        reader.readAsArrayBuffer(file);
    });
}

// class new file
class Newfile {
    constructor(name, content) {
      this.name = name;
      this.content = content;
    }
}


// prepare input file from user's csv file
// let inputFileContent = ""; // global variable for later use
let infile = "";
let outfile = "";

async function prepareInput(evt) {
	var files = evt.target.files; // FileList object
	if(files.length==0) return;
    infile = files[0].name.replace(/ |-/g, "_");
    Aioli.mount(files[0], infile, null, boxshade);// only to worker fastp
    // inputFileContent = await readTextFileAsync(files[0]);
    return 0;
}

// load the SNP file
document.getElementById("snpfile").addEventListener("change", prepareInput, false);

// loop through all fastq.gz files
async function process(){
    outfile = "/data/output";
    // if (document.getElementById("output").value) outfile = "/data/" + document.getElementById("output").value;
    let dna = "";
    let ruler = "";
    let consensus = "";
    let fraction = "-thr=" + document.getElementById("fraction").value;
    if (document.getElementById("dna").checked) dna = "-dna";
    if (document.getElementById("ruler").checked) ruler = "-ruler";
    if (document.getElementById("consensus").checked) consensus = "-cons";
    let sel1 = document.getElementById("box1"); 
    let inputtype = sel1.options[sel1.selectedIndex].value;
    if (inputtype == "guess") inputtype = "";
    let sel2 = document.getElementById("box2"); 
    let outformat = sel2.options[sel2.selectedIndex].value;
    let fileExt = ".rtf";
    if (outformat == "-dev=1") fileExt = ".ps";
    else if (outformat == "-dev=2") fileExt = ".eps";
    outfile += fileExt;

    let cmd = ["-in=/data/" + infile, "-out=" + outfile, dna, ruler, consensus, fraction, inputtype, outformat, "-def"].join(" ").replace(/  +/g, ' ');

    // fastp.setwd("/data") // set working directory
    console.log("input is", cmd);
    document.getElementById("stdout").innerHTML = "Filtering ...";
    let dd = await boxshade.exec(cmd);
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

