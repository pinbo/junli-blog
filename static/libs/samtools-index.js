// document.getElementById("reference").addEventListener("change", loadRef, false);
document.getElementById("fastq").addEventListener("change", loadFq, false);

let samtools = new Aioli("samtools/latest");
// Initialize
samtools.init()
.then(() => samtools.exec("index"))
.then(d => console.log("STDOUT", d.stdout, "STDERR", d.stderr));

// make all bams
async function makeAll(){
    let filenames = document.getElementById("demoFq").innerHTML.split("\t");
    let wd = "/data/";
    samtools.setwd(wd);
    let promises = [];
    for (i = 0; i < filenames.length; i++) {
        let ff = filenames[i];
        if (ff.endsWith(".bam")) {
            promises.push(indexBam(ff));
        }
    }
    await Promise.all(promises);
    document.getElementById("bam").innerHTML = "All the files have been processed!";
    document.getElementById("download-btn").style.visibility = "visible";
}

// make single bam
async function indexBam (filename) {
    let cmd = ["index", filename].join(' ');
    console.log(cmd);
    let std = await samtools.exec(cmd);
    document.getElementById("bam").innerHTML = "Finished indexing " + filename;
    document.getElementById("indexErr").innerHTML = std.stderr;
    console.log(std.stderr);
}

function loadFq(event)
{
    document.getElementById("demoFq").innerHTML = "";
    var files = event.target.files;
    for (var i = 0, f; f = files[i]; i++) {
        Aioli.mount(f, null, null, samtools);// only to worker samtools
        document.getElementById("demoFq").innerHTML += f.name + "\t";
    }
}

// delay before
const delay = ms => new Promise(res => setTimeout(res, ms));

// download all the output files as a zip file
// samtools.downloadBinary("/samtools/examples/out2.bam.bai").then(d => saveAs(d, "download.bam"));
async function downloadBam(){
    document.getElementById("download").innerHTML = "... Preparing files for downloading. It may take some time ...";
    let files = await samtools.ls("/data"); // an array of files
    let zip = new JSZip();
    let promises = [];
    for (let i = 0, f; f = files[i]; i++) {
        if (f.includes(".bam.bai")) {
            console.log("Prepare downloading ", f);
            let aa = samtools.downloadBinary("/data/" + f).then(d => d.arrayBuffer()).then(d => zip.file(f, d));
            promises.push(aa);
        }
    }
    const d = await Promise.all(promises);
    console.log("Finished preparing downloanding bam indexes!");
	zip.generateAsync({type:"blob"})
		.then(function(content) {
			saveAs(content, "bams_indexes.zip");
		});
}
