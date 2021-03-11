document.getElementById("reference").addEventListener("change", loadRef, false);
document.getElementById("fastq").addEventListener("change", loadFq, false);

let exactSNP = new Aioli("exactSNP/2.0.1");
// exactSNP will be Aioli.workers[1]
exactSNP
.init()
.then(() => exactSNP.exec("-v"))
.then(d => console.log("STDOUT\n", d.stdout, "STDERR\n", d.stderr));

// make all bams
async function makeAll(){
    let filenames = document.getElementById("demoFq").innerHTML.split("\t");
    document.getElementById("download-btn").style.visibility = "visible"; // in case promise.all did not finish
    let promises = [];
    for (i = 0; i < filenames.length; i++) {
        let ff = filenames[i];
        if (ff.includes(".bam")) {
            console.log("Processing: ", ff);
            promises.push(callVar(ff));
        }
    }
    await Promise.all(promises);
    document.getElementById("bam").innerHTML = "All the files have been processed!";
}

// call single bam
async function callVar (bam) {
    let ref = document.getElementById("demoRef").innerHTML;
    let wd = "/data/";
    exactSNP.setwd(wd);
    let out = bam.replace(".bam", ".variant.vcf");
    let cmd = ["-b -i", bam, "-g", ref, "-o", out].join(' ');
    console.log(cmd);
    let std = await exactSNP.exec(cmd);
    document.getElementById("bam").innerHTML = "Finished calling bam file " + out;
    console.log(std.stderr);
    console.log("Finished writing ", out);
    // document.getElementById("stderr").value += std.stderr + "\n";
}

function loadFq(event) // load bams here
{
    var files = event.target.files;
    for (var i = 0, f; f = files[i]; i++) {
        Aioli.mount(f, null, null, exactSNP);// only to worker exactSNP
        document.getElementById("demoFq").innerHTML += f.name + "\t";
    }
}

async function loadRef(event)
{
    let files = event.target.files;
    let f = files[0];
    document.getElementById("demoRef").innerHTML = f.name;
    await Aioli.mount(f, null, null, exactSNP); // only to worker buildindex
}

// delay before
const delay = ms => new Promise(res => setTimeout(res, ms));
// transfer sam files to samtools /data

// download all the output files as a zip file
// samtools.downloadBinary("/samtools/examples/out2.bam.bai").then(d => saveAs(d, "download.bam"));
async function downloadVar(){
    let files = await exactSNP.ls("/data"); // an array of files
    let zip = new JSZip();
    let promises = [];
    for (let i = 0, f; f = files[i]; i++) {
        if (f.includes(".vcf")) {
            console.log("Prepare downloading ", f);
            let aa = await exactSNP.downloadBinary("/data/" + f).then(d => d.arrayBuffer()).then(d => zip.file("vcf/" + f, d));
            promises.push(aa);
        }
    }
    let cc = merge_indels().then(d => d.arrayBuffer()).then(d => zip.file("Summary_of_SNPs_and_small_indels.txt", d));
    promises.push(cc);
    await Promise.all(promises);
    console.log("Finished preparing downloanding bams!");
	zip.generateAsync({type:"blob"})
		.then(function(content) {
			saveAs(content, "variants_calling_results.zip");
		});
}

// merge all the indel.vcf files
async function merge_indels(){
    let files = await exactSNP.ls("/data"); // an array of files
    let promises = [];
    let indelSummary = "Sample\tGene\tPOS\tREF\tALT\tQUAL\tTotalCoverage\tindelCoverage\tindelPercent\tindelSize\n";
    for (let i = 0, f; f = files[i]; i++) {
        if (f.includes(".vcf")) {
            let aa = await process_indel_vcf(f);
            promises.push(aa);
        }
    }
    await Promise.all(promises);
    indelSummary += promises.join("");
    let blob = new Blob([indelSummary], { type: "text/plain;charset=utf-8" });
    // saveAs(blob, "Summary_of_SNPs_and_indels_less_or_equal_16_bp.csv");
    return blob;
}

// process indel vcf content for only 1 file
async function process_indel_vcf(f){//filename
    let filename = f.replace(".variant.vcf", "");
    let vcf = await exactSNP.cat("/data/" + f);
    let lines = vcf.split(/\r?\n/);
    let summary = "";
    for (let line of lines){
        if (line && !line.includes("#")){
            let ss = line.split(/\t/);
            if (ss[7].includes("MM")){// SNPs
                let ee = ss[7].split(/;/);
                let DP = ee[0].replace("DP=", "");
                let SR = ee[2].replace("MM=", "");
                let pct = String(parseInt(SR) / parseInt(DP) * 100); // percent of indels
                let size = String(ss[4].length - ss[3].length);
                summary += [filename, ss[0], ss[1], ss[3], ss[4], ss[5], DP, SR, pct, size].join('\t') + "\n";
            } else { // indels
                let DP = ss[7].replace("INDEL;DP=", "").split(";SR="); // DP and SR
                let pct = String(parseInt(DP[1]) / parseInt(DP[0]) * 100); // percent of indels
                let size = String(ss[4].length - ss[3].length);
                summary += [filename, ss[0], ss[1], ss[3], ss[4], ss[5], DP[0], DP[1], pct, size].join('\t') + "\n";
            }
        }
    }
    return summary;
}

// download all files at once
function downloadAll(){
    document.getElementById("download").innerHTML = "... Preparing files for downloading. It may take some time ...";
    downloadVar();
    // merge_indels();
    // merge_sv();
    // download_stderr();
}