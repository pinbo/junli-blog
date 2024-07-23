// IPUAC codes for conversion

const CLI = await new Aioli(["samtools/1.17"]);

// trim white spaces
function mytrim(x) {
    return x.replace(/^\s+|\s+$/g,'');
}

// read fasta lines pasted in textarea with id "input"
function processGeneIDs(){
    let inputString = document.getElementById("input").value.replace(/^\s+|\s+$/g,'');
    let IDstring = ""
	let lines = inputString.split('\n');
	for(let i = 0; i < lines.length; i++){
		let ll = mytrim(lines[i]);
        if (ll){
            if (ll.includes(".")) IDstring += ll + "\n"; // already a transcript ID
            else IDstring += ll + ".1" + "\n"; // if a gene ID, use .1 stranscript ID
        }
	}
	return IDstring;
}

// give examples
function giveExamples(){
    let sel = document.getElementById("box1"); // box1 is selection of database
    let database = sel.options[sel.selectedIndex].value;
    if (database == "CS_cDNA_HC_v1.1") document.getElementById("input").value = "TraesCS5A02G391700.2\nTraesCS5B02G396600";
    else if (database == "CS_cDNA_LC_v1.1") document.getElementById("input").value = "TraesCS7D02G442000LC.1\nTraesCSU02G156900LC";
    else if (database == "Kronos_cDNA_v1.0") document.getElementById("input").value = "TrturKRN6B01G025800.1\nTrturKRN7A01G081900";
}

const testNames = [];
let nprocess = 0; // for region file numbers

async function extractSeq() {
    document.getElementById("output").value = "Fetching... Be patient...";
    nprocess += 1;
    let faURL  = "https://jzseqbucket.s3.us-east-2.amazonaws.com/Kronos.v1.0.all.cds.fa.gz";
    let testdb = "Kronos.v1.0";
    let sel = document.getElementById("box1"); // box1 is selection of database
    let database = sel.options[sel.selectedIndex].value;

    if (database == "CS_cDNA_HC_v1.1"){
        faURL = "https://jzseqbucket.s3.us-east-2.amazonaws.com/IWGSC_v1.1_HC_20170706_cds.fasta.gz";
        testdb = "CS_v1.1_HC";
    } else if (database == "CS_cDNA_LC_v1.1") {
        faURL = "https://jzseqbucket.s3.us-east-2.amazonaws.com/IWGSC_v1.1_LC_20170706_cds.fasta.gz";
        testdb = "CS_v1.1_LC";
    }
    let faiURL = faURL + ".fai";
    let gziURL = faURL + ".gzi";

    // This mounts the URLs lazily on the virtual file system. In other words, no data is downloaded yet.
    if(!(testNames.includes(testdb))){
        testNames.push(testdb);
        await CLI.mount([
            { name: testdb + ".fa.gz", url: faURL },
            { name: testdb + ".fa.gz.fai", url: faiURL },
            { name: testdb + ".fa.gz.gzi", url: gziURL },
        ]);
    }
    // input gene IDs
    // Mount a string to path filename.txt
    let geneIDs = processGeneIDs();
    let regionFile = "regions." + nprocess + ".txt";
    console.log(geneIDs);
    await CLI.mount([{
        name: regionFile,
        data: geneIDs
    }]);
    const cmd = "samtools faidx " + testdb + ".fa.gz -r " + regionFile;
    console.log(cmd);

    if (geneIDs){
        const output = await CLI.exec(cmd);
        document.getElementById("output").value = output;
    } else {
        document.getElementById("output").value = "No gene IDs provided!!!";
    }
}

// clear input box
function clearseq() {
	document.getElementById("input").value = "";
	document.getElementById("output").value = "";
};

document.getElementById("run").addEventListener("click", extractSeq);
document.getElementById("clearseq").addEventListener("click", clearseq);
document.getElementById("example").addEventListener("click", giveExamples);