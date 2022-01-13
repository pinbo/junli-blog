let primer3 = new Aioli("primer3_core/2.5.0");
// Initialize primer3 and output the version
primer3
.init();
// .then(() => primer3.exec("-a"))
// .then(d => console.log(d.stdout, "ERRRRR", d.stderr));

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

// class primer
class Primer {
    constructor() {
        this.name = "";
        this.length = "";
        this.start = "";
        this.tm = "";
        this.gc = "";
        this.anyTH = "";
        this.endTH = "";
        this.hairpin = "";
        this.seq = "";
      }
}

class PrimerPair {
    constructor() {
        this.left = new Primer();
        this.right = new Primer();
        this.product_size = "";
        this.snp1 = ""; // the 1st SNP
        this.snp2 = ""; // the 2nd SNP
        this.pairAnyTH = "";
        this.pairEndTH = "";
    }
}

// prepare input file from user's csv file
let inputFileContent = ""; // global variable for later use
async function prepareInput(evt) {
	var files = evt.target.files; // FileList object
	if(files.length==0) return;
    // Aioli.mount(files[0], null, null, primer3);// only to worker fastp
    inputFileContent = await readTextFileAsync(files[0]);
    return 0;
}




// process the output of primer3 
async function designPrimer() {
    // first, prepaire primer3 input file
    console.log("prepare primer3 input");
    let primer3input = await csv2primer3(inputFileContent);
    const p3input = new Newfile("/data/p3input", primer3input);
    console.log("writing primer3 input");
    primer3.write(p3input);
    console.log(await primer3.ls("/data"));
    // start processing
    document.getElementById("download-btn").style.display = "block";
    let dd = await primer3.exec("/data/p3input");
    if (dd.stdout) { // if successful run
        await parse_primer3output(dd.stdout);
        console.log("success!!!");
        console.log("success!!! but with err: ", dd.stderr);
    } else {
        console.log(dd.stderr);
        document.getElementById("error").innerHTML = "Failed to design primers. Please check your input file format and refresh the window to try again.";
    }
    return 0;
}

// process the output of primer3
async function parse_primer3output(primer3output) {
    let lines = primer3output.trim().split(/\n/);
    // let primerpairList = [];
    let primerpairs = {};
    let seqID = "";
    let direction = "";
    let snp2 = "";
    for (let i = 0; i < lines.length; i++) {
        let line = lines[i];
        if (line.startsWith("SEQUENCE_ID")) {
            console.log("sequence_ID is ", line);
            [seqID, direction, snp1, snp2] = line.split("=")[1].split("__");
        } else if (line.match("^PRIMER_.*_SEQUENCE")) {
            let fields = line.split("_");
            let ppn = fields[2]; // primer pair number
            let LR = fields[1]; // left or right
            let primerpairID = seqID + "-" + direction + "-" + ppn;
            console.log("init primerpairID is", primerpairID);
            if (!(primerpairID in primerpairs)) primerpairs[primerpairID] = new PrimerPair();
            primerpairs[primerpairID].snp1 = snp1;
            primerpairs[primerpairID].snp2 = snp2;
            if (LR == "LEFT") primerpairs[primerpairID].left.seq = line.split("=")[1];
            else primerpairs[primerpairID].right.seq = line.split("=")[1];
        } else if (line.match("GC_PERCENT")) {
            let fields = line.split("_");
            let ppn = fields[2]; // primer pair number
            let LR = fields[1]; // left or right
            let primerpairID = seqID + "-" + direction + "-" + ppn;
            if (LR == "LEFT") primerpairs[primerpairID].left.gc = line.split("=")[1];
            else primerpairs[primerpairID].right.gc = line.split("=")[1];
        } else if (line.match("[0-9]+_TM")) {
            let fields = line.split("_");
            let ppn = fields[2]; // primer pair number
            let LR = fields[1]; // left or right
            let primerpairID = seqID + "-" + direction + "-" + ppn;
            if (LR == "LEFT") primerpairs[primerpairID].left.tm = line.split("=")[1];
            else primerpairs[primerpairID].right.tm = line.split("=")[1];
        } else if (line.match("[0-9]+_HAIRPIN")) { // PRIMER_LEFT_0_HAIRPIN_TH=0.00
            let fields = line.split("_");
            let ppn = fields[2]; // primer pair number
            let LR = fields[1]; // left or right
            let primerpairID = seqID + "-" + direction + "-" + ppn;
            if (LR == "LEFT") primerpairs[primerpairID].left.hairpin = line.split("=")[1];
            else primerpairs[primerpairID].right.hairpin = line.split("=")[1];
        } else if (line.match("SELF_ANY_TH")) { // PRIMER_LEFT_0_SELF_ANY_TH=0.00
            let fields = line.split("_");
            let ppn = fields[2]; // primer pair number
            let LR = fields[1]; // left or right
            let primerpairID = seqID + "-" + direction + "-" + ppn;
            if (LR == "LEFT") primerpairs[primerpairID].left.anyTH = line.split("=")[1];
            else primerpairs[primerpairID].right.anyTH = line.split("=")[1];
        } else if (line.match("SELF_END_TH")) { // PRIMER_LEFT_0_SELF_END_TH=0.00
            let fields = line.split("_");
            let ppn = fields[2]; // primer pair number
            let LR = fields[1]; // left or right
            let primerpairID = seqID + "-" + direction + "-" + ppn;
            if (LR == "LEFT") primerpairs[primerpairID].left.endTH = line.split("=")[1];
            else primerpairs[primerpairID].right.endTH = line.split("=")[1];
        } else if (line.match("[0-9]+=")) { // start, length: PRIMER_LEFT_0=29,23
            // console.log(line);
            let fields = line.split("_");
            let ppn = fields[2].split("=")[0]; // primer pair number
            let LR = fields[1]; // left or right
            let primerpairID = seqID + "-" + direction + "-" + ppn;
            if (LR == "LEFT") {
                primerpairs[primerpairID].left.start = line.split("=")[1].split(",")[0];
                primerpairs[primerpairID].left.length = line.split("=")[1].split(",")[1];
            } else {
                primerpairs[primerpairID].right.start = line.split("=")[1].split(",")[0];
                primerpairs[primerpairID].right.length = line.split("=")[1].split(",")[1];
            }
        } else if (line.match("PRODUCT_SIZE=")) { // PRIMER_PAIR_0_PRODUCT_SIZE=78
            // console.log(line);
            let fields = line.split("_");
            let ppn = fields[2]; // primer pair number
            let primerpairID = seqID + "-" + direction + "-" + ppn;
            primerpairs[primerpairID].product_size = line.split("=")[1];
        } else if (line.match("COMPL_ANY_TH")) { // PRIMER_PAIR_0_COMPL_ANY_TH=0.00
            // console.log(line);
            let fields = line.split("_");
            let ppn = fields[2]; // primer pair number
            let primerpairID = seqID + "-" + direction + "-" + ppn;
            primerpairs[primerpairID].pairAnyTH = line.split("=")[1];
        } else if (line.match("COMPL_END_TH")) { // PRIMER_PAIR_0_COMPL_END_TH=0.00
            // console.log(line);
            let fields = line.split("_");
            let ppn = fields[2]; // primer pair number
            let primerpairID = seqID + "-" + direction + "-" + ppn;
            primerpairs[primerpairID].pairEndTH = line.split("=")[1];
        }
    }
    // print
    console.log("start printing");
    const FAM = "GAAGGTGACCAAGTTCATGCT";
    const VIC = "GAAGGTCGGAGTCAACGGATT";
    let out = "snpID,Direction,product_size,start,length,Tm,GC%,SELF_ANY_TH,SELF_END_TH,Hairpin,Seq,ReverseComplement,PAIR_COMPL_ANY_TH,PAIR_COMPL_END_TH\n";
    for (const [key, value] of Object.entries(primerpairs)) {
        // console.log(key, value);
        let leftSeq = value.left.seq;
        let leftSeqSNP1 = leftSeq.substring(0, leftSeq.length-1) + value.snp1.toUpperCase();
        let leftSeqSNP2 = leftSeq.substring(0, leftSeq.length-1) + value.snp2.toUpperCase();
        let leftSeqSNP1rc = reverse_complement(leftSeqSNP1);
        let leftSeqSNP2rc = reverse_complement(leftSeqSNP2);
        if (document.getElementById("addTail").checked) {
            leftSeqSNP1 = FAM + leftSeqSNP1;
            leftSeqSNP2 = VIC + leftSeqSNP2;
        }
        out += [key, "SNP1", value.product_size, value.left.start, value.left.length, value.left.tm, value.left.gc, value.left.anyTH, value.left.endTH, value.left.hairpin, leftSeqSNP1, leftSeqSNP1rc, value.pairAnyTH, value.pairEndTH].join(",") + "\n";
        out += [key, "SNP2", value.product_size, value.left.start, value.left.length, value.left.tm, value.left.gc, value.left.anyTH, value.left.endTH, value.left.hairpin, leftSeqSNP2, leftSeqSNP2rc].join(",") + "\n";
        out += [key, "Common", value.product_size, value.right.start, value.right.length, value.right.tm, value.right.gc, value.right.hairpin, value.right.anyTH, value.right.endTH, value.right.seq, reverse_complement(value.right.seq)].join(",") + "\n";
      }
      document.getElementById("stdout").innerHTML = out;
}

// function to format fasta to fixed length
async function csv2primer3 (fileContent) {
    let lines = fileContent.trim().split(/\r?\n/);

    let pickAnyway = "0";
    if (document.getElementById("PRIMER_PICK_ANYWAY").checked) pickAnyway = "1";
    let settings_common = "PRIMER_TASK=generic" + "\n" + 
                "PRIMER_PRODUCT_SIZE_RANGE=" + document.getElementById("PRIMER_PRODUCT_SIZE_RANGE").value + "\n" + 
                "PRIMER_MIN_SIZE=" + document.getElementById("PRIMER_MIN_SIZE").value + "\n" + 
                "PRIMER_OPT_SIZE=" + document.getElementById("PRIMER_OPT_SIZE").value + "\n" + 
                "PRIMER_MAX_SIZE=" + document.getElementById("PRIMER_MAX_SIZE").value + "\n" + 
                "PRIMER_MIN_TM=" + document.getElementById("PRIMER_MIN_TM").value + "\n" + 
                "PRIMER_OPT_TM=" + document.getElementById("PRIMER_OPT_TM").value + "\n" + 
                "PRIMER_MAX_TM=" + document.getElementById("PRIMER_MAX_TM").value + "\n" + 
                "PRIMER_PAIR_MAX_DIFF_TM=" + document.getElementById("PRIMER_PAIR_MAX_DIFF_TM").value + "\n" + 
                "PRIMER_MAX_HAIRPIN_TH=" + document.getElementById("PRIMER_MAX_HAIRPIN_TH").value + "\n" + 
                "PRIMER_FIRST_BASE_INDEX=1" + "\n" + 
                "PRIMER_LIBERAL_BASE=1" + "\n" + 
                "PRIMER_NUM_RETURN=" + document.getElementById("PRIMER_NUM_RETURN").value + "\n" + 
                "PRIMER_EXPLAIN_FLAG=1"  + "\n" + 
                "PRIMER_PICK_ANYWAY=" + pickAnyway + "\n";
    let newContent = "";
    for (let i = 0; i < lines.length; i++) {
        let line = lines[i].trim();
        if (line){
            let info = line.split(",");
            let snpID = info[0];
            // let chrom = info[1];
            let seq = info[2].toLowerCase(); // AAAAAAAAAAA[T/C]GGGGGGGGGGG
            let forwardSetting = await parseSNP(snpID, "forward", seq);
            let reverseSetting = await parseSNP(snpID, "reverse", reverse_complement(seq));
            // prepare primer3 input
            newContent += settings_common + forwardSetting;
            newContent += settings_common + reverseSetting;
        }

    }
    return newContent.trim();
}

// funciton to process each line
async function parseSNP(snpID, direction, seq){ // seq is AAAAAAAAAAA[T/C]GGGGGGGGGGG
    let [ll, snps, rr] = seq.split(/\[|\]/); // left flanking, [A/G], right flanking
    let [snp1, snp2] = snps.split("/");
    let template = ll + snp1 + rr;
    let snpPos = ll.length + 1; // default for SNPs
    let seqTarget = (snpPos+1).toString() + ",1"; // SEQUENCE_TARGET: <start>,<length>
    if (snp1.length != snp2.length) { // indels, look for the 1st difference to simplify the procedure
        let ss1 = (snp1 + rr).replace("-","");
        let ss2 = (snp2 + rr).replace("-","");
        let longSeq = ss1;
        let shortSeq = ss2;
        if (ss1.length < ss2.length) { // snp1 is shorter than snp2
            longSeq = ss2;
            shortSeq = ss1;
        }
        for (let j=0; j < shortSeq.length; j++){
            if (longSeq[j] != shortSeq[j]) {
                snpPos += j;
                snp1 = longSeq[j];
                snp2 = shortSeq[j];
                template = ll + longSeq;
                if (longSeq.length - shortSeq.length - j - 1 <= 0) seqTarget = (snpPos+1).toString() + ",1";
                else seqTarget = (snpPos+1).toString() + "," + (longSeq.length - shortSeq.length - j - 1).toString();
                break;
            }
        }
    }
    let settings = "SEQUENCE_ID=" + snpID + "__" + direction + "__" + snp1 + "__" + snp2 + "\n" +
    "SEQUENCE_TEMPLATE=" + template + "\n" + 
    "SEQUENCE_FORCE_LEFT_END=" + snpPos.toString() + "\n" + 
    "SEQUENCE_TARGET=" + seqTarget + "\n=\n";
    return settings;
}


// load the SNP file
document.getElementById("snpfile").addEventListener("change", prepareInput, false);


// download
async function download(){
    // document.getElementById("stdout").innerHTML = "Preparing downloading file ... It might take a while";
    let result = document.getElementById("stdout").innerHTML;
    let blob = new Blob([result], { type: "text/plain;charset=utf-8" });
    saveAs(blob, "KASP-primers.csv");
}

// IPUAC codes for conversion

var ipuac = {
    'A':'T',
    'C':'G',
    'G':'C',
    'T':'A',
    'a':'t',
    'c':'g',
    'g':'c',
    't':'a',
    // RNA
    'u':'a',
    'U':'A',
    //IPUAC: www.bioinformatics.org/sms/iupac.html
    "r":"y",
    "R":"Y",
    "y":"r",
    "Y":"R",
    "k":"m",
    "K":"M",
    "m":"k",
    "M":"K",
    "b":"v",
    "B":"V",
    "d":"h",
    "D":"H",
    "h":"d",
    "H":"D",
    "v":"b",
    "V":"B",
    "n":"n",
    "N":"N",
    "-":"-",
    "[":"]",
    "]":"[",
    "/":"/"
}
// to convert one sequence
function reverse_complement (seq, method = "RC"){
    // method: RC = reverse complement; C = just complement; R: just reverse
    var o = seq.split("");
    if (method == "R") return o.reverse().join('');
    // complement
    for (i=0; i < seq.length; i++){
        o[i] = ipuac[seq[i]];
    }

    if (method == "C") return o.join('');
    else return o.reverse().join('');
}