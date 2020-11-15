// global variables
// var seqdict = {};

// function to generate results
async function analyze(){
    document.getElementById('output').value = "";
    document.getElementById('template').value = "";
    await analyzeFiles();
    //document.getElementById('list').value = '<table id="summary">' + createTable(document.getElementById('list').value.split('\n')) + '</table>';
    //document.getElementById("demo").innerHTML = createTable(document.getElementById('list').value.split('\n'));
    if (document.getElementById('output').value){
        document.getElementById("download").style.visibility = "visible";
    }
}

// trim white spaces
function mytrim(x) {
    return x.replace(/^\s+|\s+$/g,'');
}

// read any text file
function readFile(file, onLoadCallback){
    var fr = new FileReader();
	fr.onload = function(event) {
		var result = new TextDecoder('utf-8').decode(event.target.result, {stream: true});
		onLoadCallback(result);
	}
	fr.readAsArrayBuffer(file);
}

// read fasta lines
function getfasta(fileContent){
    var lines = fileContent.split(/\r?\n/);
    var seqdict = {};
    var mykey = ">unnamed";
    for(var i = 0; i < lines.length; i++){
        //code here using lines[i] which will give you each line
        var ll = mytrim(lines[i]);
        if (ll){
            if (ll[0] == ">") {mykey = ll.split(/ +|>/)[1]; seqdict[mykey] = ""}
            else if (mykey in seqdict) seqdict[mykey] += ll;
            else seqdict[mykey] = ll;
        }
    }
    return seqdict;
}
// read fasta file
function readfasta(geneID){
	var files = document.getElementById('reference').files;
    if (!files.length) {
      alert('Please select your reference file!');
      return;
    }
	var f = files[0];
    // readFile(f, getfasta);
    var fr = new FileReader();
    fr.geneID = geneID;
	fr.onload = function (event) {
        var result = new TextDecoder('utf-8').decode(event.target.result, { stream: true });
        var seqdict = getfasta(result);
        if (!(event.target.geneID in seqdict))
            alert('Your sequence name is NOT in the fasta file!');
        else
            document.getElementById('template').value = seqdict[event.target.geneID];
    }
	fr.readAsArrayBuffer(f);
	// return seqdict;
}

// to get all the input sequences
async function getSeq() {
	// get the input sequences
    // var seqdict = readfasta();
	// console.log(seqdict);
	var sequences = document.getElementById('sequences').value.trim().split(/ +/);
	var geneID = sequences[0];
    console.log(geneID);
    readfasta(geneID)
    var leftSeq = sequences[1].toUpperCase();
    var rightSeq = sequences[2].toUpperCase();
	var gRNA = sequences[3].toUpperCase();
	// if (!(geneID in seqdict)) {
    //     alert('Your sequence name is NOT in the fasta file!');
    //     return;
    // }
    while(!document.getElementById('template').value) {
        await new Promise(r => setTimeout(r, 10));
    }
    var template = document.getElementById('template').value.toUpperCase();
    console.log("template:\n" + template);
    // template = template.replace(/(\r\n|\n|\r| )/gm,""); // wild type template
    // check whether to read on R1 or R2
    var sel = document.getElementById("box1"); // box1 is the selection ID
    var strand = sel.options[sel.selectedIndex].value;
    if (strand == "R"){
        console.log("Use Reverse Complement Sequences");
        var leftSeq0 = leftSeq;
        leftSeq = reverseComplement(rightSeq);
        rightSeq = reverseComplement(leftSeq0);
        gRNA = reverseComplement(gRNA);
        template = reverseComplement(template);
    }
    // get the intact sequences
    var n1 = template.indexOf(leftSeq); // start position of left flanking sequence
    var n2 = template.indexOf(rightSeq); // start position of the right flanking sequence
    if (n1 < 0){
        alert('Left sequence is NOT in the template!');
        return;
    }
    if (n2 < 0){
        alert('Right sequence is NOT in the template!');
        return;
    }
    var wtSeq = template.substring(n1, n2 + rightSeq.length);
    // check whether we get the correct wild type seqeunce
    //document.getElementById("demo").innerHTML = "<strong>Intact sequence</strong><br>" + wtSeq;
    return [leftSeq, rightSeq, gRNA, wtSeq];
}
// add click to start analyzing
// https://www.html5rocks.com/en/tutorials/file/dndfiles//#toc-reading-files

// function to analyze the fastq file
function checkIndel(leftSeq, rightSeq, gRNA, wtSeq, filename, fileContent){//read sequence is the R1 or reverse complement R2 reads
    var lines = fileContent.split(/\r?\n/);
    var matchLeft = 0;
    var matchRight = 0;
    var matchGRNA = 0;
    var matchBoth = 0;
    var indels = {}; // a dictionary for indel count: key is length difference; value is the counts
    var uniqMutSeq = {}; // a dictionary for unique mutSeq count: key is mutSeq; value is the counts
    var lineNum = 0; // line number count
    var nindel = 0;
    for (var readseq of lines){
        lineNum += 1;
        if (lineNum%4 != 2) {
            continue
        } // only check the sequence line
        var n1 = readseq.indexOf(leftSeq); // start position of left flanking sequence
        var n2 = readseq.indexOf(rightSeq);
        if (n1 >= 0){
            matchLeft += 1;
        }
        if (n2 >= 0){
            matchRight += 1;
        }
        if (n1 >= 0 && n2 >= 0){
            matchBoth += 1;
            if (readseq.includes(gRNA)){
                matchGRNA += 1;
            }
            var mutSeq = readseq.substring(n1, n2 + rightSeq.length);
            var indelLen = mutSeq.length - wtSeq.length;
            if (indelLen != 0){
                nindel += 1;
            }
            uniqMutSeq[mutSeq] = (uniqMutSeq[mutSeq] || 0) + 1; // if exist +=1, else 0+1=1
            indels[indelLen] = (indels[indelLen] || 0) + 1; // if exist +=1, else 0+1=1
        }
    }
    // check linenumber to make sure everything is read
    console.log(filename + " has line number " + (lineNum - 1));
    document.getElementById("demo").innerHTML += filename + " has line number " + (lineNum - 1) + "<br>";
    //return [matchLeft, matchRight, matchGRNA, matchBoth, nindel, indels, uniqMutSeq];
    pamPos = wtSeq.indexOf(gRNA) + gRNA.length; // PAM position
    console.log("matchBoth  " + matchBoth);
    if (matchBoth) {
        document.getElementById('output').value += [filename, matchBoth, matchGRNA, matchGRNA/matchBoth*100, nindel, nindel/matchBoth*100, matchLeft, matchRight, matchLeft/matchRight].join(",");
        // find the top2 mutations
        var sortedUniqSeq = sortDict(uniqMutSeq);
        var cc = 0
        for (var mutSeq of sortedUniqSeq){
            if (cc > 1) break;
            if (mutSeq[0] != wtSeq){
                cc += 1
                console.log(wtSeq)
                console.log(mutSeq[0])
                const [ref, alt, difLeft, difRight] = align(wtSeq, mutSeq[0]);
                document.getElementById('output').value += "," + [mutSeq[0].length - wtSeq.length, mutSeq[1], mutSeq[1]/matchBoth*100, mutSeq[0], ref, alt, pamPos - difLeft].join(",");
            }
        }
        document.getElementById('output').value += "\n";
    }
    document.getElementById("progress").value += 1;
}

// process multiple fastq files
//function handleFileSelect(evt){
async function analyzeFiles() {
    //var files = evt.target.files; // FileList object
    // files is a FileList of File objects. List some properties.
    var files = document.getElementById('files').files;
    document.getElementById("progress").max = files.length;
    document.getElementById("progress").value = 0;
    document.getElementById("progress").style.visibility = "visible";
    if (!files.length) {
        alert('Please select a file!');
        return;
    }
    // get the sequences
    const [leftSeq, rightSeq, gRNA, wtSeq] = await getSeq();
    document.getElementById('output').value += "Intact reference sequence," + wtSeq + "\n\n";
    document.getElementById('output').value += "fastq_file,number_of_matched_reads,number_of_reads_with_intact_gRNA,%intact_gRNA,total_indel,%indel,number_of_reads_with_leftSeq,number_of_reads_with_rightSeq,nleftSeq/nrightSeq,";
    document.getElementById('output').value += "#1_indel,#1_count,#1_%,#1_seq,#1_ref,#1_alt,#1_bp_left_of_PAM,";
    document.getElementById('output').value += "#2_indel,#2_count,#2_%,#2_seq,#2_ref,#2_alt,#2_bp_left_of_PAM\n";
    for(var i = 0, f; f = files[i]; i++) {
        // read file content
        console.log(f.name);
        //document.getElementById("demo").innerHTML += "Processing " + f.name + "<br>";
        var fr = new FileReader();
        fr.filename = f.name;
        if (f.name.split('.').pop() == "gz") {
            console.log("unzipping gz");
            fr.onload = function(event) {
                var result = pako.inflate(event.target.result, { to: 'string' });
                checkIndel(leftSeq, rightSeq, gRNA, wtSeq, event.target.filename,result);
            }
        } else {
            fr.onload = function(event) {
                var result = new TextDecoder('utf-8').decode(event.target.result, {stream: true});
                checkIndel(leftSeq, rightSeq, gRNA, wtSeq, event.target.filename,result);
            }
        }
        // let blob = f.slice(0, 200); // slice not working with file list
        // fr.readAsArrayBuffer(blob);
        fr.readAsArrayBuffer(f);
    } //end of file loop
    //document.getElementById("demo").innerHTML += "Job done: all input files have been processed!<br>";
    //console.log(document.getElementById('output').value);
    //var table1 = document.getElementById('output').value.split(/\n/);
    //document.getElementById('output').value = createTable(document.getElementById('output').value.split('\n'));
}
//document.getElementById('files').addEventListener('change', handleFileSelect, false);


// function to sort a dictionary by desending value, and return an array of [key, value]
function sortDict(dict){
    // Create items array
    var items = Object.keys(dict).map(function(key) {
        return [key, dict[key]];
    });
    // Sort the array based on the second element
    items.sort(function(first, second) {
        return second[1] - first[1];
    });
    return items
}
// function to check indel position
function align(wt, mut){//two sequences
    var shortSeq = (wt.length < mut.length ? wt : mut);
    var difLeft = 0; // the first difference from the left
    var difRight = 0; //// the first difference from the right
    for (let i = 0; i < wt.length; i++) {
        if (wt[i] != mut[i]){
            difLeft = i;
            break;
        }
    }
    wt2 = reverseString(wt);
    mut2 = reverseString(mut);
    for (let i = 0; i < wt.length; i++) {
        if (wt2[i] != mut2[i]){
            difRight = i;
            break;
        }
    }
    if (difLeft + difRight > shortSeq.length) difLeft = shortSeq.length - difRight;
    var ref = wt.substring(difLeft-1, wt.length-difRight);
    var alt = mut.substring(difLeft-1, mut.length-difRight);
    return [ref, alt, difLeft, difRight];
}
//align("ACGTGGTGCGCCT", "ACGTGGCCT")
// reverse a string
function reverseString(str) {
    return str.split("").reverse().join("");
}

// clear input sequences
function clearseq() {
    document.getElementById("sequences").value = "";
    document.getElementById("demo").innerHTML = "";
    document.getElementById('output').value = "";
    document.getElementById('template').value = "";
    document.getElementById("download").style.visibility = "hidden";
    document.getElementById("progress").style.visibility = "hidden";
};

// to convert one sequence
function reverseComplement (seq, method = "RC"){
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
    "N":"N"
}

// download file
function download(){
    var text = document.getElementById('output').value;
    var filename = "";
    if (document.getElementById('download-name').value){
        filename = document.getElementById('download-name').value;
    } else filename = "CRISPR-eidting-summary.csv";
    var blob = new Blob([text], {
        type: "text/plain;charset=utf-8"
    });
    saveAs(blob, filename);
}



//add example input
function putExample(){
    document.getElementById("sequences").value = "DELLA-4B-Kronos GAGGAGGTGGACGAGC GCGCCGCGCCCGACG GCTGGCGGCGCTCGGGTACA";
}