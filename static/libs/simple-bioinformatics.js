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

// trim white spaces
function mytrim(x) {
    return x.replace(/^\s+|\s+$/g,'');
}

// read fasta lines pasted in textarea with id "input"
function readfasta(mode){
	var lines = document.getElementById("input").value.replace(/^\s+|\s+$/g,'').split('\n');
	if (mode == "single_line") return lines; // a list
	var seqdict = {};
	var mykey = ">unnamed";
	for(var i = 0; i < lines.length; i++){
		//code here using lines[i] which will give you each line
		var ll = mytrim(lines[i]);
		if (ll[0] == ">") {mykey = ll; seqdict[mykey] = ""}
		else if (mykey in seqdict) seqdict[mykey] += ll;
		else seqdict[mykey] = ll;
	}
	return seqdict;
}


// overall function to get desired conversion
function revcomp(){
    var sel = document.getElementById("box1"); // box1 is the selection ID
    var sel2 = document.getElementById("box2"); // mode: fasta or single line: single line will treat each line as an input unit
	var method = sel.options[sel.selectedIndex].value;
	var mode = sel2.options[sel2.selectedIndex].value;
	//var seq = document.getElementById("input").value.replace(/\s/g,""); // in case multiple lines
	//var outseq = reverse_complement(seq, method);
	var seqs = readfasta(mode);
	var outseq = "";
	if (mode == "single_line"){// seqs is a list
		for(var i = 0; i < seqs.length; i++){
			//code here using lines[i] which will give you each line
			var ll = mytrim(seqs[i]);
			var rc = reverse_complement(ll, method);
			outseq += rc + "\n";
		}
	}
	else {//fasta mode, seqs is an dict
		for (var i in seqs){
			var rc = reverse_complement(seqs[i], method);
			// if users paste in one sequece wihtout name
			if (i == ">unnamed") outseq += rc + "\n\n";
			else outseq += i + " " + method + "\n" + rc + "\n\n";
		}
    }
    // output is the id for the output box
	document.getElementById("output").value = outseq;
	};

// clear input box
function clearseq() {
	document.getElementById("input").value = "";
	document.getElementById("output").value = "";
};
