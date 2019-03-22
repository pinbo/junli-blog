// trim white spaces
function mytrim(x) {
    return x.replace(/^\s+|\s+$/g,'');
}

// read fasta lines pasted in textarea with id "input"
function readfasta(){
	var lines = document.getElementById("input").value.replace(/^\s+|\s+$/g,'').split('\n');
	var seqdict = {};
	var seqnames = [];
	var mykey = "";
	for(var i = 0; i < lines.length; i++){
		//code here using lines[i] which will give you each line
		var ll = mytrim(lines[i]);
		if (ll[0] == ">") {
			mykey = ll.replace(/^[> ]+/g,"");
			seqdict[mykey] = "";
			seqnames.push(mykey);
		}
		else if (mykey in seqdict) seqdict[mykey] += ll;
		else seqdict[mykey] = ll;
	}
	return [seqdict, seqnames];
}

// concatenate tow list element by element
function myzip(list1, list2){//two lists should have equal length
	if (list1.length !== list2.length) return;
	var list3 = [];
	for (var i=0; i<list1.length; i++){
		list3.push(list1[i] + list2[i]); //string + string
	}
	return list3;
}

// merge continuous deletions: ATTT to A---
function allpos(arr, val) {
    var indexes = [], i;
    for(i = 0; i < arr.length; i++)
        if (arr[i] === val)
            indexes.push(i);
    return indexes.join(",");
}

// clear input box
function clearseq() {
	document.getElementById("input").value = "";
	document.getElementById("output").value = "";
};

// example input
function paste_example(){
	document.getElementById("input").value = ">seq1\nAATTG-GGGCCT\n\n>seq2\nAAATGAGC--CC\n\n>seq3\nAAATG-GCGCCC\n";
}

function callsnps(){
	let [fasta, seqnames] = readfasta();
	var refid = seqnames[0]; // reference name
	var refseq = fasta[refid]; // ref seq
	var refseq_nogap = refseq.replace(/-/g,"");

	var n = -1;
	var seq2 = {}; // new dictionary of rotated sequences with no-gap postions as keys

	for (var i = 0; i < refseq.length; i++){ // i is the position of seq position
		if (refseq[i] !== "-") {n += 1};
		var templist = [];
		for (var j=0; j < seqnames.length; j++){ // j is seq name
			var name = seqnames[j];
			var ss = fasta[name];
			templist.push(ss[i]);
		}
		if (!(n in seq2)) {
			seq2[n] = templist;
		}
		else {
			seq2[n] = myzip(seq2[n], templist);
		}
	}
	// output
	var outlist = [];

	for (var k=0; k < refseq_nogap.length; k++){
		// remove "-" in alleles in "A--"
		var seq = seq2[k].map(function(x){
			if (x.length > 1) return x.replace(/-/g,"");
			else return x;
		});

		//alleles = set(seq) // set: unique values
		let alleles = new Set(seq);
		if (alleles.size !== 1){// if there is variation
			var ref_allele = refseq_nogap[k]; // string
			//var alt_allele_set = alleles - set(ref_allele) // set
			let alt_allele_set = [...alleles].filter(x => x!==ref_allele); // an array without ref_alllele
			var alt_allele = alt_allele_set.join('/'); // string
			outlist.push([k+1, ref_allele, alt_allele].concat(seq)); // K+1: starting from 1
		}
	}

	var outlist2 = [["Pos", "Ref", "alt"].concat(seqnames)] // merged continuous deletions
	var mm = []; // temperary merging of lists
	var n = 0; // the start position of continious deletion
	for (var i = 0; i < outlist.length - 1; i++){
		var p1 = outlist[i];
		var p2 = outlist[i+1];
		var allpos1 = allpos(p1, "-");
		var allpos2 = allpos(p2, "-");
		if (p1[0] + 1 === p2[0] && p1.indexOf("-") >= 0 && p2.indexOf("-") >= 0 && allpos(p1, "-") === allpos(p2, "-") ){
			if (mm.length > 0) // if mm already have merged lists
				mm = myzip(mm, p1);
			else {
				mm = p1;
				n = p1[0];
			}
		}
		else {
			if (mm.length > 0){
				mm = myzip(mm, p1);
				mm[0] = n;
				var mm2 = mm.map(x => x.toString().indexOf("-") >= 0 ? "-" : x); // replace "-----" to "-"
				outlist2.push(mm2 ); 
				mm = [];
			}
			else outlist2.push(p1);
			if (i + 1 === outlist.length - 1)
				outlist2.push(p2);
		}
	}

// output is the id for the output box
	var output = [];
	for (var i =0; i < outlist2.length; i++){
		var vv = outlist2[i];
		output.push(vv.join("\t"));
	}
// output is the id for the output box
	document.getElementById("output").value = output.join("\n");
}




