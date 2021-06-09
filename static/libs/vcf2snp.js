// This is just an example of the function below.
document.getElementById('start').onclick = function() {
    var file = document.getElementById('infile').files[0];
    if (!file) {
        // console.log('No file selected.');
        alert("No file selected!");
        return;
    }
    // var maxlines = parseInt(document.getElementById('maxlines').value, 10);
    // var lineno = 1;
    // readSomeLines is defined below.
    readSomeLines(file, parseLine, function onComplete() {
        // console.log('Read all lines');
        alert("Completed!");
    });
};

/**
 * Read up to and including |maxlines| lines from |file|.
 *
 * @param {Blob} file - The file to be read.
 * @param {integer} maxlines - The maximum number of lines to read.
 * @param {function(string)} forEachLine - Called for each line.
 * @param {function(error)} onComplete - Called when the end of the file
 *     is reached or when |maxlines| lines have been read.
 */
function readSomeLines(file, forEachLine, onComplete) {
    var CHUNK_SIZE = 64 * 1024; // 64 kb, arbitrarily chosen.
    var decoder = new TextDecoder();
    var offset = 0;
    var linecount = 0;
    var linenumber = 0;
    var results = '';
    var fr = new FileReader();
    fr.onload = function() {
        // Use stream:true in case we cut the file
        // in the middle of a multi-byte character
        results += decoder.decode(fr.result, {stream: true}); // put the begining of last line into next chunck
        // results += pako.inflate(fr.result, { to: 'string' });
        
        var lines = results.split(/\r?\n/);
        results = lines.pop(); // In case the line did not end yet. This is usually the beginning of hte last line
        linecount += lines.length;
        console.log(lines.length);
    
        for (var i = 0; i < lines.length; ++i) {
            let tt = forEachLine(lines[i]);
            document.getElementById('output').value += tt;
        }
        offset += CHUNK_SIZE;
        seek();
    };
    fr.onerror = function() {
        onComplete(fr.error);
    };
    seek();
    
    function seek() {
        if (offset !== 0 && offset >= file.size) {
            // We did not find all lines, but there are no more lines.
            if (results){
                let tt = forEachLine(results); // This is from lines.pop(), before.
                document.getElementById('output').value += tt;
            }
            onComplete(); // Done
            return;
        }
        var slice = file.slice(offset, offset + CHUNK_SIZE);
        fr.readAsArrayBuffer(slice);
    }
}

// parse each line for callback
function parseLine(line) {
    let ll = line.split("\t");
    if (line.startsWith('#')) {
        if (line.startsWith('#CHROM')){
            let tt = ll.slice(0,2).concat(ll.slice(3,6), ["AC", "AN", "DP", "MQ"], ll.slice(9));
            return tt.join("\t") + "\n";
        } else return "";
    }
    // get info
    let info = {}; // dict
    ll[7].split(";").forEach(
        m => {
          var [k, v] = m.split('=')
          info[k] = v
        }
    )
    let AC = info["AC"];
    let AN = info["AN"];
    let DP = info["DP"];
    let MQ = info["MQ"];
    // convert GT from number to alleles
    let GTs = ll.slice(9).map(x => {
        return x.split(":")[0];
    })
    let ref = ll[3];
    let alt = ll[4].split(","); // there may be more than one alternative alleles
    let alleles = [ref].concat(alt); // this way, 0 will be ref, and 1, 2 ... will be alternative allele
    let SNPs = GTs.map(x => {
        return gt2snp(alleles, x);
    });
    // out.write("\t".join(ll[0:2] + ll[3:6] + [AC, AN, DP, MQ] + SNPs) + "\n")
    let tt = ll.slice(0,2).concat(ll.slice(3,6), [AC, AN, DP, MQ], SNPs);
    return tt.join("\t") + "\n";
}

// convert GT to SNPs
function gt2snp(allele_list, gt) {
	let c = ""; // snp calls
	if (gt.includes(".")) // missing "./." or "."
		c = "N";
	else {
		let [a, b] = gt.split(/\/|\|/); // "/" unphased; "|" is phased
		if (a == b) // homozygous
			c = allele_list[a];
		else
			c = "H"; // heterozygous
    }
	return c;
}
// download the formatted fasta file
function download(){
    let snps = document.getElementById("output").value;
    let blob = new Blob([snps], { type: "text/plain;charset=utf-8" });
    saveAs(blob, "converted-snps.txt");
}