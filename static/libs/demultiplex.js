// global variables
var dictBarcode = {}; // a dictionary of barcodes
var dictSample = {}; // a dictionary of demultiplexed samples: key is barcodes
//myStorage = window.sessionStorage;


//download a zipped file containing all gz files
function download(){
	// document.getElementById("progress2").max = Object.keys(dictSample).length;
	// document.getElementById("progress2").value = 0;
	// document.getElementById("gzip-progress").style.visibility = "visible";
	document.getElementById("demo3").innerHTML = "... Creating files for download. It may take a while depending on your input size";
	var zip = new JSZip();
	for (let k in dictSample){
		// document.getElementById("progress2").value += 1;
		let text = dictSample[k];
		let output = pako.gzip(text);
		zip.file(k + ".gz", output);
	}
	zip.generateAsync({type:"blob"})
		.then(function(content) {
			// see FileSaver.js
			saveAs(content, "decomplexed-samples.zip");
		});
}

// clear input sequences
function clearseq() {
    document.getElementById("left").value = "";
    document.getElementById("right").value = "";
	document.getElementById("demo1").innerHTML = "";
	document.getElementById("demo2").innerHTML = "";
    document.getElementById('output').value = "";
    document.getElementById("download-btn").style.visibility = "hidden";
	// document.getElementById("demultiplex-progress").style.visibility = "hidden";
	// document.getElementById("gzip-progress").style.visibility = "hidden";
};

//add example input
function putExample(){
    document.getElementById("left").value = "TCCTCTGTCACGGAAGCG";
    document.getElementById("right").value = "TTTAGCCTCCCCACCGAC";
}

//
function readFile(file, onLoadCallback){
    var fr = new FileReader();
	if (file.name.split('.').pop() == "gz") {
		//console.log("unzipping gz");
		fr.onload = function(event) {
			var result = pako.inflate(event.target.result, { to: 'string' });
			onLoadCallback(result);
		}
	} else {
		fr.onload = function(event) {
			var result = new TextDecoder('utf-8').decode(event.target.result, {stream: true});
			onLoadCallback(result);
		}
	}
	// let blob = f.slice(0, 200); // slice not working with file list
	// fr.readAsArrayBuffer(blob);
	fr.readAsArrayBuffer(file);
}

// get barcode from the barcode file
// 3 columns: sample name, left barcode, right barcode
function getBarcode(fileContent){
	var lines = fileContent.split(/\r?\n/);
	//var barcodes = {}; // a dictionary of barcodes
	for (var line of lines){
		if (line){
			var ss = line.split(/\t/);
			dictBarcode[ss[1] + "-" + ss[2]] = ss[0];
		}
	}
	//return barcodes;
	document.getElementById('demo1').innerHTML = "Finished reading barcodes"
}

// read barcode file
function readBarcode(){
	var files = document.getElementById('barcode').files;
    if (!files.length) {
      alert('Please select a barcode file!');
      return;
    }
	var f = files[0];
	readFile(f, getBarcode)
}

// read fastq file
function readFastq(){
	var files = document.getElementById('fastq').files;
    if (!files.length) {
      alert('Please select a fastq or fastq.gz file!');
      return;
    }
	var f = files[0];
	readFile(f, demultiplex)
}

// function process
function startAnalyze(){
	readBarcode();
	readFastq();
}
// sample class
class sample {
	constructor() {
		this.R1 = []; // 4 lines of read1
		this.R2 = []; // 4 lines of read2
	}
}

// process fastq file
function demultiplex(fileContent){
	var lines = fileContent.split(/\r?\n/);
	// document.getElementById("progress1").max = lines.length;
	// document.getElementById("progress1").value = 0;
	// document.getElementById("demultiplex-progress").style.visibility = "visible";
	//var samples = {}; // a dictionary of demultiplexed samples
	var leftAdapter = document.getElementById('left').value;
	var rightAdapter = document.getElementById('right').value;
	var leftBarcodeLen = document.getElementById('leftBarcodeLen').value;
	var rightBarcodeLen = document.getElementById('rightBarcodeLen').value;
	//var barcodeLenToCheck = 8;// whole barcode
	var dictPairReads = {};
	//var dictPairReads = {} # dictionary of sample reads R1 and R2
	var n1 = 0; 
	var n2 = 1;
	var readID = "";
	for (var line of lines){
		// document.getElementById("progress1").value += 1;
		if (line){
			let ss = line.split(/\t/);
			if (line.startsWith("@") && line.includes(" ")){ // I found sometimes quality line (the 4th line) also starts with @, but they have no space
			// example: @M02850:171:000000000-J54GV:1:1101:19874:1325 1:N:0:1
				readID = line.split(" ")[0]; // read is R1 or R2, here 1:N:0:1
				n1 = 2; // counter of the 4 lines of each read
				if (readID in dictPairReads){ // should be R2
					dictPairReads[readID].R2.push(line);
					n2 = 2; // now reads are R2
				} else {
					ss = new sample();
					ss.R1.push (line);
					dictPairReads[readID] = ss;
				}
			} else if (n2 == 1){// input to R1
				dictPairReads[readID].R1.push(line);
				n1 += 1;
			} else if (n2 == 2 && n1 < 4){
				dictPairReads[readID].R2.push(line);
				n1 += 1;
			} else if (n2 == 2 && n1 == 4){// # ready to write
				dictPairReads[readID].R2.push(line);
				n2 = 1; // reset to read 1
				var leftbarcode = "";
				var rightbarcode = "";
				var R1First30bp = dictPairReads[readID].R1[1].substring(0,30); // the first 30 bp of R1
				var R2First30bp = dictPairReads[readID].R2[1].substring(0,30); // the first 30 bp of R2
				var P1 = R1First30bp.indexOf(leftAdapter); // position of the left adpator in R1
				var P2 = R2First30bp.indexOf(leftAdapter);
				var P3 = R1First30bp.indexOf(rightAdapter);
				var P4 = R2First30bp.indexOf(rightAdapter);
				if (P1 >= 0 && P4 >= 0){// # if there are still at least 5 bp on the left
					leftbarcode = R1First30bp.substring(P1-leftBarcodeLen, P1);
					rightbarcode = R2First30bp.substring(P4-rightBarcodeLen,P4);
				}
				if (P2 >= 0 && P3 >= 0){
					leftbarcode = R2First30bp.substring(P2-leftBarcodeLen, P2);
					rightbarcode = R1First30bp.substring(P3-rightBarcodeLen, P3);
					// switch R1 and R2, so R1 always has the left adapter
					let tmp = dictPairReads[readID].R1;
					dictPairReads[readID].R1 = dictPairReads[readID].R2;
					dictPairReads[readID].R2 = tmp;
				}
				barcode = leftbarcode + "-" + rightbarcode;
				if (barcode in dictBarcode){
					out1 = dictBarcode[barcode] + "_R1_001.fastq";
					out2 = dictBarcode[barcode] + "_R2_001.fastq";
					dictSample[out1] = (dictSample[out1] || "") + dictPairReads[readID].R1.join("\n") + "\n";
					dictSample[out2] = (dictSample[out2] || "") + dictPairReads[readID].R2.join("\n") + "\n";
				} else {
					out = "unassigned.fastq";
					dictSample[out] = (dictSample[out] || "") + (dictPairReads[readID].R1.concat(dictPairReads[readID].R2)).join("\n") + "\n";
				}
				delete dictPairReads[readID];
			}
		}
	}
	document.getElementById('demo2').innerHTML = "Finished demultiplexing";
	document.getElementById("download-btn").style.visibility = "visible";
}