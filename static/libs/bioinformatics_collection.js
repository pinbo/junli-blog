// function to format fasta to fixed length
function formatFasta (fileContent, width = 60) {// 60 bp per line
    let lines = fileContent.split(/\r?\n/);
    let dna = "";   // all the bases of 1 seq (removed new lines)
    let title = ""; // title of seq
    let line = "";  // content of each line
    let i = 0;      // line number
    let newContent = "";
    for (;;){
        if(i==lines.length || (line=lines[i].trim())[0]=='>'){// trim to avoid space at the beginning
            if(dna.length != 0){
                newContent += title + "\n";
                while(dna.length != 0){
                    let n = Math.min(dna.length, width);
                    newContent += dna.substring(0,n) + "\n";
                    dna = dna.substring(n);
                }
            }
            if(i===lines.length) break;
            title = line.replace(/> +/, '>');// in case something like "> seq1"
            dna = "";
        } else dna += line.trim().replace(/ +/g, ""); // in case space inside
        ++i;
    }
    return newContent;
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

// function to read a text file async, for loadRef
function readTextFileAsync(file) {
    return new Promise((resolve, reject) => {
        let reader = new FileReader();
        if (file.name.split('.').pop() == "gz") {
            reader.onload  = () => {
                // resolve(pako.inflate(new Uint8Array(reader.result), { to: 'string' }));
                resolve(pako.inflate(reader.result, { to: 'string' }));
            }
        } else {
            reader.onload = () => {
                resolve(new TextDecoder('utf-8').decode(reader.result, {stream: true}));
            };
        }
        reader.onerror = reject;
        reader.readAsArrayBuffer(file);
    })
}