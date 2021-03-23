document.getElementById('fasta').addEventListener('change', formatFastaFile, false);

async function formatFastaFile(evt) {
	var files = evt.target.files; // FileList object
	if(files.length==0) return;
	var pre = document.getElementById("formatted");
	// while(pre.hasChildNodes()) pre.removeChild( pre.firstChild );
    // pre.innerText = "";
    let fileContent = await readTextFileAsync(files[0]);
    let width = parseInt(document.getElementById("line_width").value);
    pre.innerText = formatFasta(fileContent, width);
}

function formatFastaPaste() {
    let fileContent = document.getElementById("paste").value;
    let width = parseInt(document.getElementById("line_width").value);
    var pre = document.getElementById("formatted");
    pre.innerText = formatFasta(fileContent, width);
}

// download the formatted fasta file
function download(){
    fasta = document.getElementById("formatted").innerText;
    let blob = new Blob([fasta], { type: "text/plain;charset=utf-8" });
    saveAs(blob, "formatted_fasta.fa");
}


// clear input sequences
function clearseq() {
    document.getElementById("paste").value = "";
    document.getElementById("formatted").innerText = "";
};

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