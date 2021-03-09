function formatFasta(evt) {
	var files = evt.target.files; // FileList object
	if(files.length==0) return;
	var pre = document.getElementById("formatted");
	// while(pre.hasChildNodes()) pre.removeChild( pre.firstChild );
    pre.innerText = "";
	
    var reader = new FileReader();
    reader.onprogress=function(evt) {}
    reader.onloadend=function(e) {
        if(e.target.result==null) return;
        var lines=e.target.result.split(/\r?\n/);
        var dna="";
        var title="";
        var i=0;
        var line;
        let width = parseInt(document.getElementById("line_width").value);
        for(;;){
            if(i==lines.length || (line=lines[i].trim())[0]=='>'){// trim to avoid space at the beginning
                if(dna.length != 0){
                    pre.appendChild(document.createTextNode(title+"\n"));
                    while(dna.length != 0){
                        var n=Math.min(dna.length, width);
                        pre.appendChild(document.createTextNode(dna.substring(0,n)+"\n"));
                        dna= dna.substring(n);
                    }
                }
                if(i===lines.length) break;
                title=line.replace(/> +/, '>');// in case something like "> seq1"
                dna="";
            }
            else dna+=line.trim().replace(/ +/g,""); // in case space inside
            ++i;
        }
    };
    reader.readAsText(files[0]);
}

// download the formatted fasta file
function download(){
    fasta = document.getElementById("formatted").innerText;
    let blob = new Blob([fasta], { type: "text/plain;charset=utf-8" });
    saveAs(blob, "formatted_fasta.fa");
}

document.getElementById('fasta').addEventListener('change', formatFasta, false);