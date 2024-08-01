// v5: from v3, example to use sqlite3 website wasm package
// need to download wasm packages from sqlite3 officical website

// function to import database from URL
const sqlite3 = await globalThis.sqlite3InitModule({
    /* We can redirect any stdout/stderr from the module like so, but
        note that doing so makes use of Emscripten-isms, not
        well-defined sqlite APIs. */
    // print: log,
    // printErr: error
  });

const makedb = function(buf){
    const db = new sqlite3.oo1.DB();
    const bytes = new Uint8Array(buf);
    const p = sqlite3.wasm.allocFromTypedArray(bytes);
    db.onclose = {after: function(){sqlite3.wasm.dealloc(p)}};
    const rc = sqlite3.capi.sqlite3_deserialize(
        db.pointer, 'main', p, bytes.length, bytes.length,
        0
    );
    db.checkRc(rc);
    return db;
}

const dbs = {}; //databases downloaded

const dataPromise1 = fetch("https://jzseqbucket.s3.us-east-2.amazonaws.com/Kronos_Homeolog_with_function4.db.gz").then(res => res.arrayBuffer()).then(raw => pako.inflate(raw));
const dataPromise2 = fetch("https://jzseqbucket.s3.us-east-2.amazonaws.com/CS_Homeolog_with_function4.db.gz").then(res => res.arrayBuffer()).then(raw => pako.inflate(raw));
const [buf1, buf2] = await Promise.all([dataPromise1, dataPromise2])
dbs.Kronos = buf1;
dbs.CSv1HC = buf2;

// trim white spaces
function mytrim(x) {
    return x.replace(/^\s+|\s+$/g,'');
}

// get homeolog
async function getHomeolog(){
    // document.getElementById("output").value = "";
    document.getElementById("alert").innerText = "Working on it... Be patient...";
    let inputString = document.getElementById("input").value.replace(/^\s+|\s+$/g,'');
    let lines = inputString.split('\n');
    let sel = document.getElementById("box1"); // box1 is selection of database
    let database = sel.options[sel.selectedIndex].value;
    let genePrefix = "TrturKRN";
    let db = makedb(dbs.Kronos);
    if (database == "CS_cDNA_HC_v1.1") {
        genePrefix = "TraesCS";
        db = makedb(dbs.CSv1HC);
    }
    // load db
    let whereStr = "";
    for(let i = 0; i < lines.length; i++){
		let ll = mytrim(lines[i]);
        if (ll){
            let gene = ll.replace(genePrefix, "");
            if (gene == ll){
                document.getElementById("alert").innerText = ll + ": wrong gene names found; please check your gene names or SELECT proper database!";
                document.getElementById("alert").style.color = "red";
                return 1;
            }
            whereStr += "'" + gene + "',";
	    }
    }
    whereStr = whereStr.replace(/,$/, '');
    // console.log(whereStr);
    let sqlstr = "Select \
  d1.name As c1, \
  d2.name As c2, \
  t.pctIdent As c3, \
  d3.name As c4, \
  d2.At_ident As c5, \
  d3.description As c6, \
  d4.name As c7, \
  d2.Os_ident As c8, \
  d4.description As c9 \
From wheatID_and_hits d1 \
join wheat_hit t on ( d1.gene = t.gene ) \
Join wheatID_and_hits d2 On ( d2.gene = t.hit ) \
Left Join AtOsID d3 On ( d3.gene = d2.AtID ) \
Left Join AtOsID d4 On ( d4.gene = d2.OsID ) \
WHERE c1 IN (" + whereStr + ") \
ORDER BY c1;";

    let resultRows = [];
    try {
        db.exec({
        sql: sqlstr,
        rowMode: 'object',
        resultRows: resultRows
        });
    }catch (error) {
        console.log("Junli error: ", error);
    } finally{
        console.log("SQL running is DONE!");
        db.close();
    }
    console.log("Result rows:",JSON.stringify(resultRows,undefined,2));

    document.getElementById('tbody').innerHTML = "";
    // inner.HTML is much much slower than appendChild
    // for (let row of resultRows) { //
    //     let k = '<tr>';
    //     row.c1 = genePrefix + row.c1;
    //     row.c2 = genePrefix + row.c2;
    //     for (let x in row) {
    //         k+= '<td>' + row[x] + '</td>';
    //     }
    //     k+= '</tr>';
    //     document.getElementById('tbody').innerHTML += k;
    // }
    // const newTable = document.createElement("table");
    // newTable.innerHTML = "<thead><th>Player</th><th>Score</th></thead>";
    const newTable = document.getElementById("tbody");
    for(let row of resultRows){
        row.c1 = genePrefix + row.c1;
        row.c2 = genePrefix + row.c2;
        const newRow = document.createElement("tr");
        for (let x in row) {
            const newtd = document.createElement("td");
            newtd.textContent = row[x];
            newRow.appendChild(newtd);
        }
        newTable.appendChild(newRow);
    }
    document.getElementById("alert").innerText = "Done!";
    return 0;
}


// give examples
function giveExamples(){
    let sel = document.getElementById("box1"); // box1 is selection of database
    let database = sel.options[sel.selectedIndex].value;
    if (database == "Kronos_cDNA_v1.0") document.getElementById("input").value = "TrturKRN6B01G025800\nTrturKRN7A01G081900";
    else if (database == "CS_cDNA_HC_v1.1") document.getElementById("input").value = "TraesCS5A02G391700\nTraesCS5B02G396600";
}

function copy2clipboard() {
    // Get the text field
    let tableElement = document.getElementById("datatable");
  
    // Convert the table element to an HTML string
    const tableHTMLString = tableElement.outerHTML;

    // Convert the table element to plain text string
    const plainTextString = tableElement.innerText;

    navigator.clipboard.write([
        new ClipboardItem({
            'text/html': new Blob([tableHTMLString], {
                type: 'text/html',
            }),
            'text/plain': new Blob([plainTextString], {
                type: 'text/plain',
            }),
        }),
    ]);
}

// clear input box
function clearseq() {
	document.getElementById("input").value = "";
    document.getElementById('tbody').innerHTML = "";
    document.getElementById("alert").innerText = "";
};

document.getElementById("run").addEventListener("click", getHomeolog);
document.getElementById("clearseq").addEventListener("click", clearseq);
document.getElementById("example").addEventListener("click", giveExamples);
document.getElementById("copytable").addEventListener("click", copy2clipboard);
