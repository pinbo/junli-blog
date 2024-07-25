const initSqlJs = window.initSqlJs;
const SQL = await initSqlJs({
    // Required to load the wasm binary asynchronously. Of course, you can host it wherever you want
    // You can omit locateFile completely when running in node
    locateFile: file => `/tools/sqljs/v1.10.3/sql-wasm.wasm`
  });

// const sqlPromise = initSqlJs({
// locateFile: file => `/tools/sqljs/v1.10.3/sql-wasm.wasm`
// });
// const dataPromise = fetch("https://dl.dropboxusercontent.com/scl/fi/f40srzyyp8mw5kn6zaqh6/Kronos_Homeolog_with_function2.db?rlkey=pyh3t76nroonxol1uueytekwx&st=mf42zky3").then(res => res.arrayBuffer());
// const [SQL, buf] = await Promise.all([sqlPromise, dataPromise])
// const db = new SQL.Database(new Uint8Array(buf));
// let sqlstr = "CREATE TABLE hello (gene char);";
// db.run(sqlstr);

// trim white spaces
function mytrim(x) {
    return x.replace(/^\s+|\s+$/g,'');
}


// Create an Object
const buffers = {}; //databases downloaded

// get homeolog
async function getHomeolog(){
    // document.getElementById("output").value = "";
    document.getElementById("alert").innerText = "Working on it... Be patient...";
    let inputString = document.getElementById("input").value.replace(/^\s+|\s+$/g,'');
    let lines = inputString.split('\n');
    let sel = document.getElementById("box1"); // box1 is selection of database
    let database = sel.options[sel.selectedIndex].value;
    let genePrefix = "";
    let dbURL = "";
    let buf;
    if (database == "Kronos_cDNA_v1.0") {
        genePrefix = "TrturKRN";
        dbURL = "https://jzseqbucket.s3.us-east-2.amazonaws.com/Kronos_Homeolog_with_function3.db"
        if ("Kronos" in buffers){
            buf = buffers.Kronos;
        } else {
            buf = await fetch(dbURL).then(res => res.arrayBuffer());
            buffers.Kronos = buf;
        }
    }
    else if (database == "CS_cDNA_HC_v1.1") {
        genePrefix = "TraesCS";
        dbURL = "https://jzseqbucket.s3.us-east-2.amazonaws.com/CS_Homeolog_with_function3.db";
        if ("CSv1HC" in buffers){
            buf = buffers.CSv1HC;
        } else {
            buf = await fetch(dbURL).then(res => res.arrayBuffer());
            buffers.CSv1HC = buf;
        }
    }
    // load db
    // const buf = await fetch(dbURL).then(res => res.arrayBuffer());
    const db = new SQL.Database(new Uint8Array(buf));
    let sqlstr = "CREATE TABLE hello (gene char);";
    db.run(sqlstr);
    // prepare input db
    for(let i = 0; i < lines.length; i++){
		let ll = mytrim(lines[i]);
        if (ll){
            let gene = ll.replace(genePrefix, "");
            console.log(gene);
            // db.run("INSERT INTO hello VALUES ('?')", gene);
            db.run("INSERT INTO hello VALUES ('" + gene + "')");
	    }
    }
    sqlstr = "Select \
  d1.name As c1, \
  d2.name As c2, \
  t.pctIdent As c3, \
  d3.name As c4, \
  d2.At_ident As c5, \
  d3.description As c6, \
  d4.name As c7, \
  d2.Os_ident As c8, \
  d4.description As c9 \
From wheat_hit t \
Join wheatID_and_hits d1 On ( d1.gene = t.gene ) \
Join wheatID_and_hits d2 On ( d2.gene = t.hit ) \
Left Join AtOsID d3 On ( d3.gene = d2.AtID ) \
Left Join AtOsID d4 On ( d4.gene = d2.OsID ) \
Join hello ON hello.gene = d1.name \
ORDER BY d1.name;";
    // const stmt = db.prepare("SELECT * FROM Kronos_Homeolog WHERE gene = $gene");
    const stmt = db.prepare(sqlstr);
//   stmt.getAsObject({$gene:gene});
    // const row = stmt.getAsObject({});
    // console.log('Here is a row: ' + JSON.stringify(row));
    document.getElementById('tbody').innerHTML = "";
    while(stmt.step()) { //
        const row = stmt.getAsObject();
        console.log('Here is a row: ' + JSON.stringify(row));
        let text = "";
        let k = '<tr>';
        row.c1 = genePrefix + row.c1;
        row.c2 = genePrefix + row.c2;
        for (let x in row) {
            text += row[x] + "\t";
            k+= '<td>' + row[x] + '</td>';
        }
        k+= '</tr>';
        // document.getElementById("output").value += text + "\n";
        document.getElementById('tbody').innerHTML += k;
    }
    stmt.free(); // free the memory used by the statement
    document.getElementById("alert").innerText = "Done!";

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
	// document.getElementById("output").value = "";
    document.getElementById('tbody').innerHTML = "";
    document.getElementById("alert").innerText = "";
};

document.getElementById("run").addEventListener("click", getHomeolog);
document.getElementById("clearseq").addEventListener("click", clearseq);
document.getElementById("example").addEventListener("click", giveExamples);
document.getElementById("copytable").addEventListener("click", copy2clipboard);
