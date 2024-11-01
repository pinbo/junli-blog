// v6: based on v3, but use left join to get the same order as input even with duplicates

const initSqlJs = window.initSqlJs;
// const SQL = await initSqlJs({
//     // Required to load the wasm binary asynchronously. Of course, you can host it wherever you want
//     // You can omit locateFile completely when running in node
//     locateFile: file => `/tools/sqljs/v1.10.3/sql-wasm.wasm`
//   });

  // Create an Object
const dbs = {}; //databases downloaded

const sqlPromise = initSqlJs({
    locateFile: file => `/tools/sqljs/v1.10.3/sql-wasm.wasm`
});
const dataPromise1 = fetch("https://jzseqbucket.s3.us-east-2.amazonaws.com/Kronos_Homeolog_with_function6.db.gz").then(res => res.arrayBuffer()).then(raw => pako.inflate(raw));
const dataPromise2 = fetch("https://jzseqbucket.s3.us-east-2.amazonaws.com/CS_Homeolog_with_function6.db.gz").then(res => res.arrayBuffer()).then(raw => pako.inflate(raw));
const [SQL, buf1, buf2] = await Promise.all([sqlPromise, dataPromise1, dataPromise2])
// const db = new SQL.Database(new Uint8Array(buf));
dbs.Kronos = new SQL.Database(new Uint8Array(buf1));
dbs.CSv1HC = new SQL.Database(new Uint8Array(buf2));
dbs.Kronos.run("CREATE TABLE hello (name TEXT);");
dbs.CSv1HC.run("CREATE TABLE hello (name TEXT);");

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
    let genePrefix = "";
    let db;
    if (database == "Kronos_cDNA_v1.0") {
        genePrefix = "TrturKRN";
        db = dbs.Kronos;
    }
    else if (database == "CS_cDNA_HC_v1.1") {
        genePrefix = "TraesCS";
        db = dbs.CSv1HC;
    }
    // prepare input db
    // db.run("DROP TABLE IF EXISTS hello;");
    // db.run("CREATE TABLE hello (name TEXT);");
    // prepare input db
    let sqlstr = "INSERT INTO hello (name) VALUES";
    let whereStr = "";
    for(let i = 0; i < lines.length; i++){
		let ll = mytrim(lines[i]);
        if (ll){
            let gene = ll.replace(genePrefix, "");
            // console.log(gene);
            // db.run("INSERT INTO hello VALUES ('" + gene + "')");
            sqlstr += " ('" + gene + "'),";
            whereStr += "'" + gene + "',";
	    }
    }
    sqlstr = sqlstr.replace(/,$/, ';');
    whereStr = whereStr.replace(/,$/, '');
    console.log(sqlstr);
    db.run(sqlstr);
    
    sqlstr = "Select \
  d1.name As c1, \
  d2.name As c2, \
  t.pctIdent As c3, \
  d3.name As c4, \
  d2.At_ident As c5, \
  d2.At_alignment, \
  d3.description As c6, \
  d4.name As c7, \
  d2.Os_ident As c8, \
  d2.Os_alignment, \
  d4.description As c9 \
From wheatID_and_hits d1 \
join wheat_hit t on ( d1.gene = t.gene ) \
Join wheatID_and_hits d2 On ( d2.gene = t.hit ) \
Left Join AtOsID d3 On ( d3.gene = d2.AtID ) \
Left Join AtOsID d4 On ( d4.gene = d2.OsID ) \
WHERE d1.name IN (" + whereStr + ") \
ORDER BY d1.name;";
    // const stmt = db.prepare("SELECT * FROM Kronos_Homeolog WHERE gene = $gene");
    // console.log(sqlstr);
    document.getElementById("thead").innerHTML = "<th>WheatGeneID</th> \
        <th>Best Wheat matches</th> \
        <th>Wheat %identity</th> \
	    <th>Best At matches</th> \
	    <th>At %identity</th> \
	    <th>At align length</th> \
	    <th>At description</th> \
	    <th>Best Os matches</th> \
	    <th>Os %identity</th> \
	    <th>Os align length</th> \
	    <th>Os description</th>";
    if (document.getElementById("check1").checked){// only output At/Os hits
        sqlstr = "Select \
hello.name As c1, \
d3.name As c4, \
d2.At_ident As c5, \
d2.At_alignment, \
d3.description As c6, \
d4.name As c7, \
d2.Os_ident As c8, \
d2.Os_alignment, \
d4.description As c9 \
From hello \
Left Join wheatID_and_hits d2 on (hello.name = d2.name) \
Left Join AtOsID d3 On ( d3.gene = d2.AtID ) \
Left Join AtOsID d4 On ( d4.gene = d2.OsID );";
        // also change table title
        document.getElementById("thead").innerHTML = "<th>WheatGeneID</th> \
<th>Best At matches</th> \
<th>At %identity</th> \
<th>At align length</th> \
<th>At description</th> \
<th>Best Os matches</th> \
<th>Os %identity</th> \
<th>Os align length</th> \
<th>Os description</th>";
    }
    if (document.getElementById("check2").checked){ // over-write sqlstr from check1 
        sqlstr = "Select \
d2.name As c1, \
d3.name As c4, \
d2.At_ident As c5, \
d2.At_alignment, \
d3.description As c6, \
d4.name As c7, \
d2.Os_ident As c8, \
d2.Os_alignment, \
d4.description As c9 \
From wheatID_and_hits d2 \
Join AtOsID d3 On ( d3.gene = d2.AtID ) \
Join AtOsID d4 On ( d4.gene = d2.OsID ) \
WHERE c4 in (" + whereStr + ") \
OR c7 in (" + whereStr + ");";
    // also change table title
    document.getElementById("thead").innerHTML = "<th>WheatGeneID</th> \
<th>Best At matches</th> \
<th>At %identity</th> \
<th>At align length</th> \
<th>At description</th> \
<th>Best Os matches</th> \
<th>Os %identity</th> \
<th>Os align length</th> \
<th>Os description</th>";
    }
    const stmt = db.prepare(sqlstr);
    document.getElementById('tbody').innerHTML = "";
    const newTable = document.getElementById("tbody");
    while(stmt.step()) { //
        const row = stmt.getAsObject();
        row.c1 = genePrefix + row.c1;
        if ('c2' in row) row.c2 = genePrefix + row.c2;
        const newRow = document.createElement("tr");
        for (let x in row) {
            const newtd = document.createElement("td");
            newtd.textContent = row[x];
            newRow.appendChild(newtd);
        }
        newTable.appendChild(newRow);
    }
    stmt.free(); // free the memory used by the statement
    document.getElementById("alert").innerText = "Done!";
    // db.run("DROP TABLE IF EXISTS hello;");
    db.run("DELETE FROM hello;");
    return 0;
}


// give examples
function giveExamples(){
    let sel = document.getElementById("box1"); // box1 is selection of database
    let database = sel.options[sel.selectedIndex].value;
    if (database == "Kronos_cDNA_v1.0") document.getElementById("input").value = "TrturKRN1A01G056670\nTrturKRN5B01G032410";
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
