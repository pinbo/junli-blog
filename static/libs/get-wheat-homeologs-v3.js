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
const dataPromise1 = fetch("https://jzseqbucket.s3.us-east-2.amazonaws.com/Kronos_Homeolog_with_function4.db.gz").then(res => res.arrayBuffer()).then(raw => pako.inflate(raw));
const dataPromise2 = fetch("https://jzseqbucket.s3.us-east-2.amazonaws.com/CS_Homeolog_with_function4.db.gz").then(res => res.arrayBuffer()).then(raw => pako.inflate(raw));
const [SQL, buf1, buf2] = await Promise.all([sqlPromise, dataPromise1, dataPromise2])
// const db = new SQL.Database(new Uint8Array(buf));
dbs.Kronos = new SQL.Database(new Uint8Array(buf1));
dbs.CSv1HC = new SQL.Database(new Uint8Array(buf2));

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
    // load db
    // const buf = await fetch(dbURL).then(res => res.arrayBuffer());
    // const db = new SQL.Database(new Uint8Array(buf));
    // let sqlstr = "CREATE TABLE hello (name char);";
    // db.run(sqlstr);
    // prepare input db
    let whereStr = "";
    for(let i = 0; i < lines.length; i++){
		let ll = mytrim(lines[i]);
        if (ll){
            let gene = ll.replace(genePrefix, "");
            if (gene == ll){
                document.getElementById("alert").innerText = "Wrong gene names found; please check your gene names or SELECT proper database!";
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
WHERE d1.name IN (" + whereStr + ") \
ORDER BY d1.name;";
    // const stmt = db.prepare("SELECT * FROM Kronos_Homeolog WHERE gene = $gene");
    // console.log(sqlstr);
    const stmt = db.prepare(sqlstr);
    // stmt.bind({$gene:whereStr});
    // const row = stmt.getAsObject({});
    // console.log('Here is a row: ' + JSON.stringify(row));
    document.getElementById('tbody').innerHTML = "";
    while(stmt.step()) { //
        const row = stmt.getAsObject();
        // console.log('Here is a row: ' + JSON.stringify(row));
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
	// document.getElementById("output").value = "";
    document.getElementById('tbody').innerHTML = "";
    document.getElementById("alert").innerText = "";
};

document.getElementById("run").addEventListener("click", getHomeolog);
document.getElementById("clearseq").addEventListener("click", clearseq);
document.getElementById("example").addEventListener("click", giveExamples);
document.getElementById("copytable").addEventListener("click", copy2clipboard);
