const initSqlJs = window.initSqlJs;
// const SQL = await initSqlJs({
//     // Required to load the wasm binary asynchronously. Of course, you can host it wherever you want
//     // You can omit locateFile completely when running in node
//     locateFile: file => `/tools/sqljs/v1.10.3/sql-wasm.wasm`
//   });

const sqlPromise = initSqlJs({
locateFile: file => `/tools/sqljs/v1.10.3/sql-wasm.wasm`
});
const dataPromise = fetch("https://jzseqbucket.s3.us-east-2.amazonaws.com/cDNA_Kronos_v1_vs_CS_HC_LC_v2.db.gz").then(res => res.arrayBuffer()).then(raw => pako.inflate(raw));
const [SQL, buf] = await Promise.all([sqlPromise, dataPromise])
const db = new SQL.Database(new Uint8Array(buf));


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
    let sel = document.getElementById("box1"); // box1 is selection of todo
    let todo = sel.options[sel.selectedIndex].value;
    // let genePrefix = "TrturKRN";
    let wheat = "Kronos";
    if (todo == "C2K") {
        // genePrefix = "TraesCS";
        wheat = "CS";
    }
    // prepare input cmd
    let whereStr = "";
    for(let i = 0; i < lines.length; i++){
		let ll = mytrim(lines[i]);
        if (ll){
            let gene = ll; //ll.replace(genePrefix, "");
            console.log(gene);
            whereStr += "'" + gene + "',";
	    }
    }
    whereStr = whereStr.replace(/,$/, '');
    console.log(whereStr);
    let sqlstr = "SELECT * FROM hits WHERE hits." + wheat + " IN (" + whereStr + ") ORDER BY " + wheat + ";";
    console.log(sqlstr);
    const stmt = db.prepare(sqlstr);
    // stmt.bind({$gene:whereStr});
    // const row = stmt.getAsObject({});
    // console.log('Here is a row: ' + JSON.stringify(row));
    document.getElementById('tbody').innerHTML = "";
    const newTable = document.getElementById("tbody");
    while(stmt.step()) { //
        const row = stmt.getAsObject();
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

}


// give examples
function giveExamples(){
    let sel = document.getElementById("box1"); // box1 is selection of todo
    let todo = sel.options[sel.selectedIndex].value;
    if (todo == "K2C") document.getElementById("input").value = "TrturKRN6B01G025800\nTrturKRN7A01G081900";
    else if (todo == "C2K") document.getElementById("input").value = "TraesCS5A02G391700\nTraesCS5B02G396600";
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
