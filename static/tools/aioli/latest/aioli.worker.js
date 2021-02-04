// =============================================================================
// Config
// =============================================================================

// State
DEBUG = false;
MSG_UUID = null;

// Stdout/stderr indexed by message uuid
STDOUT = {};
STDERR = {};

// Files mounted and paths
FILES = [];
DIR_DATA_FILES = "/data";
DIR_DATA_URLS = "/urls";

// JZ: function to remove elemnt by value from an arry
function arrayRemove(arr, value) { 
    return arr.filter(function(ele){ 
        return ele != value; 
    });
}
// Initialization -- two conditions for this worker to be ready:
//   1) Got UUID from Main Thread that it sent with the "init" message
//   2) Wasm module is initialized
resolveInitWasm = null;
resolveInitWorker = null;
promiseInitWasm = new Promise(resolve => resolveInitWasm = resolve);
promiseInitWorker = new Promise(resolve => resolveInitWorker = resolve);
Promise.all([ promiseInitWasm, promiseInitWorker ])
       .then(() => send(MSG_UUID, "ready"));

// WebAssembly Module config
Module = {
    // When the module is initialized, resolve the initWasm promise
    onRuntimeInitialized: () => {
        // Setup folders
        FS.mkdir(DIR_DATA_FILES);
        FS.mkdir(DIR_DATA_URLS);
        // Resolve promise
        resolveInitWasm();
    },

    // Load .wasm/.data files from a custom path
    locateFile: (path, dir) => {
        var dirRoot = "";

        // Use hardcoded path if `BIOWASM_URL` was defined when creating WebWorker script
        if(typeof BIOWASM_URL !== 'undefined')
            dirRoot = BIOWASM_URL;
        // Or infer it from the path to the JS file
        else {
            var dirJS = self.location.href;
            dirRoot = dirJS.substring(0, dirJS.lastIndexOf("/") + 1);
        }
        return dirRoot + path;    
    },

    // Setup print functions to store stdout/stderr based on id
     print: text => STDOUT[MSG_UUID] += `${text}\n`,
     printErr: text => STDERR[MSG_UUID] += `${text}\n`
}


// =============================================================================
// Handle messages from the outside
// =============================================================================

// Format: "message": d => { <do stuff>; return <response>; }
API = {
    // -------------------------------------------------------------------------
    // Initialize WebWorker: resolve initWorker promise when get UUID from main thread
    // -------------------------------------------------------------------------
    init: (id, data) => {
        resolveInitWorker();
    },

    // -------------------------------------------------------------------------
    // File system operations
    // -------------------------------------------------------------------------
    ls: (id, path) => {
        return FS.readdir(path);
    },

    cat: (id, path) => {
        return FS.readFile(path, { encoding: "utf8" });
    },

    downloadText: (id, path) => {
        let file = FS.readFile(path, { encoding: "utf8" });
        let blob = new Blob([ file ]);
        // return URL.createObjectURL(blob);
        return blob;
    },
    downloadBinary: (id, path) => {
        let file = FS.readFile(path);//, { encoding: "utf8" });
        let blob = new Blob([ file ]);
        // return URL.createObjectURL(blob);
        return blob;
    },

    fs: (id, config) => {
        let fn = config.fn;
        let args = config.args;
        
        try {
            if(!(fn in FS))
                throw `Invalid function ${fn}. See <https://emscripten.org/docs/api_reference/Filesystem-API.html> for valid functions.`;
            let response = FS[fn](...args);
            if(response == null)
                response = "ok";
            return response;    
        } catch(err) {
            console.error(`[AioliWorker] Failed to run FS.${fn}(${args}): ${err}`);
            return "error";
        }
    },

    // added by junli
    symlink: (id, oldpath) => {
    // symlink: (id, path) => {
        let newpath = oldpath.replace("/data/", "/bwa2/examples/");
        FS.symlink(oldpath, newpath);
        return newpath;
        // return FS.readdir(path);
    },
    // added by junli
    write: (id, file) => {
        return FS.writeFile(file.name, file.content);
        // const buffer = file.content;
        // const stream = FS.open(file.name, "w+");
        // FS.write(stream, buffer, 0, buffer.length, 0);
        // FS.close(stream);
    },
    // added by junli: remove a file
    rm: (id, path) => {
        return FS.unlink(path);
    },
    // added by junli: set working diretory 
    setwd: (id, path) => {
        return FS.chdir(path);
    },

    // -------------------------------------------------------------------------
    // Call main function with custom command
    // -------------------------------------------------------------------------
    exec: (id, command) => {
        // Initialize stdout/stderr
        STDOUT[id] = "";
        STDERR[id] = "";

        // Call main function with command
        Module.callMain(command.split(" "));

        // Re-open stdout/stderr (fix error "error closing standard output: -1")
        FS.streams[1] = FS.open("/dev/stdout", "w");
        FS.streams[2] = FS.open("/dev/stderr", "w");

        return {
            stdout: STDOUT[id],
            stderr: STDERR[id]
        };
    },

    // -------------------------------------------------------------------------
    // Mount files
    // -------------------------------------------------------------------------
    mount: (id, file) => {
        // Support File objects
        // console.log(file.name);
        // console.log(file.source);
        // console.log(file.file);
        // console.log(file.file instanceof File);
        // console.log(file.file instanceof Blob);
        // file.path = DIR_DATA_FILES + "/" + file.name;
        console.log("File path " + file.path);
        if(file.source == "file")
        {
            // Unmount & remount all files (can only mount a folder once)
            // try {
            //     FS.unmount(DIR_DATA_FILES);
            // } catch(e) {}
            // FILES.push(file);
            // // Handle File and Blob objects
            // FS.mount(WORKERFS, {
            //     files: FILES.filter(f => f.file instanceof File).map(f => f.file),
            //     blobs: FILES.filter(f => f.file instanceof Blob).map(f => ({ name: f.name, data: f.file }))
            // }, DIR_DATA_FILES);

            // junli added for adding a file in the default filesystem
            let reader= new FileReader();
            reader.addEventListener('loadend', function(){
                let result=reader.result;
                const uint8_view = new Uint8Array(result);
                FS.writeFile(file.path, uint8_view);
            });
            reader.readAsArrayBuffer(file.file);
        }

        // Support URLs
        else if(file.source == "url")
            FS.createLazyFile(DIR_DATA_URLS, file.name, file.url, true, true);

        // Otherwise invalid input
        else throw "Only accept File objects or URL strings.";
        // FS.symlink(oldpath, newpath): create symlink node
        // FS.symlink(file.path, "/bwa2/example/" + file.name);

        return file.path;
    },

    // -------------------------------------------------------------------------
    // Transfer files from one worker to another. This is useful when one Worker
    // creates a file on their file system (i.e. a file that is not mounted) and
    // another Worker needs access to that file.
    // -------------------------------------------------------------------------
    transfer: (id, data) => {
        const role = data.role;
        const port = data.port;
        const path = data.path;

        // If this is the WebWorker that is sending the file to the other worker,
        // first read the file and then *transfer* (not copy!) the ArrayBuffer over
        if(role == "sender") {
            const file = FS.readFile(path);
            port.postMessage(file, [file.buffer]);
            FS.unlink(path); // delete the file to save memory
        }

        // If this is the WebWorker receiving files, write the ArrayBuffer to a file
        else if(role == "receiver") {
            port.onmessage = d => {
                const buffer = d.data;
                const stream = FS.open(path, "w+");
                FS.write(stream, buffer, 0, buffer.length, 0);
                FS.close(stream);
            }
        }
        return 0;
    },
};

// -------------------------------------------------------------------------
// On message handler
// -------------------------------------------------------------------------
onmessage = message => {
    // Parse message
    const id = message.data.id;
    const data = message.data.data;
    const action = message.data.action || "default";
    log(`MainThread Says: Action=%c${action}%c; Data=%c${JSON.stringify(data)}%c [id=${id}]`, "color:red; font-weight:bold", "", "color:red; font-weight:bold");

    // Figure out what to do and return in response
    MSG_UUID = id;
    const response = API[action](id, data);
    if(response != null)
        send(id, response);
    
    // Clean up stdout/stderr after sending message
    delete STDOUT[id];
    delete STDERR[id];
}


// =============================================================================
// Utility functions
// =============================================================================

function send(id, response, action="resolve", transferables=[])
{
    log(`Sending: Action=%c${action}%c; Data=%c${JSON.stringify(response)}%c [id=${id}]`, "color:deepskyblue; font-weight:bold", "", "color:deepskyblue; font-weight:bold");
    postMessage({
        id: id,
        action: action,
        data: response
    }, transferables);
}

function log(message)
{
    if(!DEBUG)
        return;

    // Get all arguments except `message`
    let args = [...arguments];
    args.shift();
    console.warn(`%c[Worker]%c ${message}`, "font-weight:bold", "", ...args);
}
