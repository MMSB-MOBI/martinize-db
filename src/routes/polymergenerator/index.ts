import { Router } from 'express';
import glob from 'glob';
import ItpFile from 'itp-parser';
import SocketIo from 'socket.io';
import ShellManager, { JobInputs } from '../../Builders/ShellManager';
import TmpDirHelper from '../../TmpDirHelper';
import * as fs from 'fs';
import { POLYPLYPATHDATA, POLYPLY_VENV } from "../../constants";
import checkError from './errorParser';
import { Readable } from 'stream';
//import {jobFS} from 'ms-jobmanager';

interface Blob {
    readonly size: number;
    readonly type: string;
    arrayBuffer(): Promise<ArrayBuffer>;
    slice(start?: number, end?: number, contentType?: string): Blob;
    stream(): NodeJS.ReadableStream;
    text(): Promise<string>;
}


const polymer = Router();

polymer.get('/data', async (req, res) => {
    let avaibleData: any = {}
    let listfile = glob.sync(POLYPLYPATHDATA + "/polyplydata/*/*.+(itp|ff)").map(f => { return f })
    // console.log(listfile)

    for (let file of listfile) {
        let forcefield = file.split('/')[file.split('/').length - 2]

        const itp = await ItpFile.read(file);
        for (let e of itp.getField('moleculetype')) {
            if (!e.startsWith(';')) {
                const mol = e.split(' ')[0]
                if (Object.keys(avaibleData).includes(forcefield)) {
                    avaibleData[forcefield].push(mol)
                }
                else {
                    avaibleData[forcefield] = [mol]
                }
            }
        }
    }

    // for (let ff of Object.keys(avaibleData)) {
    //     let list2 = glob.sync(PATHDATA + ff + "/*.itp").map(f => { return f })
    //     for (let file of list2) {
    //         let forcefield = file.split('/')[6]

    //         const itp = await ItpFile.read(file);
    //         for (let e of itp.getField('moleculetype')) {
    //             e = e.split('\t')[0]
    //             if (!e.startsWith(';')) {
    //                 const mol = e.split(' ')[0]
    //                 if (Object.keys(avaibleData).includes(forcefield)) {
    //                     avaibleData[forcefield].push(mol)
    //                 }
    //                 else {
    //                     avaibleData[forcefield] = [mol]
    //                 }
    //             }
    //         }
    //     }
    // }

    //remove duplkiicate 
    for (let forcefield of Object.keys(avaibleData)) {
        avaibleData[forcefield] = [...new Set(avaibleData[forcefield])];
    }

    res.send(avaibleData);
    console.log("Sending forcefields and residues data")
});

// polymer.get("/customdata/:forcefield", (req, res) => {
//     let avaibleData: any = {}
//     let testdico: any = {}
//     console.log(req.params)
//     let ff = req.params.forcefield
//     const pathPolyplyData =
//         glob(PATHDATA + ff + "/*.itp", async function (er, files) {
//             for (let file of files) {
//                 let forcefield = file.split('/')[7]
//                 const itp = await ItpFile.read(file);
//                 for (let e of itp.getField('moleculetype')) {
//                     e = e.split('\t')[0]
//                     if (!e.startsWith(';')) {
//                         const mol = e.split(' ')[0]
//                         if (Object.keys(avaibleData).includes(forcefield)) {
//                             avaibleData[forcefield].push(mol)
//                         }
//                         else {
//                             avaibleData[forcefield] = [mol]
//                         }
//                         if (Object.keys(testdico).includes(forcefield)) {
//                             testdico[forcefield].push(mol)
//                         }
//                         else {
//                             testdico[forcefield] = [mol]
//                         }
//                     }
//                 }
//             }
//             res.send(avaibleData);
//         })
// });

polymer.get("/fastaconversion", (req, res) => {
    const fasta = {
        'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
        'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
        'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
        'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'
    }
    res.send(fasta);
})

export async function SocketIoPolymerizer(socket: SocketIo.Socket) {
    const WORKDIR = "/data3/rmarin/projet_polyply/job"
    console.log("je suis dans SocketIoPolymerizer")

    // JobManager.start({ 'port': 6001, 'TCPip': "localhost" })
    //     .then(() => {
    //         console.log("JobManager start")


    socket.on("runpolyply", async (dataFromClient: any) => {

        const tmp_dir = await TmpDirHelper.get();
        console.log("Run polyply gen itp in ", tmp_dir)

        //Get forcefield 
        const ff = dataFromClient['polymer']['forcefield']

        const name = dataFromClient['name']

        const density = dataFromClient['density']

        let additionalfile = ""
        if (dataFromClient['customITP'] !== undefined) {
            additionalfile = dataFromClient['customITP'].join("\n")
        }

        // Stream 
        const json_stream: Readable = new Readable();
        json_stream.push(JSON.stringify(dataFromClient.polymer));
        json_stream.push(null)

        const exportVar = {
            polyplyenv: POLYPLY_VENV,
            ff: ff,
            density: density,
            name: name,
            action: "itp",
            customfile: additionalfile,
        }

        const inputs = {
            "martiniForceField.itp": POLYPLYPATHDATA + "/martini_v3.0.0.itp",
            "polymer.json": json_stream,
        }

        //Need to add module quand exectuer sur child instead of jobmanager
        // },
        //     "modules": ["polyply"]
        //         }

        let result: string = ""
        try {
            console.log(ShellManager.mode)
            if (ShellManager.mode === "child") {
                console.log("### Running polyply with CHILD")
                await ShellManager.run('polyply', { exportVar, inputs/* , "modules": ["polyply"]*/ }, tmp_dir, "create_itp", undefined);
                result = fs.readFileSync(tmp_dir + "/create_itp.stdout").toString();
            }
            else {
                console.log("### Running polyply")
                const { jobFS, stdout } = await ShellManager.run('polyply', { exportVar, inputs /* , "modules": ["polyply"]*/ }, tmp_dir, "create_itp", undefined);

                result = stdout
            }
        }
        catch (e) {
            console.log(e)
            // Handle error and throw the right error
            console.error("ShellManager.run crash");
        }


        const itp = result.split("STOP\n")[0]
        const error = result.split("STOP\n")[1]

        console.log("error", error)

        const errorParsed = checkError(error)
        if (errorParsed.ok == true) {
            //Then on fait la requete pour le gro
            console.log("yes on passe au gro")
            socket.emit("itp", itp)
            socket.on("continue", async () => {
                //Liste les fichier itp du champs de force 
                //Tres belle ligne (1h de travail)
                //let bordelitp = glob.sync(PATHDATA + ff + "/*.itp").map(f => { return "#include " + f + "\r" }).join('')
                //super stupid variable pour avoir un retour a la ligne
                const topfilestr = `#include "${POLYPLYPATHDATA + "/martini_v3.0.0.itp"}"
         #include "polymere.itp"
         [ system ]
         ; name
         mylovelypolymer
         [ molecules ]
         ; name  number
         ${name} 1`

                //const topFile = tmp_dir + "/system.top"
                //fs.writeFileSync(topFile, topfilestr)


                // Stream itp 
                const itp_stream: Readable = new Readable();
                itp_stream.push(itp);
                itp_stream.push(null)


                const top_stream: Readable = new Readable();
                top_stream.push(topfilestr);
                top_stream.push(null)


                const exportVar = {
                    polyplyenv: POLYPLY_VENV,
                    density: density,
                    name: name,
                    action: "gro"
                }

                const inputs = {
                    "itp": itp_stream,
                    "martiniForceField": POLYPLYPATHDATA + "/martini_v3.0.0.itp",
                    "top": top_stream
                }

                try {
                    console.log(ShellManager.mode)
                    if (ShellManager.mode === "child") {
                        console.log("### Running polyply with CHILD")
                        await ShellManager.run('polyply', { exportVar, inputs /* , "modules": ["polyply"]*/ }, tmp_dir, "create_gro", undefined);
                        result = fs.readFileSync(tmp_dir + "/create_itp.stdout").toString();
                        socket.emit("gro", result)
                    }
                    else {
                        console.log("### Running polyply")
                        const { jobFS, stdout } = await ShellManager.run('polyply', { exportVar, inputs /* , "modules": ["polyply"]*/ }, tmp_dir, "create_gro", undefined);
                        socket.emit("gro", stdout);
                    }
                    
                }
                catch (e) {
                    console.log( 'ERROR WITH GRO')
                    console.log(e)
                    // Handle error and throw the right error
                    console.error("ShellManager.run crash");
                }


            })
        }
        else {
            console.log("oups", errorParsed)
            socket.emit("oups", errorParsed)
        }

    })
}

export default polymer;
