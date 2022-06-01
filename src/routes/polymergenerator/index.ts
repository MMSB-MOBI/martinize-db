import { Router } from 'express';
import glob from 'glob';
import ItpFile from 'itp-parser';
import SocketIo from 'socket.io';
//import * as jobmanagerClient from 'ms-jobmanager';
// import { PromiseManager } from '/data3/rmarin/projet_polyply/msjob-aspromise';
import * as JobManager from 'ms-jobmanager';
import ShellManager, { JobInputs } from '../../Builders/ShellManager';
import TmpDirHelper from '../../TmpDirHelper';
import * as fs from 'fs';
import {POLYPLYPATHDATA} from "../../constants";

interface ErrorToClient {
    disjoint: boolean,
    errorlinks: any[]
}
 
//const PATHDATA = "/data3/rmarin/projet_polyply/PolymerGeneratorServerDev/data/"
const polymer = Router();

//Build a dictionnary with all molecule avaible and a dictionnary with linking rules
// Settings file "settings.json" at project root
//const jobmanagerClient = new PromiseManager("localhost", 6001)


polymer.get('/hello', async (req, res) => {
    console.log("hello")
    res.json("hello")
})

polymer.get('/data', async (req, res) => {
    let avaibleData: any = {}
    let listfile = glob.sync(POLYPLYPATHDATA + "/polyply_1.0/polyply/data/*/*.+(itp|ff)").map(f => { return f })
     
    for (let file of listfile) {
        let forcefield = file.split('/')[ file.split('/').length - 2 ]

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
    console.log("Sending forcefields and residues data" )
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

const PORT = process.env.PORT || 4123;

const server = require('http').createServer(polymer)
const io = SocketIo(server, { path: '/socket' });

// server.listen(PORT, () => {
//     console.log(`Server listening on ${PORT}`);
// });



function checkResult(result: string): boolean {
    let restemp = result.split("STOP\n")
    //process error str
    let processerror = restemp[1].split('\n').filter((line: string) => line.startsWith('WARNING - general ')).filter((line: string) => !line.startsWith('WARNING - general - Node '))
    console.log(processerror)
    if (processerror.length === 0) {
        return false
    }
    else {
        return true
    }
}

function parseError(result: string) {
    let processerror = result.split('STOP')[1].split('\n').filter((line: string) => line.startsWith('WARNING - general ')).filter((line: string) => !line.startsWith('WARNING - general - Node '))

    //init le dico d'erreur potentiel
    let dicErreur: ErrorToClient = { disjoint: false, errorlinks: [] }
    if (processerror.filter((line: string) => line.startsWith('WARNING - general - Your molecule consists of disjoint parts.Perhaps links were not applied correctly.'))) {
        dicErreur.disjoint = true
    }
    let wronglinkmerror = processerror.filter((line: string) => line.startsWith("WARNING - general - Missing link between"))
    for (let i of wronglinkmerror) {
        let splitline = i.split(' ')
        let resname1 = splitline[9]
        let idname1 = parseInt(splitline[8]) - 1
        let resname2 = splitline[13]
        let idname2 = parseInt(splitline[12]) - 1
        dicErreur.errorlinks.push([resname1, idname1, resname2, idname2])
    }
    return dicErreur
}

export async function SocketIoPolymerizer(socket: SocketIo.Socket) {
    const WORKDIR = "/data3/rmarin/projet_polyply/job"
    console.log("je suis dans SocketIoPolymerizer")

    // JobManager.start({ 'port': 6001, 'TCPip': "localhost" })
    //     .then(() => {
    //         console.log("JobManager start")


    socket.on("runpolyply", async (dataFromClient: any) => {
        //const data = { polymer: jsonpolymer, density: density, name: name }

        //console.log(dataFromClient)
        //Get forcefield 
        const ff = dataFromClient['polymer']['forcefield']
        //-f martini_v3.0.0_phospholipids_v1.itp
        //const additionalfile = PATHDATA + ff + "/martini_v3.0.0_phospholipids_v1.itp"
        const jsonInStr = JSON.stringify(dataFromClient.polymer)
        const name = dataFromClient['name']
        const density = dataFromClient['density']
        const jobOpt1: JobInputs = {
            "exportVar": {
                "polyplyenv": POLYPLYPATHDATA + "/polyply_1.0/venv/bin/activate",
                "ff": ff,
                "density": density,
                "name": name,
                //"file": additionalfile,
                "action": "itp"
            },
            "inputs": {
                "json": jsonInStr,
                "martiniForceField": POLYPLYPATHDATA + "/martini_v3.0.0.itp",
            }
        }

        const tmp_dir = await TmpDirHelper.get();
        
        console.log("Run polyply gen itp in ", tmp_dir)
        console.log("POLYPLYPATHDATA",POLYPLYPATHDATA )
        try {

            //Run gen itp 
            await ShellManager.run('polyply', jobOpt1, tmp_dir, ShellManager.mode = "jm",);

            let result = fs.readFileSync(tmp_dir + "/jm.stdout").toString();

            console.log(checkResult(result))
            if (checkResult(result)) {
                console.log("Oups pas de fichier gro")
                socket.emit("oups", parseError(result))
            }
            else {

                //Then on fait la requete pour le gro
                console.log("yes on passe au gro")
                const itp = result.split('STOP')[0]
                socket.emit("itp", itp)
                socket.on("continue", () => {
                    //Liste les fichier itp du champs de force 
                    //Tres belle ligne (1h de travail)
                    //let bordelitp = glob.sync(PATHDATA + ff + "/*.itp").map(f => { return "#include " + f + "\r" }).join('')
                   
                    //super stupid variable pour avoir un retour a la ligne
                    const stupidline = '\r'
                    const topfilestr = `#include "${POLYPLYPATHDATA + "/martini_v3.0.0.itp"}"${stupidline}
            #include "polymere.itp"${stupidline}
            [ system ]${stupidline}
            ; name${stupidline}
            mylovelypolymer${stupidline}
            [ molecules ]${stupidline}
            ; name  number${stupidline}
            ${name} 1 ${stupidline}`

                    const jobOpt2: JobInputs = {
                        "exportVar": {
                            "polyplyenv": POLYPLYPATHDATA + "/polyply_1.0/venv/bin/activate",
                            "density": density,
                            "name": name,
                            "top": topfilestr,
                            "action": "gro"
                        },
                        "inputs": {
                            "itp": itp,
                            "martiniForceField": POLYPLYPATHDATA + "/martini_v3.0.0.itp",
                        }
                    }
                    try {
                        ShellManager.run('polyply', jobOpt2, tmp_dir, ShellManager.mode = "jm",);
                        let groJob = fs.readFileSync(tmp_dir + "/jm.stdout").toString();
                        socket.emit("gro", groJob)
                    }
                    catch (e) {
                        console.log(e)
                        // Handle error and throw the right error
                        console.error("ShellManager.run crash");
                    }
                })
            }
        }
        catch (e) {
            console.log(e)
            // Handle error and throw the right error
            console.error("ShellManager.run crash");

        }

    }
    )


}

export default polymer;
