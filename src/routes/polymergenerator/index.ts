import { Router } from 'express';
import glob from 'glob';
import ItpFile from 'itp-parser';
import SocketIo from 'socket.io';
import ShellManager, { JobInputs } from '../../Builders/ShellManager';
import TmpDirHelper from '../../TmpDirHelper';
import * as fs from 'fs';
import { MINIMIZEPDB, POLYPLYPATHDATA, POLYPLY_PATH, POLYPLY_VENV } from "../../constants";
import checkError from './errorParser';
import { JOB_MANAGER_SETTINGS, CONECT_MDP_PATH, CREATE_MAP_PY_SCRIPT_PATH, CREATE_GO_PY_SCRIPT_PATH, DSSP_PATH } from '../../constants';
import jmClient from 'ms-jobmanager'
import JMSurcouche from '../../Builders/JMSurcouche';
import { str_to_stream } from '../../Builders/JMSurcouche';


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


    //remove duplkiicate 
    for (let forcefield of Object.keys(avaibleData)) {
        avaibleData[forcefield] = [...new Set(avaibleData[forcefield])];
    }

    res.send(avaibleData);
    console.log("Sending forcefields and residues data")
});

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
    socket.on("runpolyply", async (dataFromClient: any) => {

        const tmp_dir = await TmpDirHelper.get();
        console.log("Run polyply gen itp in ", tmp_dir)

        //Get forcefield 
        const ff = dataFromClient['polymer']['forcefield']

        const name = dataFromClient['name']

        const boxsize = dataFromClient['box']

        let additionalfile = ""
        if (dataFromClient['customITP'] !== undefined) {
            additionalfile = dataFromClient['customITP'].join(";NEWITP\n")
        }

        const exportVar = {
            ff: ff,
            name: name,
            action: "itp",
        }

        let inputs = {
            'monfichier.itp': str_to_stream(additionalfile),
            "martiniForceField.itp": POLYPLYPATHDATA + "/martini_v3.0.0.itp",
            "polymer.json": str_to_stream(JSON.stringify(dataFromClient.polymer)),
        }

        let result: string = ""

        try {
            const { stdout, jobFS } = await JMSurcouche.run('polyply', {exportVar, inputs})
            result = stdout
        }
        catch (e) {
            //Socket emit
            throw new Error(`Error with job manager : ${e}`)
        }


        const itp = result.split("STOP\n")[0]
        const error = result.split("STOP\n")[1]

        const errorParsed = checkError(error)
         
        if (errorParsed.ok == true) {
            //Then on fait la requete pour le gro
            console.log("yes on passe au gro")
            socket.emit("itp", itp)
            socket.on("continue", async () => {


                const topfilestr = `
#include "${POLYPLYPATHDATA + "/martini_v3.0.0.itp"}"
#include "${POLYPLYPATHDATA + "/martini_v3.0.0_solvents_v1.itp"}"
#include "polymere.itp"
[ system ]
; name
mylovelypolymer
[ molecules ]
; name  number
${name} 1
`
                const exportVar = {
                    box: boxsize,
                    name: name,
                    action: "gro"
                }

                const inputs = {
                    "polymere.itp": str_to_stream(itp),
                    "martiniForceField": POLYPLYPATHDATA + "/martini_v3.0.0.itp",
                    "system.top": str_to_stream(topfilestr)
                }

                try {
                    let resultatGro = " "
                    try {
                        const { stdout, jobFS } = await JMSurcouche.run('polyply', { exportVar, inputs })
                         
                        resultatGro = stdout
                    }
                    catch (e) {
                        throw new Error(`Error with job manager : ${e}`)
                    }

                    const gro = resultatGro.split("STOP\n")[0]
                    const error = resultatGro.split("STOP\n")[1]

                    /////////CHECK IF GRO IS EMPTY 

                    const errorParsed = checkError(error)
                    if (errorParsed.ok !== true) {
                        console.log("######## Error ITP ##########", errorParsed)
                        socket.emit("oups", errorParsed)
                    }
                    else {
                        console.log('socket.emit("gro", gro);')
                        console.log("pouet", gro)
                        socket.emit("gro", gro);

                        console.log('Lancement du gmx')
                        
                        const exportVar = {
                            "basedir": '',
                            "DEL_WATER_BOOL": "NO",
                            "MDP_FILE": CONECT_MDP_PATH
                        }

                        const inputs = {
                            "polymere.itp": str_to_stream(itp),
                            "water.gro": POLYPLYPATHDATA + "/water.gro",
                            "em.mdp": POLYPLYPATHDATA + "/em.mdp",
                            "file.gro": str_to_stream(gro),
                            "file.top": str_to_stream(topfilestr),
                        }
                        try {
                             
                            const { stdout, jobFS } = await JMSurcouche.run('convert', {exportVar, inputs })

                            const fileContent = await jobFS.readToString('output-conect.pdb');
                            console.log("Bravo monsieur! PDB done")
                            socket.emit("pdb", fileContent);

                        }
                        catch (e) {
                            console.log('ERROR WITH GMX CONVERSION')
                            throw new Error(`Error with job manager : ${e}`)
                        }

                    }

                }
                catch (e) {
                    console.log('ERROR WITH GRO')
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
