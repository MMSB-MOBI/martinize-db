import { Router } from 'express';
import ItpFile from 'itp-parser';
import SocketIo from 'socket.io';
import { POLYPLYPATHDATA } from "../../constants";
import checkError from './errorParser';
import { JOB_MANAGER_SETTINGS, CONECT_MDP_PATH, CREATE_MAP_PY_SCRIPT_PATH, CREATE_GO_PY_SCRIPT_PATH, DSSP_PATH } from '../../constants';
import jmClient from 'ms-jobmanager'
import JMSurcouche from '../../Builders/JMSurcouche';
import { str_to_stream } from '../../Builders/JMSurcouche';


const router = Router();

let polyplyData: any = {}
const f = async () => {
    console.log("init residue avaible")
    const { stdout, jobFS } = await JMSurcouche.run("get_residue_avaible", { exportVar: {}, inputs: {} })
    return stdout

}; (async () => {
    const res = await f()
    let tempff = ''

    for (let line of res.split("\n")) {

        if (line.startsWith("FORCEFIELD :")) {
            tempff = line.replace("FORCEFIELD :", "")
            polyplyData[tempff] = []
        }
        else polyplyData[tempff].push(line)
    }

})()

router.get('/data', async (req, res) => {
    console.log("Sending forcefields and residues data")
    //Select only martini forcefield
    let MARTINIpolyplyData = Object.fromEntries(Object.entries(polyplyData).filter(([key]) => key.includes('martini')));
    res.send(MARTINIpolyplyData);
});

router.get("/fastaconversion", (req, res) => {
    const fasta = {
        'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
        'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
        'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
        'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'
    }
    res.send(fasta);
})



export async function SocketIoPolymerizer(socket: SocketIo.Socket) {
    socket.on("runpolyply", async (dataFromClient: any) => {
        console.log("Run polyply gen itp")

        //Get forcefield 
        const ff = dataFromClient['polymer']['forcefield']
        const name = dataFromClient['name']
        const boxsize = dataFromClient['box']
        const numberpolymer = dataFromClient['number']

        let additionalfile = ""
        if (dataFromClient['customITP'] !== undefined) {
            for (let itpname of Object.keys(dataFromClient['customITP'])) {
                additionalfile = additionalfile + dataFromClient['customITP'][itpname]
                additionalfile = additionalfile + ";NEWITP\n"
            }
        }

        let ffpath = ''
        if (ff == "martini2") {
            ffpath = POLYPLYPATHDATA + "/" + ff + "/martini_v2.3P.itp"
        }
        else ffpath = POLYPLYPATHDATA + "/" + ff + "/martini_v3.0.0.itp"


        const exportVar = {
            ff: ff,
            name: name,
            action: "itp",
        }

        let inputs = {
            'monfichier.itp': str_to_stream(additionalfile),
            "martiniForceField.itp": ffpath,
            "polymer.json": str_to_stream(JSON.stringify(dataFromClient.polymer)),
        }

        let result: string = ""

        try {
            const { stdout, jobFS } = await JMSurcouche.run('polyply', { exportVar, inputs })
            result = stdout
        }
        catch (e) {
            //Socket emit
            throw new Error(`Error with job manager : ${e}`)
        }

        const itp = result.split("STOP\n")[0]
        const error = result.split("STOP\n")[1]

        const errorParsed = checkError(error)

        if (errorParsed.ok !== true) {
            console.log("######## Error ITP ##########", errorParsed)
            errorParsed['itp'] = itp
            socket.emit("oups", errorParsed)
        }
        else {
            //Then on fait la requete pour le gro
            console.log("Let's GrOOOO")
            socket.emit("itp", itp)
        }

        socket.on("continue", async (itp) => {

            const topfilestr = `
#include "${ffpath}"
#include "polymere.itp"
[ system ]
; name
mylovelypolymer
[ molecules ]
; name  number
${name} ${numberpolymer}
`
            const exportVar = {
                box: boxsize,
                name: name,
                action: "gro"
            }

            const inputs = {
                "polymere.itp": str_to_stream(itp),
                "martiniForceField": ffpath,
                "system.top": str_to_stream(topfilestr)
            }


            let resultatGro = ""
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
                console.log("######## Error GRO ##########", errorParsed)
                socket.emit("oups", errorParsed)
            }
            else {
                socket.emit("gro", gro);
                console.log('Minimize GRO to PDB')

                const exportVar = {
                    "basedir": '',

                    "MDP_FILE": CONECT_MDP_PATH
                }

                const inputs = {
                    "polymere.itp": str_to_stream(itp),

                    "em.mdp": POLYPLYPATHDATA + "/em.mdp",
                    "file.gro": str_to_stream(gro),
                    "file.top": str_to_stream(topfilestr),
                }
                // try {

                const { stdout, jobFS } = await JMSurcouche.run('convert', { exportVar, inputs })
                const fileContent = await jobFS.readToString('output-conect.pdb');
                console.log("Good job Sir! PDB done")
                socket.emit("top", topfilestr);
                socket.emit("pdb", fileContent);

                // }
                // catch (e) {
                //     console.log('ERROR WITH GMX CONVERSION', e)
                //     socket.emit("oups", { ok: false, message: 'ERROR WITH GMX CONVERSION', errorlinks: [] })
                //     throw new Error(`Error with job manager : ${e}`)
                // }

            }
        })
    })
}

export default router;
