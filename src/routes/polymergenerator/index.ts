import { Router } from 'express';
import SocketIo from 'socket.io';
import { POLYPLYPATHDATA } from "../../constants";
import checkError from './errorParser';
import { JOB_MANAGER_SETTINGS, CONECT_MDP_PATH, CREATE_MAP_PY_SCRIPT_PATH, CREATE_GO_PY_SCRIPT_PATH, DSSP_PATH } from '../../constants';
import JMSurcouche from '../../Builders/JMSurcouche';
import { str_to_stream } from '../../Builders/JMSurcouche';
import HistoryOrganizer from '../../HistoryOrganizer';
import logger from '../../logger';
import { PolyplyJobToSave } from '../molecule/molecule.types';
import { dateFormatter, generateSnowflake } from '../../helpers';


const router = Router();

let polyplyData: any = {}
const f = async () => {
    console.log("JMSurcouche.mode", JMSurcouche.mode)
    console.log("Init residue available")
    let ff = ["martini2", "martini3"]
    const { stdout, jobFS } = await JMSurcouche.run("get_residue_avaible", { exportVar: { forcefields: ff.toString() }, inputs: {} })
    return stdout
}

const get_truc = async () => {
    const res = await f()
    let tempff = ''

    for (let line of res.split("\n")) {

        if (line.startsWith("FORCEFIELD :")) {
            tempff = line.replace("FORCEFIELD :", "")
            polyplyData[tempff] = []
        }
        else polyplyData[tempff].push(line)
    }

}

router.get('/data', async (req, res) => {
    if (Object.keys(polyplyData).length === 0) get_truc()
    console.log("Sending forcefields and residues data")
    //Select only martini forcefield
    let MARTINIpolyplyData = Object.fromEntries(Object.entries(polyplyData));
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
    console.log("je suis dans SocketIoPolymerizer ")
    socket.on("get_polyply_data", async () => {
        if (Object.keys(polyplyData).length === 0) await get_truc()
        console.log("Sending forcefields and residues data")
        //Select only martini forcefield
        let MARTINIpolyplyData = Object.fromEntries(Object.entries(polyplyData).filter(([key]) => key.includes('martini')));
        socket.emit("polyply_data", MARTINIpolyplyData);
    }
    )

    socket.on("run_itp_generation", async (dataFromClient: any) => {
        console.log("Run polyply gen itp", dataFromClient['name'])

        //Get forcefield 
        const ff = dataFromClient['polymer']['forcefield']
        const name = dataFromClient['name']

        let additionalfile = ""
        if (dataFromClient['customITP'] !== undefined) {
            for (let itpname of Object.keys(dataFromClient['customITP'])) {
                //console.log("################", dataFromClient['customITP'][itpname])
                additionalfile = additionalfile + dataFromClient['customITP'][itpname]
                additionalfile = additionalfile + ";NEWITP\n"
            }
        }

        let coordfile = ""
        if (dataFromClient['proteinGRO'] !== "") {
            coordfile = dataFromClient['proteinGRO']
        }


        let ffpath = ''
        if (ff == "martini2") {
            ffpath = POLYPLYPATHDATA + "/" + ff + "/martini_v2.1-dna.itp"
        }
        else ffpath = POLYPLYPATHDATA + "/" + ff + "/martini_v3.0.0.itp"


        const exportVar = {
            ff: ff,
            name: name,
            action: "itp",
        }

        let inputs = {
            "monfichier.itp": str_to_stream(additionalfile),
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
            socket.emit("error_itp", e.toString())
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
            console.log('emit itp')
            socket.emit("itp", itp)
        }

    })
    socket.on("run_gro_generation", async (data) => {
        console.log("Let's GrOOOO")

        let coordfile = ""
        if (data['proteinGRO'] !== "") {
            coordfile = data['proteinGRO']
        }


        //Get forcefield 
        const ff = data['polymer']['forcefield']
        const name = data['name']
        const boxsize = data['box']
        const numberpolymer = data['number']
        const itp = data['itp']

        let inputpdb = ""
        if (data['inputpdb'] !== undefined) {
            inputpdb = data['inputpdb']
        }

        let exportVar = {}
        let inputs = {}
        let topfilestr = ''
        let ffpath = ''
        if (ff == "martini2") {
            ffpath = POLYPLYPATHDATA + "/" + ff + "/martini_v2.1-dna.itp"
        }
        else ffpath = POLYPLYPATHDATA + "/" + ff + "/martini_v3.0.0.itp"

        topfilestr = `
#include "${ffpath}"
#include "polymere.itp"
[ system ]
; name
mylovelypolymer
[ molecules ]
; name  number
${name} ${numberpolymer}
`
        exportVar = {
            box: boxsize,
            name: name,
            action: "gro"
        }

        inputs = {
            "coord.gro": str_to_stream(coordfile),
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
            console.log("A l'aide je veux mourir", e)
            socket.emit("error_gro", e.stderr)
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
            console.log("emit go et top")
            socket.emit("gro_top", { gro: gro, top: topfilestr });
        }
    })


    socket.on("run_pdb_generation", async (data) => {
        console.log('Convert and minimize GRO to PDB')
        const exportVar = {
            "basedir": '',
            "MDP_FILE": CONECT_MDP_PATH
        }
        const inputs = {
            "polymere.itp": str_to_stream(data['itp']),
            "em.mdp": POLYPLYPATHDATA + "/em.mdp",
            "file.gro": str_to_stream(data['gro']),
            "file.top": str_to_stream(data['top']),
        }
        try {
            const { stdout, jobFS } = await JMSurcouche.run('convert', { exportVar, inputs })
            const fileContent = await jobFS.readToString('output-conect.pdb');
            console.log("Good job Sir! PDB done")
            socket.emit("pdb", fileContent);
        }
        catch (e) {
            console.log('ERROR WITH GMX CONVERSION', e)
            socket.emit("oups", { ok: false, message: 'ERROR WITH GMX CONVERSION', errorlinks: [] })
            throw new Error(`Error with job manager : ${e}`)
        }

    })

    socket.on("add_to_history", async (d) => {
        console.log("Add to history", d)

        const gro = d["gro"]
        const pdb = d["pdb"]
        const top = d["top"]
        const itp = d["itp"]
        const name = d["name"]

        console.log(d["polymer"]["forcefield"])

        let settings
        if (d["polymer"]["forcefield"] === "martini3") settings = { ff: "martini3001", position: "none", cter: "COOH-ter", nter: "NH2-ter", sc_fix: false, cystein_bridge: "auto", builder_mode: "classic", send_mail: false, user_id: d["userId"] }
        else if (d["polymer"]["forcefield"] === "martini2") settings = { ff: "martini22", position: "none", cter: "COOH-ter", nter: "NH2-ter", sc_fix: false, cystein_bridge: "auto", builder_mode: "classic", send_mail: false, user_id: d["userId"] }
        else {
            console.log("Hello sir, Il y a une erreur ici avec le forcefield !", d["polymer"]["forcefield"])
            settings = { ff: "martini3001", position: "none", cter: "COOH-ter", nter: "NH2-ter", sc_fix: false, cystein_bridge: "auto", builder_mode: "classic", send_mail: false, user_id: d["userId"] }
        }

        const job: PolyplyJobToSave = {
            jobId: generateSnowflake(),
            userId: d["userId"],
            type: "polyply",
            date: dateFormatter("Y-m-d H:i"),
            //@ts-ignore
            settings,
            name: name,

        }

        try {
            await HistoryOrganizer.saveToHistoryFromPolyply(job, itp, gro, top, pdb)
            console.log("job.jobId", job.jobId)
            socket.emit("add_to_history_answer", job.jobId)
        } catch (e) {
            logger.warn("error save to history", e)
            socket.emit("add_to_history_answer", false)
        }
    })

}

export default router;
