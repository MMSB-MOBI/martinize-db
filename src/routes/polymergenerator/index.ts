import { Router } from 'express';
import SocketIo from 'socket.io';
import { POLYPLYPATHDATA, POLYPLY_VERSION } from "../../constants";
import checkError from './errorParser';
import { JOB_MANAGER_SETTINGS, CONECT_MDP_PATH, CREATE_MAP_PY_SCRIPT_PATH, CREATE_GO_PY_SCRIPT_PATH, DSSP_PATH } from '../../constants';
import JMSurcouche from '../../Builders/JMSurcouche';
import { str_to_stream } from '../../Builders/JMSurcouche';
import HistoryOrganizer from '../../HistoryOrganizer';
import logger from '../../logger';
import { PolyplyJobToSave } from '../molecule/molecule.types';
import { dateFormatter, generateSnowflake } from '../../helpers';
import { createReadStream, ReadStream } from 'fs';
import ItpFile from 'itp-parser-forked';

const router = Router();

const ffDefReadStream = (ff_name:string) : ReadStream => {
    logger.info(`[polymer_generator:ffDefReadStream] ${POLYPLYPATHDATA} and ${ff_name}`);
    const fpath = ff_name == "martini2"
    ? POLYPLYPATHDATA + "/" + ff_name + "/martini_v2.1-dna.itp"
    : POLYPLYPATHDATA + "/" + ff_name + "/martini_v3.0.0.itp";
    logger.info(`[polymer_generator:ffDefReadStream] opening ff itp @ ${fpath}`);
    const s = createReadStream(fpath);
    return s;
}



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

// The above 2 routes should be deprecated
router.get('/data', async (req, res) => {
    console.log("WESHH")
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

const createVirtBoundPolymerMix = (list_graph_component:[string], itp:string) : string =>  {
    //Need to add fake links 

    let previous_res = list_graph_component[0][0];
    let itpparsed = ItpFile.readFromString(itp);
    const atoms = itpparsed.getField('atoms', true)

    const splititp = itp.split("[ bonds ]")
    let itpSTART = splititp[0] + "[ bonds ]\n"

    for (let i of list_graph_component.slice(1)) {
        let next_res = i[0]
        logger.info(`[polymeGenerator:createVirtBoundPolymerMix] need to link  ${previous_res} with ${next_res}`);
        let previousbead = ''
        let nextbead = ''
        for (let i of atoms) {
            if (i.split(' ').filter((e) => { return e !== "" })[2] == previous_res) {
                previousbead = i.split(' ').filter((e) => { return e !== "" })[0]
            }
            if (i.split(' ').filter((e) => { return e !== "" })[2] == next_res) {
                //console.log( "nextbead", i.split(' ').filter((e) => { return e !== "" })[0] )
                nextbead = i.split(' ').filter((e) => { return e !== "" })[0]
            }
        }
        //console.log(i)
        //1  3 1 0.350 4000
        //create a new link 
        let new_link = previousbead + ' ' + nextbead + ' 6 1 1000 ;FAKE LINK\n'
        console.log("new_link",new_link)
        itpSTART = itpSTART + new_link

        previous_res = next_res
    }
    const copy_itp = itpSTART + splititp[1];
    logger.info(`[polymeGenerator:createVirtBoundPolymerMix] rebuildt itp:\n${copy_itp}`);
    return copy_itp;
}


export async function SocketIoPolymerizer(socket: SocketIo.Socket) {
    logger.info(`[route:SocketIoPolymerizer] connected`);
    logger.info(POLYPLYPATHDATA);

    socket.on('version', () => {
        logger.info(POLYPLY_VERSION)
        socket.emit("version_answer", POLYPLY_VERSION)
    })

    socket.on("get_polyply_data", async () => {
    /*
     Returns the dictionnary of available polyply building block as 3 letters code
     Required to initialize client, called once per session
    */
        if (Object.keys(polyplyData).length === 0) await get_truc()
        console.log("Sending forcefields and residues data")
        //Select only martini forcefield
        let MARTINIpolyplyData = Object.fromEntries(Object.entries(polyplyData).filter(([key]) => key.includes('martini')));
        socket.emit("polyply_data", MARTINIpolyplyData);
    }
    )

    socket.on("run_itp_generation", async (dataFromClient: any) => {
        /*
            Returns ITP file of the desired polymer
            the polymer is passed as a graph encoded in JSON.
            Additional connection rules or custom molecules can pe passed
            under the 'customITP' field.
        */
       
        logger.info("Run polyply gen itp", dataFromClient['name'])

        //Get forcefield 
        const ff = dataFromClient['polymer']['forcefield']
        const name = dataFromClient['name']
        /* 
            Write all the eventual additional molecules
            itps in a dedicated "custom" file
        */
        let additionalfile = ""
        if (dataFromClient['customITP'] !== undefined) {
            for (let itpname of Object.keys(dataFromClient['customITP'])) {
                additionalfile = additionalfile + dataFromClient['customITP'][itpname]
                additionalfile = additionalfile + ";NEWITP\n"
            }
        }

        let coordfile = ""
        if (dataFromClient['proteinGRO'] !== "") {
            coordfile = dataFromClient['proteinGRO']
        }


        const exportVar = {
            ff: ff,
            name: name,
            action: "itp",
        }

        let inputs = {
            "monfichier.itp": str_to_stream(additionalfile),
          //  "martiniForceField.itp": ffDefReadStream(ff),
            "polymer.json": str_to_stream(JSON.stringify(dataFromClient.polymer)),
        }

        let result: string = ""

        try {
            const { stdout, jobFS } = await JMSurcouche.run('polyply', { exportVar, inputs })
            result = stdout
        }
        catch (e:any) {
            //Socket emit
            socket.emit("error_itp", e.toString())
            throw new Error(`Error with job manager : ${e}`)
        }


        const itp = result.split("STOP\n")[0]
        const error = result.split("STOP\n")[1]

        const errorParsed = checkError(error)

        if (errorParsed.ok !== true) {
            logger.error(`[route:polymer_generator::run_itp_generation] Error: \"${errorParsed}\"`);
            errorParsed['itp'] = itp
            socket.emit("oups", errorParsed)
        }
        else {
            logger.info('emit itp')
            socket.emit("itp", itp)
        }
    })

    socket.on("run_gro_generation", async (data) => {
        /*
        Second stage, following 'emit itp' reception by client
        The ITPs, polymer number, box dimensions are passed under the data parameter
        */
        logger.info(`[route:polymer_generator::run_gro_generation] running`);

        let coordfile = ""
        if (data['proteinGRO'] !== "") {
            coordfile = data['proteinGRO']
        }


        //Get forcefield 
        const ff = data['polymer']['forcefield']
        const name = data['name']
        const boxsize = data['box']
        const numberpolymer = data['number']
        let itp = data['itp']

        let inputpdb = ""
        if (data['inputpdb'] !== undefined) {
            inputpdb = data['inputpdb']
        }

        let exportVar  = {}
        let inputs     = {}
        let topfilestr = ''

        topfilestr     = `
#include "martini_force_field.itp"
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
        
        if (data["list_graph_component"].length > 1) 
            itp = createVirtBoundPolymerMix(data["list_graph_component"], itp);

        inputs = {
            "coord.gro": str_to_stream(coordfile),
            "polymere.itp": str_to_stream(itp),
            //"martini_force_field.itp": ffDefReadStream(ff),
            "martini_force_field.itp": ff == "martini2"
                ? POLYPLYPATHDATA + "/" + ff + "/martini_v2.1-dna.itp"
                : POLYPLYPATHDATA + "/" + ff + "/martini_v3.0.0.itp",
            "system.top": str_to_stream(topfilestr)
        }

        let resultatGro = ""
        try {
            const { stdout, jobFS } = await JMSurcouche.run('polyply', { exportVar, inputs })
            resultatGro = stdout
        }
        catch (e:any) {
            logger.error(`[route:polymer_generator::run_gro_generation] ${e}`);
            socket.emit("error_gro", e.stderr)
           // throw new Error(`Error with job manager : ${e}`) // Romuald said this should not happen anymore
        }

        const gro = resultatGro.split("STOP\n")[0]
        const error = resultatGro.split("STOP\n")[1]

        /////////CHECK IF GRO IS EMPTY 

        const errorParsed = checkError(error)
        if (errorParsed.ok !== true) {
            logger.error(`[route:polymer_generator::run_gro_generation] ${errorParsed}`);
            socket.emit("oups", errorParsed)
        }
        else {
            logger.info(`[route:polymer_generator::run_gro_generation] success, emittting .gro and .top`);
            socket.emit("gro_top", { gro: gro, top: topfilestr });
        }
    });

    socket.on("run_pdb_generation", async (data) => {
        /*
        Converting GRO into PDB for 3D vizu, running a quick minimization
        */

        

        const ff = data['polymer']['forcefield']; // Nice dts !!
        logger.info(`[route:polymer_generator::run_pdb_generation] starting GRO to PDB: conversion and quick minimize `);
        const exportVar = {
            "basedir": '',
            "MDP_FILE": CONECT_MDP_PATH
        }
        const inputs = {
            "polymere.itp": str_to_stream(data['itp']),
            "em.mdp": POLYPLYPATHDATA + "/em.mdp",
            "martini_force_field.itp": ff == "martini2"
                ? POLYPLYPATHDATA + "/" + ff + "/martini_v2.1-dna.itp"
                : POLYPLYPATHDATA + "/" + ff + "/martini_v3.0.0.itp",
            "file.gro": str_to_stream(data['gro']),
            "file.top": str_to_stream(data['top']),
        }
        try {
            const { stdout, jobFS } = await JMSurcouche.run('convert', { exportVar, inputs })
            const fileContent = await jobFS.readToString('output-conect.pdb');
            logger.info(`[route:polymer_generator::run_pdb_generation] Success, Returning minimized PDB`);
            socket.emit("pdb", fileContent);
        }
        catch (e) {
            logger.info(`[route:polymer_generator::run_pdb_generation] GMX, conversion/minimization Error.`);
            socket.emit("oups", { ok: false, message: 'Error during Gromacs convertion. Please try to increase the box size.', errorlinks: [] })
          //  throw new Error(`Error with job manager : ${e}`) // Romuald said this should not happen anymore
        }

    })

    socket.on("add_to_history", async (d) => {
        logger.info(`[route:polymer_generator::add_to_history] Trying to save ${d["name"]} ${d["polymer"]["forcefield"]}`);

        const gro = d["gro"]
        const pdb = d["pdb"]
        const top = d["top"]
        const itp = d["itp"]
        const name = d["name"]

        let settings = {};
        if (d["polymer"]["forcefield"] === "martini3") 
            settings = { ff: "martini3001", position: "none", cter: "COOH-ter", nter: "NH2-ter", sc_fix: false, cystein_bridge: "auto", builder_mode: "classic", send_mail: false, user_id: d["userId"] }
        else if (d["polymer"]["forcefield"] === "martini2") 
            settings = { ff: "martini22", position: "none", cter: "COOH-ter", nter: "NH2-ter", sc_fix: false, cystein_bridge: "auto", builder_mode: "classic", send_mail: false, user_id: d["userId"] }
        else {
            logger.warn(`[route:polymer_generator::add_to_history] Unregistred forcefield \"${d["polymer"]["forcefield"]}\"`);
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
            logger.info(`[route:polymer_generator::add_to_history] Job ${job.jobId} successfully saved.`);
            socket.emit("add_to_history_answer", job.jobId)
        } catch (e) {
            logger.error(`[route:polymer_generator::add_to_history] Failed to save ${job.jobId}: \"${e}\" `);
            socket.emit("add_to_history_answer", false)
        }
    })

}

export default router;
