


// WARNING - general - Missing link between residue 1 SER and residue 4 SER

interface ErrorToClient {
    ok: boolean,
    disjoint: boolean,
    errorlinks: any[],
    message: string[],
    itp? : string,
}


export default checkError;

function checkError(output: string) {
    //Init dico error
    let dicErreur: ErrorToClient = { ok: true, disjoint: false, message: [], errorlinks: [] }
    let pythonerror: boolean = false
    let oserror: boolean = false
    //Parse every line 
    for (let l of output.split('\n')) {
        if (l == '') continue

        if (pythonerror === true) {
            dicErreur.message.push(l)
        }

        if (oserror === true) {
            dicErreur.message.push(l)
        }

        if ((l.includes('Traceback (most recent call last):'))) {
            console.log("error disconnected parts", l)
            dicErreur.ok = false
            pythonerror = true
        }

        if ((l.includes('disconnected parts. ')) || (l.includes('disjoint parts'))) {
            console.log("error disconnected parts", l)
            dicErreur.disjoint = true
            dicErreur.ok = false
        }

        if ((l.includes('unrecognized arguments'))) {
            console.log("unrecognized arguments", l)
            dicErreur.message.push(l)
            dicErreur.ok = false
        }

        if (l.includes('disjoint parts')) {
            console.log("error disjoint parts", l)
            dicErreur.disjoint = true
            dicErreur.ok = false
        }
        // WARNING - general - Missing link between residue 1 SER and residue 4 SER

        if (l.includes('Missing link')) {
            console.log("error", l)
            let splitline = l.split(' ')
            let resname1 = splitline[9]
            let idname1 = parseInt(splitline[8]) - 1
            let resname2 = splitline[13]
            let idname2 = parseInt(splitline[12]) - 1
            dicErreur.errorlinks.push([resname1, idname1, resname2, idname2])
            console.log(splitline)
        }

        if (l.includes('OSError:')) {
            dicErreur.ok = false
            dicErreur.message.push(l)
            oserror = true
        }

    }
    return dicErreur
}
