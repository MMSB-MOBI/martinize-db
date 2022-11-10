import logger from '../logger';
import { HISTORY_ROOT_DIR } from '../constants';
import fs, { promises as FsPromise } from 'fs';
import path, { resolve } from 'path';
import { Database } from '../Entities/CouchHelper';
import { Job } from '../Entities/entities'
import { generateSnowflake, getFormattedFile } from "../helpers";
import { isCouchNotFound, notFoundOnFileSystem } from '../Errors';
import { JobFilesNames, JobReadedFiles } from '../types'
import { dateFormatter } from '../helpers'


function getUserJobObject(jobsDoc: Job[]) {
    let obj: { [userId: string]: string[] } = {}
    jobsDoc.forEach(job => {
        if (!obj.hasOwnProperty(job.userId)) obj[job.userId] = [job.id]
        else obj[job.userId].push(job.id)
    })
    return obj
}

export const HistoryOrganizer = new class HistoryOrganizer {
    constructor() {
        logger.info("HistoryOrganizer constructor")
        try {
            fs.mkdirSync(HISTORY_ROOT_DIR)
        } catch { }

    }

    async _saveToFileSystem(job_id: string, files: string[]): Promise<any> {
        const dirPath = HISTORY_ROOT_DIR + "/" + job_id
        fs.mkdirSync(dirPath)
        files.forEach(file => {
            logger.verbose(`copy ${file}`)
            const name = path.basename(file);
            fs.copyFileSync(file, dirPath + "/" + name)
        })
    }

    async _saveToFileSystemFromString(job: any, itp: string, gro: string, top: string, pdb: string): Promise<any> {
        let copyjob = job
        copyjob["files"] = {}
        console.log("_saveToFileSystemFromString", job.jobId)
        const dirPath = HISTORY_ROOT_DIR + "/" + job.jobId
        fs.mkdirSync(dirPath)
        fs.appendFileSync(dirPath + "/" + 'my.itp', itp)

        console.log('Itp Saved!');
        copyjob["files"]['itp_files'] = [['my.itp']]

        fs.appendFileSync(dirPath + "/" + 'my.gro', gro)
        console.log('Gro Saved!');
        copyjob["files"]['gro'] = 'my.gro'

        fs.appendFileSync(dirPath + "/" + 'my.top', top)
        console.log('Top Saved!');
        copyjob["files"]['top_file'] = 'my.top'

        fs.appendFileSync(dirPath + "/" + 'my.pdb', pdb)
        console.log('Pdb Saved!');
        copyjob["files"]['coarse_grained'] = 'my.pdb'


        fs.appendFileSync(dirPath + "/" + 'warning.log', '')
        console.log('warnings Saved!');
        copyjob["files"]['warnings'] = 'warning.log'

        console.log("mon nouveau job", copyjob)
        logger.verbose(`create ${job.jobId}`)
        return copyjob
        // const name = path.basename(file);
        // fs.copyFileSync(file, dirPath + "/" + name)
    }

    async updateJobInFileSystem(jobId: string, itp_files: Express.Multer.File[]) {
        const jobDir = HISTORY_ROOT_DIR + "/" + jobId;
        if (!fs.existsSync(jobDir)) {
            throw new Error("Job directory doesn't exist")
        }
        const newUuid = generateSnowflake()
        const newDir = HISTORY_ROOT_DIR + "/" + newUuid;
        fs.mkdirSync(newDir)
        const files = await FsPromise.readdir(jobDir)
        logger.debug(`copy ${jobId} to ${newUuid}`)
        await Promise.all(files.map(async (currentFile) => {
            await FsPromise.copyFile(jobDir + "/" + currentFile, newDir + "/" + currentFile)
        }))
        logger.debug(`move new itp files to new job directory ${newUuid}`)
        await Promise.all(itp_files.map(async (file) => {
            await FsPromise.rename(file.path, newDir + "/" + file.originalname)
        }))

        return newUuid

    }

    async replaceJobInFileSystem(jobId: string, itp_files: Express.Multer.File[]) {
        const jobDir = HISTORY_ROOT_DIR + "/" + jobId;
        if (!fs.existsSync(jobDir)) {
            throw new Error("Job directory doesn't exist")
        }
        logger.debug(`move new itp files to job directory ${jobId}`)
        return await Promise.all(itp_files.map(async (file) => {
            await FsPromise.rename(file.path, jobDir + "/" + file.originalname)
        }))

    }

    async updateJobForSavedBonds(jobId: string, itp_files_names: string[][]) {
        return await Database.job.updateManuallySavedBonds(jobId, itp_files_names);
    }

    async updateJobAndCreateANewOne(jobId: string, newId: string, newItpFiles: string[][], comment?: string) {
        const updateFnc = (doc: Job) => {

            for (const [idx, mol_files] of newItpFiles.entries()) {

                if (doc.files.itp_files.length <= idx) {
                    doc.files.itp_files.push(mol_files)
                }
                else {
                    const newMolFiles = mol_files.filter(itp => !doc.files.itp_files[idx].includes(itp))
                    doc.files.itp_files[idx] = [...doc.files.itp_files[idx], ...newMolFiles]
                }
            }

            if (doc.manual_bonds_edition) return doc
            doc.manual_bonds_edition = true
            return doc
        }

        const originalJob = await Database.job.get(jobId);
        const newDoc = updateFnc(originalJob)
        delete newDoc._id
        delete newDoc._rev
        newDoc.id = newId
        newDoc.jobId = newId
        newDoc.comment = comment
        newDoc.date = dateFormatter("Y-m-d H:i")

        this.saveToCouch(newDoc)

    }

    async deleteFromFileSystem(jobId: string) {
        logger.debug(`Delete ${jobId} from file system`)
        const dirPath = HISTORY_ROOT_DIR + "/" + jobId
        await FsPromise.rmdir(dirPath, { recursive: true });
    }

    async deleteMultipleFromFileSystem(jobIds: string[]) {
        logger.debug(`Delete ${jobIds} from file system`)
        return await Promise.all(jobIds.map(id => FsPromise.rmdir(HISTORY_ROOT_DIR + "/" + id, { recursive: true })))
    }

    async _deleteFromFileSystemIfExists(jobId: string) {
        const dirPath = HISTORY_ROOT_DIR + "/" + jobId
        try {
            await FsPromise.rmdir(dirPath, { recursive: true })
            logger.debug(`${jobId} deleted from file system`)
        } catch (e) {
            if (notFoundOnFileSystem(e)) logger.debug(`job ${jobId} doesn't exist on file system, no deletion`)
            else throw (e)
        }
    }

    async deleteFromCouch(jobId: string) {
        const job = await Database.job.get(jobId)
        const user = job.userId
        return await Promise.all([Database.job.delete(job), Database.history.deleteJobs(user, [jobId])]);
    }

    async deleteMultipleFromCouch(jobIds: string[]) {
        logger.debug(`Delete multipe ${jobIds} from couch`)
        const jobs = await Database.job.bulkGet(jobIds)
        let notFoundIdx: number[] = [];
        const filteredJobs = jobs.filter((job, idx) => {
            if (job === null) {
                notFoundIdx.push(idx);
                return false
            }
            return true
        })

        const usersRelatedToJobs = getUserJobObject(filteredJobs)
        let promises = jobs.map(job => Database.job.delete(job))
        for (const [user, jobIds] of Object.entries(usersRelatedToJobs)) {
            promises.push(Database.history.deleteJobs(user, jobIds))
        }
        return await Promise.all(promises);

    }

    async _deleteFromCouchIfExists(jobId: string) {
        try {
            const job = await Database.job.get(jobId);
            const user = job.userId
            await Promise.all([Database.job.delete(job), Database.history.deleteJobs(user, [jobId])])
            logger.debug(`${jobId} deleted from couch`)
        } catch (e) {
            if (isCouchNotFound(e)) logger.debug(`job ${jobId} doesn't exist on couch, no deletion`)
            else throw (e)
        }
    }

    async saveToCouch(job: any) {
        const jobDoc = { id: job.jobId, ...job }
        await Database.job.addToJob(jobDoc)
        await Database.history.addToHistory(job.userId, job.jobId)
    }

    async saveToHistory(job: any, files: string[]) {
        return new Promise((res, rej) => {
            this._saveToFileSystem(job.jobId, files).then(() => {
                this.saveToCouch(job).then(() => res(job.jobId)).catch(e => {
                    this.deleteFromFileSystem(job.jobId)
                    rej(e)
                })
            }).catch(e => rej(e))
        })

    }

    async saveToHistoryFromPolyply(job: any, itp: string, gro: string, top: string, pdb: string) {
        console.log("saveToHistoryFromPolyply", job)
        return new Promise((res, rej) => {
            this._saveToFileSystemFromString(job, itp, gro, top, pdb).then((job_new) => {
                console.log("new job", job_new)
                this.saveToCouch(job_new).then(() => res(job_new.jobId)).catch(e => {
                    this.deleteFromFileSystem(job_new.jobId)
                    rej(e)
                })
            }).catch(e => rej(e))
        })
    }


    async getHistory(userId: string) {
        const jobIds = await Database.history.getAllJobs(userId)
        return await Database.job.getJobsDetails(jobIds, userId)
    }

    async getJob(jobId: string): Promise<Job> {
        return new Promise(async (res, rej) => {
            try {
                const job = await Database.job.get(jobId)

                res(job)

            } catch (e) {
                if (isCouchNotFound(e)) {
                    this._deleteFromFileSystemIfExists(jobId)
                    rej("not_found")
                }
                else rej(e)
            }
        })

    }

    async readFiles(jobId: string, files: JobFilesNames): Promise<JobReadedFiles> {
        return new Promise(async (res, rej) => {
            console.log("hellllo", files)
            let readedFiles = {
                all_atom: (Object.keys(files).includes('all_atom')) ? await getFormattedFile(`${HISTORY_ROOT_DIR}/${jobId}/${files.all_atom}`) : { name: jobId, type: "", content: "" },
                top_file: await getFormattedFile(`${HISTORY_ROOT_DIR}/${jobId}/${files.top_file}`),
                coarse_grained: await getFormattedFile(`${HISTORY_ROOT_DIR}/${jobId}/${files.coarse_grained}`),
                itp_files: await Promise.all(files.itp_files.map(async mol_itp => await Promise.all(mol_itp.map(i => getFormattedFile(`${HISTORY_ROOT_DIR}/${jobId}/${i}`))))),
                warnings: await getFormattedFile(`${HISTORY_ROOT_DIR}/${jobId}/${files.warnings}`),
            }

            res(readedFiles)
            // } catch (e) {
            //     if (notFoundOnFileSystem(e)) {
            //         this._deleteFromCouchIfExists(jobId)
            //         rej("not_found")
            //     }
            //     else rej(e)
            // }
        })
    }

    async deleteJobs(jobIds: string[]) {
        return Promise.all([this.deleteMultipleFromCouch(jobIds), this.deleteMultipleFromFileSystem(jobIds)])

    }

    //TO DO : delete history function that delete everything (couch jobs, couch history, file system)



}


export default HistoryOrganizer;