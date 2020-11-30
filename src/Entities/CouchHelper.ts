import nano from 'nano';
import TokenDatabase from './TokenDatabase';
import MoleculeDatabase from './MoleculeDatabase';
import UserDatabase from './UserDatabase';
import StashedMoleculeDatabase from './StashedMoleculeDatabase';
import { URLS, DB_PREFIX } from '../constants';
import logger from '../logger';
import RadiusDatabase from './RadiusDatabase';
import LipidDatabase from './LipidDatabase';

export class CouchDatabase<T> {
  constructor(protected collection: nano.DocumentScope<T>) {}

  get db() {
    return this.collection;
  }
}

interface Databases {
  molecule: MoleculeDatabase;
  stashed: StashedMoleculeDatabase;
  user: UserDatabase;
  token: TokenDatabase;
  radius: RadiusDatabase;
  lipid: LipidDatabase;
}

export default class CouchHelper {
  public link!: nano.ServerScope;
  protected dbs!: Databases;

  static readonly MOLECULE_COLLECTION = `${DB_PREFIX}molecule`;
  static readonly STASHED_MOLECULE_COLLECTION = `${DB_PREFIX}stashed`;
  static readonly USER_COLLECTION = `${DB_PREFIX}user`;
  static readonly TOKEN_COLLECTION = `${DB_PREFIX}token`;
  static readonly RADIUS_COLLECTION = `${DB_PREFIX}vanderwaalsradii`;
  static readonly LIPID_COLLECTION = `${DB_PREFIX}lipid`;
  static readonly DBS = [
    `${DB_PREFIX}molecule`,
    `${DB_PREFIX}stashed`,
    `${DB_PREFIX}user`,
    `${DB_PREFIX}token`,
    `${DB_PREFIX}vanderwaalsradii`,
    `${DB_PREFIX}lipid`,
  ];

  constructor(url: string) {
    this.refresh(url);
  }

  /**
   * Link the given url to collections.
   */
  refresh(url: string) {
    this.link = nano({ url, requestDefaults: { proxy: null } });

    this.dbs = {
      molecule: new MoleculeDatabase(this.link.use(CouchHelper.MOLECULE_COLLECTION)),
      stashed: new StashedMoleculeDatabase(this.link.use(CouchHelper.STASHED_MOLECULE_COLLECTION)),
      user: new UserDatabase(this.link.use(CouchHelper.USER_COLLECTION)),
      token: new TokenDatabase(this.link.use(CouchHelper.TOKEN_COLLECTION)),
      radius: new RadiusDatabase(this.link.use(CouchHelper.RADIUS_COLLECTION)),
      lipid: new LipidDatabase(this.link.use(CouchHelper.LIPID_COLLECTION)),
    };
  }

  async ping() {
    try {
      await this.link.db.use(CouchHelper.USER_COLLECTION).list();

    } catch (e) {
      if (!e.message.startsWith('Database does not exist')) {
        throw e;
      }
    }
  }

  get molecule() {
    return this.dbs.molecule;
  }

  get token() {
    return this.dbs.token;
  }

  get user() {
    return this.dbs.user;
  }

  get stashed() {
    return this.dbs.stashed;
  }

  get radius() {
    return this.dbs.radius;
  }

  get lipid() {
    return this.dbs.lipid;
  }

  /** Create a database */
  create(name: string) {
    return this.link.db.create(name);
  }

  async createAll() {
    for (const db of CouchHelper.DBS) {
      logger.silly(JSON.stringify(await this.create(db)));
    }
  }

  /** Delete a database */
  delete(name: string) {
    return this.link.db.destroy(name).catch(e => e);
  }

  async deleteAll() {
    for (const db of CouchHelper.DBS) {
      await this.delete(db);
    }
  }

  /** Wipe all databases and recreate them all */
  async wipeAndCreate() {
    await this.deleteAll();
    await this.createAll();
  }
}

export const Database = new CouchHelper(URLS.COUCH);
