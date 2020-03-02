import nano from 'nano';
import TokenDatabase from './TokenDatabase';
import MoleculeDatabase from './MoleculeDatabase';
import UserDatabase from './UserDatabase';
import StashedMoleculeDatabase from './StashedMoleculeDatabase';
import { COUCH_URL } from '../constants';
import logger from '../logger';

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
}

export default class CouchHelper {
  // @ts-ignore
  public link: nano.ServerScope;
  // @ts-ignore
  protected dbs: Databases;

  static readonly MOLECULE_COLLECTION = "molecule";
  static readonly STASHED_MOLECULE_COLLECTION = "stashed";
  static readonly USER_COLLECTION = "user";
  static readonly TOKEN_COLLECTION = "token";
  static readonly DBS = [
    "molecule",
    "stashed",
    "user",
    "token"
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
    };
  }

  ping() {
    return this.link.db.list();
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

export const Database = new CouchHelper(COUCH_URL);
