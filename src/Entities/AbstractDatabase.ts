import nano from "nano";
import logger from "../logger";

export default abstract class AbstractDatabase<T extends { id: string, _id?: string, _rev?: string }> {
  constructor(protected _db: nano.DocumentScope<T>) {}

  /** Get all documents in the database */
  async all() {
    const res = await this._db.list({ include_docs: true });
    return res.rows.map(d => d.doc).filter(d => d) as (T & nano.Document)[];
  }

  /** Bulk create */
  bulkCreate(docs: T[]) {
    const original_length = docs.length;
    if ((docs = docs.filter(d => !d._rev)).length !== original_length) {
      logger.warn("bulkCreate: One element has already been inserted, this should not happen. It has been filtered.");
    }

    return this.bulk({ docs });
  }

  /** Bulk update */
  bulkUpdate(docs: T[]) {
    const original_length = docs.length;
    if ((docs = docs.filter(d => d._rev)).length !== original_length) {
      logger.warn("bulkUpdate: One element has no revision information, this should not happen. It has been filtered.");
    }

    return this.bulk({ docs });
  }

  /** Bulk delete */
  bulkDelete(documents: T[]) {
    const original_length = documents.length;
    if ((documents = documents.filter(d => d._rev)).length !== original_length) {
      logger.warn("bulkDelete: One element has no revision information, this should not happen. It has been filtered.");
    }

    const docs: nano.BulkModifyDocsWrapper = { docs: [] };
    for (const doc of documents) {
      docs.docs.push({ _deleted: true, _id: doc._id, _rev: doc._rev });
    }

    return this.bulk(docs);
  }

  /** Bulk create, update or delete */
  bulk(docs: nano.BulkModifyDocsWrapper) {
    return this._db.bulk(docs);
  }

  /**
   * Warning: By default, {query.limit} is 25 !
   */
  async find(query: nano.MangoQuery) {
    const res = await this._db.find(query);
    return res.docs;
  }

  async findOne(query: nano.MangoQuery) {
    const res = await this._db.find({ ...query, limit: 1 });
    return res.docs;
  }

  save(element: T) {
    return this._db.insert({ _id: element.id, ...element });
  }

  get(id: string) {
    return this._db.get(id);
  }

  delete(element: T) {
    return this._db.destroy(element._id!, element._rev!);
  }

  async count() {
    const info = await this._db.info();
    return info.doc_count;
  }

  async exists(element: T | string) {
    try {
      await this.get(typeof element === 'string' ? element : element.id);
    } catch (e) {
      return false;
    }
    return true;
  }

  async isCreated() : Promise<boolean> {
    try {
      await this._db.info();
    } catch (e) {
      return false;
    }
    return true;
  }

  get db() {
    return this._db;
  }
}
