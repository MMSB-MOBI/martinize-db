import nano from "nano";
import { User } from "./entities";
import bcrypt from 'bcrypt';
import { Database } from "./CouchHelper";

export default class UserDatabase {
  constructor(protected _db: nano.DocumentScope<User>) {}

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

  save(user: User) {
    return this._db.insert({ _id: user.id, ...user });
  }

  async get(id: string) {
    return this._db.get(id);
  }

  async exists(user: User | string) {
    try {
      await this.get(typeof user === 'string' ? user : user.id);
    } catch (e) {
      return false;
    }
    return true;
  }

  /** Get the user related to given token. undefined if it doesn't exists. */
  async fromToken(jti: string) {
    try {
      const token = await Database.token.get(jti);
      return await this.get(token.user_id);
    } catch (e) {}

    return undefined;
  }

  async fromUsername(username: string) : Promise<User | undefined> {
    const res = await this.findOne({ selector: { name: username } });
    return res[0];
  }
  
  async fromEmail(email: string) : Promise<User | undefined> {
    const res = await this.findOne({ selector: { email } });
    return res[0];
  }

  /**
   * Unhashed password is required in {password}.
   */
  async setOrChangePassword(user: User | string, new_password: string) {
    if (typeof user === 'string') {
      user = await this.get(user);
    }

    user.password = await bcrypt.hash(new_password, 10);

    const update = await this.save(user);
    user._rev = update.rev;
    user._id = update.id;

    return user;
  }

  async verifyPassword(user: User | string, input_password: string) {
    if (typeof user === 'string') {
      user = await this.get(user);
    }

    return bcrypt.compare(input_password, user.password);
  }

  get db() {
    return this._db;
  }
}
