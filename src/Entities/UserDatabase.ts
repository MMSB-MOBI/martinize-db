import nano from "nano";
import { User } from "./entities";
import bcrypt from 'bcrypt';

export default class UserDatabase {
  constructor(protected _db: nano.DocumentScope<User>) {}

  async find(query: nano.MangoQuery) {
    const res = await this._db.find(query);
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
