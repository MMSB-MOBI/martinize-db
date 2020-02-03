import nano from "nano";
import { Token } from "./entities";

export default class TokenDatabase {
  constructor(protected _db: nano.DocumentScope<Token>) {}

  async find(query: nano.MangoQuery) {
    const res = await this._db.find(query);
    return res.docs;
  }

  save(token: Token) {
    return this._db.insert({ _id: token.id, ...token });
  }

  async get(id: string) {
    return this._db.get(id);
  }

  tokensOf(user_id: string) {
    return this.find({ selector: { user_id } });
  }

  async userFromToken(token: Token | string) {
    try {
      const t = await this.get(typeof token === 'string' ? token : token.id);
      return t.user_id;
    } catch (e) { }
    return undefined;
  }

  async exists(token: Token | string) {
    try {
      await this.get(typeof token === 'string' ? token : token.id);
    } catch (e) {
      return false;
    }
    return true;
  }

  get db() {
    return this._db;
  }
}
