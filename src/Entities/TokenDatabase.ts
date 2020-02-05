import { Token } from "./entities";
import AbstractDatabase from "./AbstractDatabase";

export default class TokenDatabase extends AbstractDatabase<Token> {
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
}
