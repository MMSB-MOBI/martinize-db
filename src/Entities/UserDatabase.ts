import { User } from "./entities";
import bcrypt from 'bcrypt';
import { Database } from "./CouchHelper";
import AbstractDatabase from "./AbstractDatabase";
import { withRegex } from "../helpers";

export default class UserDatabase extends AbstractDatabase<User> {
  /** Get the user related to given token. undefined if it doesn't exists. */
  async fromToken(jti: string) {
    try {
      const token = await Database.token.get(jti);
      return await this.get(token.user_id);
    } catch (e) {}

    return undefined;
  }

  async fromUsername(username: string) : Promise<User | undefined> {
    return this.findOne({ selector: { name: withRegex(username, false, "i", true) } });
  }
  
  async fromEmail(email: string) : Promise<User | undefined> {
    return this.findOne({ selector: { email: withRegex(email, false, "i", true) } });
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
}
