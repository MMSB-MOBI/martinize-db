import { CliListener } from "interactive-cli-helper";
import CouchHelper, { Database } from "../Entities/CouchHelper";
import { User } from "../Entities/entities";

const USER_CLI = new CliListener("Available commands are list, get, wipe");

USER_CLI.addSubListener('list', async () => {
  const users = await Database.user.all();

  if (!users.length) {
    return `This server does not contain any user.`;
  }

  return `Available users are \n- ${users.map(m => m._id).join('\n- ')}`;
});

USER_CLI.addSubListener('get', rest => {
  rest = rest.trim();

  if (!rest) {
    return `Please specify a user id.`;
  }

  if (rest === "all") {
    return Database.user.all();
  }

  try {
    BigInt(rest);
  } catch (e) {
    return `ID ${rest} is not valid, please enter a valid number.`;
  }

  return Database.user.get(rest);
});

USER_CLI.addSubListener('wipe', async rest => {
  rest = rest.trim();
  
  if (!rest) {
    return `Please specify a user id or "all".`;
  }

  if (rest === "all") {
    await Database.delete(CouchHelper.USER_COLLECTION);
    await Database.create(CouchHelper.USER_COLLECTION);
    return `User database is wiped`;
  }

  try {
    BigInt(rest);
  } catch (e) {
    return `ID ${rest} is not valid, please enter a valid number.`;
  }

  let user: User;
  try {
    user = await Database.user.get(rest);
  } catch (e) {
    return `Unable to get user (${rest})`;
  }

  if (user) {
    return Database.user.delete(user);
  }
  return `Unable to find user.`
});

export default USER_CLI;
