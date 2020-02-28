import CliHelper, { CliListener } from "interactive-cli-helper";
import CouchHelper, { Database } from "../Entities/CouchHelper";
import { User } from "../Entities/entities";
import readline from 'readline';
import { USERNAME_REGEX, EMAIL_REGEX } from "../constants";
import { generateSnowflake } from "../helpers";
import { UserRole } from "../types";
import { CLI } from "./cli";

const USER_CLI = new CliListener(
  CliHelper.formatHelp("user", {
    list: 'List registred users',
    create: 'Create a new user',
    'get <id>/all': 'Get details about user <id> / about all users',
    'wipe <id>/all': 'Delete registred user <id> / all users',
  })
);

USER_CLI.addSubListener('list', async () => {
  const users = await Database.user.all();

  if (!users.length) {
    return `This server does not contain any user.`;
  }

  return `Available users are \n- ${users.map(m => m._id).join('\n- ')}`;
});

USER_CLI.addSubListener('create', async () => {
  console.log("You're about to create a new user.")

  let name = "";
  let email = "";
  let role = "admin";
  let password = "";

  console.log("To exit user creation, type \".exit\"");

  // Name
  while (true) {
    name = await CLI.question("New user name: ");

    if (name === ".exit") {
      return "User creation exited.";
    }

    if (!name.match(USERNAME_REGEX)) {
      console.log("Username is not correct. Please enter a valid user name. To exit user creation: .exit");
      continue;
    }

    const exists = await Database.user.fromUsername(name);
    if (exists) {
      console.log("This user already exists. Please choose another user name.");
      continue;
    }

    break;
  }

  // Email
  while (true) {
    email = await CLI.question("New user email address: ");

    if (email === ".exit") {
      return "User creation exited.";
    }

    if (!email.match(EMAIL_REGEX)) {
      console.log("Email is not correct. Please enter a valid email address. To exit user creation: .exit");
      continue;
    }

    const exists = await Database.user.fromEmail(email);
    if (exists) {
      console.log("This email is already registred. Please choose another email address.");
      continue;
    }

    break;
  }

  // Role
  while (true) {
    role = await CLI.question("New user role. Available roles: \"curator\" or \"admin\": ");

    if (role === ".exit") {
      return "User creation exited.";
    }

    if (role !== "admin" && role !== "curator") {
      console.log("Role is not correct. Available roles: \"curator\" or \"admin\". To exit user creation: .exit");
      continue;
    }

    break;
  }

  // Password
  while (true) {
    password = await CLI.question("New user password: ");

    if (!password.trim()) {
      console.log("Password can't be empty.");
      continue;
    }

    break;
  }

  try {      
    // Create user
    const user: User = {
      id: generateSnowflake(),
      name,
      email,
      role: role as UserRole,
      created_at: new Date().toISOString(),
      password: "",
      approved: true,
    };
  
    const res = await Database.user.save(user);

    if (!res.ok) {
      return "Unable to save user.";
    }

    user._id = res.id;
    user._rev = res.rev;

    await Database.user.setOrChangePassword(user, password);
  
    return "User has been successfully created.";
  } catch (e) {
    return "An error occured during user save (" + (e as Error).stack + ")";
  }
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
