import CliHelper, { CliListener } from "interactive-cli-helper";
import CouchHelper, { Database } from "../Entities/CouchHelper";
import { User } from "../Entities/entities";
import { USERNAME_REGEX, EMAIL_REGEX } from "../constants";
import { generateSnowflake, withRegex } from "../helpers";
import { UserRole } from "../types";
import { CLI } from "./cli";

const USER_CLI = new CliListener(
  CliHelper.formatHelp("user", {
    commands: {
      list: 'List registred users',
      create: 'Create a new user',
      'get <id>/all': 'Get details about user <id> / about all users',
      'grant <id>': 'Make user <id> an administrator',
      'revoke <id>': 'Make user <id> a curator',
      'wipe <id>/all': 'Delete registred user <id> / all users',
      'lookup <username>/<email>': 'Find user(s) with the following username/email',
    },
    onNoMatch: "Command is incorrect. Type \"user\" for help.",
  })
);

USER_CLI.command('list', async () => {
  const users = await Database.user.all();

  if (!users.length) {
    return `This server does not contain any user. Create a user with user create.`;
  }

  return `Available users are \n- ${users.map(m => m._id).join('\n- ')}`;
});

USER_CLI.command('grant', async rest => {
  if (!rest) {
    return "Please enter a user ID. You can search users with lookup.";
  }

  try {
    var user = await Database.user.get(rest);
  } catch {
    return "User not found.";
  }

  user.role = 'admin';
  await Database.user.save(user);

  return `User ${user.name} has been successfully updated.`;
});

USER_CLI.command('revoke', async rest => {
  if (!rest) {
    return "Please enter a user ID. You can search users with lookup.";
  }

  try {
    var user = await Database.user.get(rest);
  } catch {
    return "User not found.";
  }

  user.role = 'curator';
  await Database.user.save(user);

  return `User ${user.name} has been successfully updated.`;
});

USER_CLI.command('lookup', async rest => {
  if (!rest) {
    return "Please enter a query.";
  }

  const users = await Database.user.find({
    selector: {
      $or: [{
        name: withRegex(rest, true, "i"),
      }, {
        email: withRegex(rest, true, "i"),
      }],
    },
  });

  if (!users.length) {
    return `No user matched your search.`;
  }

  return `Matched users are \n- ${users.map(m => `ID ${m._id}: ${m.name} (${m.email})`).join('\n- ')}`;
});

USER_CLI.command('create', async () => {
  console.log("You're about to create a new user.")

  let name = "";
  let email = "";
  let role = "admin";
  let password = "";
  let fullname = ""; 
  let affiliation = ""; 

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

  //Full name
  while (true) {
    fullname = await CLI.question("New user full name: ");

    if (fullname === ".exit") {
      return "User creation exited.";
    }

    break;
  }

  //Affiliation 
  while (true) {
    affiliation = await CLI.question("New user affiliation: ");

    if (affiliation === ".exit") {
      return "User creation exited.";
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
      fullname, 
      affiliation,
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

USER_CLI.command('get', rest => {
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

USER_CLI.command('wipe', async rest => {
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
