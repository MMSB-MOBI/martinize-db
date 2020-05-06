# Deployment

Here is a guide to configure and deploy this server on a Linux system, Docker or not.

## Requirements

The following software must be installed or available:

- **Node.js** 12+ (this server has been developed for Node 13/14)
- **TypeScript** 3.8+ installed
- **CouchDB** server, can be on another machine
- **GROMACS** 2016+, see part `Configuring GROMACS engine`
- **Python 2** and **Python 3**, with the possibility to create `virtualenv` environnements 
- **OpenSSL**, usually included in all environnements
- **The createGoVirt.py** script

## Getting started

Clone the git repository and run `npm install`

```bash
git clone https://github.com/alkihis/martinize-db.git
cd martinize-db
npm i
```

## Building the server

Git-ed files does not include the compiled version of the server. You must compile TypeScript before running the server.

```bash
# With global TS (must be TS 3.8 or later)
tsc
# With project typescript
npm run tsc
```

## Include front-end website to this server

Front-end compiled website must be situed in `static/www` directory.

```bash
# Delete the possible existing compiled website in martinize-db server
rm -rf static/www/*

# Clone the front-end repo
cd
git clone https://github.com/alkihis/martinize-db-client.git
cd martinize-db-client

# Install dependencies
npm i

# Build the website
npm run build

# Compiled site is in build
cp -R build/* {path/to/martinize-db}/static/www

# Go back to martinize-db server
cd ~/martinize-db
```

## Setup server

The server have a bunch of configuration files and dependencies you need to setup.


### Directories

Ensure the following directories are created and write-access-friendly on server root:

- `uploads`
- `molecules`

Create a directory `.keys` for the next part.

The following directories must exist and be filled:

- `force_fields`: with Martini ITP files
- `lipids/2_2` and `lipids/3_0`: With lipid ITP files for given force field version


### JSON Web Tokens RSA Keys

You must setup JWT keys in order to use the server. Those keys helps server to sign and generate unique token nobody can decode.

**You must have OpenSSL installed to generate them**.

```bash
# If you haven't done it, create .keys server
mkdir .keys
cd .keys

# When asked, do not enter passphrase, press enter
ssh-keygen -t rsa -b 4096 -m PEM -f key_new.pem
openssl rsa -in key_new.pem -pubout -outform PEM -out key_new
```


### Constants

This guide will help you to properly setup constants.

In file `src/constants.ts`:

- `URLS`:
  - `SERVER`: Specify here server **public** URL. This will be used in sent emails.
  - `COUCH`: Specify here CouchDB URL.
- `PYTHON_3_PATH`: Path to Python 3.6+ binary.
- `CREATE_GO_PATH`: Path to `createGoVirt.py` script.
- `DSSP_PATH`: (optional) Path to anaconda's DSSP.
- `MAILER_ENFORCE_RECIPIENT`: Set this to `false`.
- `DEFAULT_MAILER_ADDRESS`: Put here the e-mail address that will send e-mails to clients.

Usually, you don't need to modify other constants.

After configuring those constants, you need to rebuild the server.

```bash
npm run tsc
```


### Python requirements

The following server parts are using Python script/utilities:
- `ccmap`
- `insane`
- `martinize2`

#### ccmap

`CCMap` is started from script `utils/get_map.sh`.

Inspect this script to include the virtualenv `source` in order to use a Python binary that have `ccmap` package installed.

If you want to change this, call to this script is done from `src/Builders/Martinizer.ts` file, method `getCcMap` of singleton `Martinizer`.


#### insane

[INSANE](https://github.com/Tsjerk/Insane) is started from script `utils/insane.sh`.

If you want to include a `source` to a virtualenv that contain INSANE, do it at the start of the script.
**INSANE can be installed via PiP, with a `pip2 install --user insane`**.

Call to this script is done from `src/Builders/MembraneBuilder.ts` file, method `run`.


#### martinize2

[Martinize 2](https://github.com/marrink-lab/vermouth-martinize) is started through `utils/martinize.sh`.

As for INSANE, inspect this script to include the virtualenv `source` in order to use a Python binary that have vermouth-martinize2 package installed.


### Configuring GROMACS engine

GROMACS is used to create PDB files from a set of PDB+ITP files. This operation is **very fast** (up to 5 seconds).

Currently, script `utils/create_conect_pdb.sh` is in charge to handle the pipeline.
This script is called in `createPdbWithConect` method of `Martinizer` singleton (`src/Builders/Martinizer.ts`).

At the start of the script, you'll find a `GROMACS_LOADER` sh variable that you can configure to start GROMACS module.


## Start the server

When all of this is ready, you can start the server.

### Start developement server

```bash
npm run start-dev
```

If everything goes well, it should be configured properly.

### Server parameters

You can set some CLI args to give parameters to server.

```
Options:
  -c, --couchdb-url <url>      Couch DB URL (default: process.env.COUCHDB_HOST || constants.URLS.COUCH)
  --server-url <url>           Server URL (default: process.env.SERVER_URL || constants.URLS.SERVER)
  -p, --port <port>            Emit port (default: 4123)
  --wipe-init                  Destroy the CouchDB bases and re-create them.
  --init-db                    Init the required CouchDB databases.
  --quit-after-init            If --wipe-init or --init-db are set, this will quit server after init.
  -l, --log-level <logLevel>   Log level [debug|silly|verbose|info|warn|error] (default: "info")
  --file-log-level <logLevel>  Log level (written to file) [debug|silly|verbose|info|warn|error] (default: "info")
  --log-file <logFile>         File log level
```

Example:
```bash
# Set env variables, SERVER_URL and COUCHDB_HOST
SERVER_URL=https://martinize-db.ibcp.fr
COUCHDB_HOST=http://localhost:5984

# Run the server on part 4123, with CLI log level verbose, 
# also log in file server.log with minimum level of info
npm run start -- -p 4123 -l verbose --file-log-level info --log-file server.log
```

### Command line interface

Once started, the server have a CLI to interact with.
You can inspect quickly which molecules are inserted, create an user, init a database,
send test mails and more. Type `help` then enter to have informations.


### Start production server / First start

Just simply run:
```bash
npm run start
```

At server first start, you must set up the environnement in order to make user system and database work.

The console will remind you to do it:
```
Welcome to Martinize server CLI. For help, type "help".

WARN: The database seems to be un-initialized. Please create all databases by entering "database create all".
WARN: Once database is created, you can create an administrator account with "user create".
```

If you need to, create the databases inside the CLI by typing: `database create all`.

Then, you must create the first administrator user in order to accept accounts.

```
> user create
You're about to create a new user.
To exit user creation, type ".exit"
New user name: <enter admin username>
New user email address: <enter desired e-mail address>
New user role. Available roles: "curator" or "admin": admin
New user password: <enter here the desired password>
cli: User has been successfully created.
```

To be sure the user has been created, use the CLI:
```
> user get all
cli: [
  {
    _id: '<here is the user ID>',
    _rev: '<here is the couch rev ID>',
    id: '<here is the user ID>',
    name: 'username',
    email: 'user@its.email',
    role: 'admin',
    created_at: '2020-03-09T09:44:02.851Z',
    password: '$2b$10$lJLO7NFdYFVXVEX1VbtcCOH9WnFK35J/1efSFJQSoTM.Dnot7y0pm',
    approved: true
  }
]
```

You can now login on the [website](http://localhost:4123/) !




