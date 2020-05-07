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

## Getting started

### Node.js

Node 14 can be installed on RHEL/CentOS/Fedora with the following command:
```bash
curl -sL https://rpm.nodesource.com/setup_14.x | sudo bash -
sudo dnf install nodejs gcc g++ make cmake  # Some pkg needs compilation tools
sudo npm install -g typescript
```

### Repository

Clone the git repository and run `npm install`

```bash
# sudo dnf install git
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
# Create the static/www folder in martinize-db
mkdir static
mkdir static/www

# Clone the front-end repo
cd ..
git clone https://github.com/alkihis/martinize-db-client.git
cd martinize-db-client
# Install dependencies
npm i

cd ..

# you will also need a ngl copy for type declarations (it's annoying but ngl bundle is made in a obscure manner..)
git clone https://github.com/alkihis/ngl.git
cd ngl
npm i
npm run dts
cp -R -f declarations ../martinize-db-client/node_modules/@mmsb/ngl/

cd ../martinize-db-client

# Build the website
npm run build

# Compiled site is in build
cp -R build/* ../martinize-db/static/www

# Go back to martinize-db server
cd ../martinize-db
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
cd ..
```


### Constants

This guide will help you to properly setup constants.

In file `src/constants.ts`:

- `URLS`:
  - `SERVER`: Specify here server **public** URL. This will be used in sent emails.
  - `COUCH`: Specify here CouchDB URL.
- `DSSP_PATH`: (optional) Path to anaconda's DSSP.
- `MAILER_ENFORCE_RECIPIENT`: Set this to `false`.
- `DEFAULT_MAILER_ADDRESS`: Put here the e-mail address that will send e-mails to clients.

Usually, you don't need to modify other constants.

After configuring those constants, you need to rebuild the server.

```bash
npm run tsc
```


### Python requirements

#### Installing venvs on CentOS

```bash
sudo dnf install git python2 python3-devel*  # Needed to install ccmap from source
sudo pip3 install virtualenv

# venv for insane (Python 2)
virtualenv --python=$(which python2) insanevenv
source insanevenv/bin/activate
pip install insane
deactivate

# venv for martinize2/ccmap/create_goVirt (Python 3)
python3 -m venv martinize2venv 
source martinize2venv/bin/activate
pip install --upgrade pip
pip install ccmap pyproteins pyproteinsExt
pip install git+https://github.com/marrink-lab/vermouth-martinize.git#vermouth

# Copy the martini304 force field files into site packages (server should be cloned!)
target=martinize2venv/lib/python3.6/site-packages/vermouth/data
origin=martinize-db/vermouth-data

cp -R $origin/force_fields/* $target/force_fields
cp -R $origin/mappings/* $target/mappings

deactivate
```

The following server parts are using Python script/utilities:
- `ccmap`
- `insane`
- `martinize2`
- `govirt`

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
If you want to change this, call to this script is done from `src/Builders/Martinizer.ts` file, method `run` of singleton `Martinizer`.


### govirt

The Go virt script (found in `utils/create_goVirt.py`) is supersetted by `utils/create_go_virt.sh` to allow usage of virtualenv.

As for INSANE and Martinize2, nspect this script to include the virtualenv `source` in order to use the same virtualenv as martinize2 (numpy is needed).

If you want to change this, call to this script is done from `src/Builders/Martinizer.ts` file, method `run` of singleton `Martinizer`.


### Configuring GROMACS engine

#### Installation guide for a Docker/Virtual machine

```bash
# If needed, install c/c++/make/cmake
bnf install cmake make gcc g++ tar wget

# Download from gromacs server the gromacs source
wget http://ftp.gromacs.org/pub/gromacs/gromacs-2020.2.tar.gz

# Follow the gromacs compile guide
tar xfz gromacs-2020.2.tar.gz
cd gromacs-2020.2
mkdir build
cd build
cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON
make
make check
sudo make install

# To each time you want to use gmx commands: You can create an alias where you want
# This line or a remplacement should be present in `utils/create_conect_pdb.sh` in GROMACS_LOADER variable. 
source /usr/local/bin/GMXRC # or source /usr/local/gromacs/bin/GMXRC
```

#### For every machine

GROMACS is used to create PDB files from a set of PDB+ITP files. This operation is **very fast** (up to 5 seconds).

Currently, script `utils/create_conect_pdb.sh` is in charge to handle the pipeline.
This script is called in `createPdbWithConect` method of `Martinizer` singleton (`src/Builders/Martinizer.ts`).

At the start of the script, you'll find a `GROMACS_LOADER` sh variable that you can configure to start GROMACS module.


## Start the server

When all of this is ready, you can start the server.

### Test if everything is okay

First, start the server in development mode.
```bash
npm run start-dev
```

In the server CLI, type the following command.
```
test dry-run
```

If everything is okay, the following message must appear at the end.
```
cli: Run end successfully, you can inspect out files in {a_temporary_directory}.
```

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
SERVER_URL=https://martinize-db.ibcp.fr
COUCHDB_HOST=http://localhost:5984

# Run the server on part 4123, with CLI log level verbose, 
# also log in file server.log with minimum level of info
npm run start -- -p 4123 -l verbose --file-log-level info --log-file server.log -c $COUCHDB_HOST --server-url $SERVER_URL
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




