# Martinize Database Server

> Provides JWT-authentification-based REST API for MArtinize Database, Molecule Coarse-Grainer and Membrane Builder.
 
This program serves the MArtinize Database server (**MAD**).

Requirements:
- A running Couch Database server 
- Free disk space and permissions to write on server directory in order to store sended files
- Node 12+ and NPM 6+
- A working SMTP server to send e-mails

## Deploy

An advanced deploy procedure has been written to [Deploy.md](./Deploy.md) file. Please ignore this file and see it instead.

## Setting up the project

The following procedure assume GROMACS is installed,
all your virtual environnements are set properly or 
your python executables are aware of martinize, insane and other requirements.

First, clone the project with `git`.
```bash
git clone https://github.com/alkihis/martinize-db.git
cd martinize-db
```

### Create JWT keys

This server use JSON Web Tokens, which requires to sign components with RSA keys. If you've cloned the project, you
must create the RSA keys with the section below.

```bash
mkdir .keys
cd .keys

# do not enter passphrase, press enter
ssh-keygen -t rsa -b 4096 -m PEM -f key_new.pem
openssl rsa -in key_new.pem -pubout -outform PEM -out key_new
```

### Get started

Install all the required dependencies with `npm`.
```bash
npm i
```

### Install TypeScript

Ensure you have TypeScript installed.
```bash
tsc -v # 3.8+
```

You can also use project local TypeScript compiler.
```bash
npm run tsc -- -v # 3.8+
```

**Make sure you have rights to create files and folders in project directory, otherwise the server can't work.**

### Start developement server

```bash
tsc  # Compile: will create the build/ folder
npm run start-dev
```

If everything goes well, it should be configured properly.

### Server parameters

The server obey to some parameters in order to run properly.
Take a look to `src/constants.ts` to know more.

You also have CLI args to give to server.

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


## First start

At server first start, you must set up the environnement in order to make user system and database work.

The console will remind you to do it:
```
Welcome to Martinize server CLI. For help, type "help".

WARN: The database seems to be un-initialized. Please create all databases by entering "database create all".
WARN: Once database is created, you can create an administrator account with "user create".
```

First, create the databases inside the CLI by typing: `database create all`.

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


## Deploy

### Compile and prepare deployment

You must have the `martinize-db-client` project on your computer.

Please see in [Deploy.md](./Deploy.md) part `Include front-end website to this server` to compile and type properly martinize db client.

#### In this project

***This project must have been built at least one time with `tsc`. Make sure `build` folder exists.***

Compile the required files in a directory, and zip them.

```bash

# Copy the required files to a tmp folder
mkdir deploy
cp -R .keys build static templates package.json settings.json deploy/

# Zip the folder
zip -r deploy.zip deploy

# Remove the temporary folder
rm -rf deploy
```

You now have the required file `deploy.zip` in order to deploy with Docker! Follow instuctions of **`Deploy with Docker`** part.


### Deploy with Docker

> Now, you must rename `_old-docker-compose.yml` to `docker-compose.yml`.

> In order to deploy, `docker` and `docker-compose` must be installed and configured properly.
> Docker is not configured/installed? See [Configure_Docker_IBCP.md](./Configure_Docker_IBCP.md) file.

#### Create the containers

`docker-compose` will do the job, following rules defined in the `docker-compose.yml` file.

```bash
docker-compose up -d --build
```

Containers are ready to be used after `docker-compose`. Martinize Database server default port 4123 will be redirected to host port 5000.

#### Initialize the server

Two containers are started: `martinize_couch_db` and `martinize_db`.

In `martinize_db`, server is not started automatically, we must initialize the server inside the container.

TODO: inside the container, must install Python3, Python2, insane+martinize2, compile GROMACS.

Enter inside the container, and start the shell script who automatically set a screen with the server.
```bash
docker exec -it martinize_db bash  # Open a interactive bash shell in container
source docker_start_server.sh  # Run the script
```

#### Stop or remove the deployment

To stop/remove (f.e., to re-deploy it), remove container and remove its image.

```bash
# Stop martinize_db then delete it.
docker stop martinize_db && docker container rm martinize_db

# Do the same with its couchdb counterpart
docker stop martinize_couch_db && docker container rm martinize_couch_db
```



