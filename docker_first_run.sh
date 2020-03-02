#!/bin/bash

cd server
npm run init

cd ..
./docker_start_server.sh
