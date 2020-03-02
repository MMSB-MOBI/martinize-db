#!/bin/bash

# Create the deploy zip
martinize_db_client_dir="../martinize-db-client/build"

rm -rf ./static/www/*
cp -R $martinize_db_client_dir/* ./static/www

mkdir deploy

cp -R .keys build static templates package.json settings.json deploy/

zip -r deploy.zip deploy

rm -rf deploy

# Create docker image
docker-compose up -d --build

rm deploy.zip

# docker start martinize_db

echo "Container is ready to be used. Server default port 4123 is redirected to host port 5000."

echo ""

echo "Initialize the server inside the container with \"docker exec -it martinize_db bash\", and type \"source docker_start_server.sh\" !"
echo "To stop, remove container and remove its image, run \"docker stop martinize_db && docker container rm martinize_db\"."
