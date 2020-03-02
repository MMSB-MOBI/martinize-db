#!/bin/bash

# Create the deploy zip
./deploy.sh

# Create docker image
docker-compose up -d --build

# docker start martinize_db

echo "Container is ready to be used. Server default port 4123 is redirected to host port 5000."

echo ""

echo "Initialize the server inside the container with \"docker exec -it martinize_db bash\", and type \"source docker_start_server.sh\" !"
echo "To stop, remove container and remove its image, run \"docker stop martinize_db && docker container rm martinize_db\"."
