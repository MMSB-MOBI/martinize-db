#!/bin/bash

# Create the deploy zip
./deploy.sh

# Remove possible existing image
docker rmi alkihis/martinize-db:1.0 2> /dev/null
# Create docker image
docker build -t alkihis/martinize-db:1.0 .

# Run the docker container
docker run --name martinize_db -p 5000:4123 -d alkihis/martinize-db:1.0

echo "Container is ready to be used. Server default port 4123 is redirected to host port 5000."

echo ""

echo "Initialize the server inside the container with \"docker exec -it martinize_db bash docker_first_run.sh\"."
echo "To restart container later, use \"docker exec -it martinize_db bash docker_start_server.sh\"."
echo "To stop, remove container and remove its image, run \"docker stop martinize_db && docker container rm martinize_db && docker rmi alkihis/martinize-db:1.0\"."
