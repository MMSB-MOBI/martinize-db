FROM node:13

RUN apt-get update && apt-get install -y unzip screen

COPY deploy.zip /
COPY docker_start_server.sh /

WORKDIR /

RUN unzip /deploy.zip && \
  mv deploy server && \
  cd server && \
  npm i

