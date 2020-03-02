FROM couchdb:2.3

RUN curl -sL https://deb.nodesource.com/setup_13.x | bash -

RUN apt-get update && apt-get install -y unzip screen nodejs

COPY deploy.zip /
COPY docker_start_server.sh /

WORKDIR /

RUN unzip /deploy.zip && \
  mv deploy server && \
  cd server && \
  npm i

